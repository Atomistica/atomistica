// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
//
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ======================================================================

#include "atomistica/core/neighbor_list.hpp"
#include "atomistica/core/atomic_system.hpp"

#include <algorithm>
#include <cmath>

namespace atomistica {

void NeighborList::set_cutoff(Scalar cutoff) {
    if (cutoff != cutoff_) {
        cutoff_ = cutoff;
        invalidate();
    }
}

void NeighborList::set_verlet_shell(Scalar shell) {
    verlet_shell_ = shell;
}

void NeighborList::invalidate() {
    cached_position_revision_ = -1;
    cached_cell_revision_ = -1;
}

bool NeighborList::needs_update(const AtomicSystem& system) const {
    // Always rebuild if revisions don't match
    if (cached_position_revision_ != system.position_revision() ||
        cached_cell_revision_ != system.cell_revision()) {
        return true;
    }

    // Check if any atom moved more than half the Verlet shell
    if (verlet_shell_ > 0 && saved_positions_.cols() == static_cast<Eigen::Index>(system.num_atoms())) {
        Scalar max_displacement_sq = 0.0;
        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            Vec3 dr = system.position(i).matrix() - saved_positions_.col(i).matrix();
            max_displacement_sq = std::max(max_displacement_sq, dr.squaredNorm());
        }
        Scalar threshold = verlet_shell_ / 2.0;
        if (max_displacement_sq > threshold * threshold) {
            return true;
        }
    }

    return false;
}

void NeighborList::update(const AtomicSystem& system) {
    if (needs_update(system)) {
        build(system);
        cached_position_revision_ = system.position_revision();
        cached_cell_revision_ = system.cell_revision();

        // Save positions for Verlet shell checking
        if (verlet_shell_ > 0) {
            saved_positions_ = system.positions();
        }
    }
}

std::pair<NeighborList::NeighborIterator, NeighborList::NeighborIterator>
NeighborList::neighbors(std::size_t i) const {
    if (i >= seed_.size() - 1) {
        return {neighbors_.end(), neighbors_.end()};
    }
    return {neighbors_.begin() + seed_[i], neighbors_.begin() + seed_[i + 1]};
}

std::size_t NeighborList::num_neighbors(std::size_t i) const {
    if (i >= seed_.size() - 1) {
        return 0;
    }
    return seed_[i + 1] - seed_[i];
}

void NeighborList::build_cell_list(const AtomicSystem& system) {
    const Scalar total_cutoff = cutoff_ + verlet_shell_;

    // Get cell parameters
    const Mat3& cell = system.cell();
    inverse_cell_ = cell.inverse();

    // Determine number of cells in each direction
    // Cell size should be at least the cutoff
    for (int d = 0; d < 3; ++d) {
        Vec3 axis = cell.col(d);
        Scalar length = axis.norm();

        if (system.pbc()[d]) {
            num_cells_[d] = std::max(1, static_cast<int>(std::floor(length / total_cutoff)));
        } else {
            // For non-periodic, use single cell (or could use multiple for efficiency)
            num_cells_[d] = std::max(1, static_cast<int>(std::floor(length / total_cutoff)));
        }
        cell_size_(d) = length / num_cells_[d];
    }

    // Initialize cell list
    int total_cells = num_cells_[0] * num_cells_[1] * num_cells_[2];
    cell_head_.assign(total_cells, -1);
    cell_list_.resize(system.num_atoms());

    // Assign atoms to cells
    for (std::size_t i = 0; i < system.num_atoms(); ++i) {
        // Convert to fractional coordinates
        Vec3 s = inverse_cell_ * system.position(i).matrix();

        // Wrap into [0, 1) for periodic directions
        for (int d = 0; d < 3; ++d) {
            if (system.pbc()[d]) {
                s(d) -= std::floor(s(d));
            } else {
                s(d) = std::clamp(s(d), 0.0, 1.0 - 1e-10);
            }
        }

        // Compute cell indices
        int ix = std::clamp(static_cast<int>(s(0) * num_cells_[0]), 0, num_cells_[0] - 1);
        int iy = std::clamp(static_cast<int>(s(1) * num_cells_[1]), 0, num_cells_[1] - 1);
        int iz = std::clamp(static_cast<int>(s(2) * num_cells_[2]), 0, num_cells_[2] - 1);

        int cell_idx = ix + num_cells_[0] * (iy + num_cells_[1] * iz);

        // Insert at head of linked list
        cell_list_[i] = cell_head_[cell_idx];
        cell_head_[cell_idx] = static_cast<int>(i);
    }
}

void NeighborList::build(const AtomicSystem& system) {
    const std::size_t num_atoms = system.num_atoms();
    const Scalar total_cutoff = cutoff_ + verlet_shell_;
    const Scalar cutoff_sq = total_cutoff * total_cutoff;

    // Build cell list
    build_cell_list(system);

    // Clear and prepare storage
    seed_.resize(num_atoms + 1);
    neighbors_.clear();
    neighbors_.reserve(num_atoms * 20);  // Estimate ~20 neighbors per atom

    const Mat3& cell = system.cell();
    const auto& pbc = system.pbc();

    // For each atom, find neighbors
    for (std::size_t i = 0; i < num_atoms; ++i) {
        seed_[i] = neighbors_.size();

        Vec3 ri = system.position(i).matrix();

        // Convert to fractional coordinates
        Vec3 si = inverse_cell_ * ri;

        // Wrap to [0, 1) and determine cell of atom i
        for (int d = 0; d < 3; ++d) {
            if (pbc[d]) {
                si(d) -= std::floor(si(d));
            } else {
                si(d) = std::clamp(si(d), 0.0, 1.0 - 1e-10);
            }
        }

        int cix = std::clamp(static_cast<int>(si(0) * num_cells_[0]), 0, num_cells_[0] - 1);
        int ciy = std::clamp(static_cast<int>(si(1) * num_cells_[1]), 0, num_cells_[1] - 1);
        int ciz = std::clamp(static_cast<int>(si(2) * num_cells_[2]), 0, num_cells_[2] - 1);

        // Search neighboring cells
        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    // Cell indices with wrapping
                    int jx = cix + dx;
                    int jy = ciy + dy;
                    int jz = ciz + dz;

                    // Cell shift for periodic images
                    std::array<int, 3> shift = {0, 0, 0};

                    // Handle periodic boundaries
                    if (pbc[0]) {
                        if (jx < 0) { jx += num_cells_[0]; shift[0] = -1; }
                        else if (jx >= num_cells_[0]) { jx -= num_cells_[0]; shift[0] = 1; }
                    } else {
                        if (jx < 0 || jx >= num_cells_[0]) continue;
                    }

                    if (pbc[1]) {
                        if (jy < 0) { jy += num_cells_[1]; shift[1] = -1; }
                        else if (jy >= num_cells_[1]) { jy -= num_cells_[1]; shift[1] = 1; }
                    } else {
                        if (jy < 0 || jy >= num_cells_[1]) continue;
                    }

                    if (pbc[2]) {
                        if (jz < 0) { jz += num_cells_[2]; shift[2] = -1; }
                        else if (jz >= num_cells_[2]) { jz -= num_cells_[2]; shift[2] = 1; }
                    } else {
                        if (jz < 0 || jz >= num_cells_[2]) continue;
                    }

                    int neighbor_cell = jx + num_cells_[0] * (jy + num_cells_[1] * jz);

                    // Iterate through atoms in neighboring cell
                    int j = cell_head_[neighbor_cell];
                    while (j >= 0) {
                        std::size_t ju = static_cast<std::size_t>(j);

                        // Skip self (full neighbor list: store both i->j and j->i)
                        if (ju != i) {
                            // Compute distance with periodic shift
                            Vec3 rj = system.position(ju).matrix();
                            Vec3 dr = rj - ri;

                            // Add periodic shift
                            dr += cell.col(0) * shift[0];
                            dr += cell.col(1) * shift[1];
                            dr += cell.col(2) * shift[2];

                            Scalar dist_sq = dr.squaredNorm();

                            if (dist_sq < cutoff_sq) {
                                neighbors_.push_back(Neighbor{ju, shift});
                            }
                        }

                        j = cell_list_[j];
                    }
                }
            }
        }
    }

    seed_[num_atoms] = neighbors_.size();
}

} // namespace atomistica
