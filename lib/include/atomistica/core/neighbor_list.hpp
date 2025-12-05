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

#pragma once

#include <array>
#include <cstddef>
#include <vector>

#include "../config.hpp"

namespace atomistica {

class AtomicSystem;

/**
 * @brief Information about a neighboring atom
 */
struct Neighbor {
    std::size_t index;              // Index of neighbor atom
    std::array<int, 3> cell_shift;  // Periodic image shift (in cell vector units)
};

/**
 * @brief Cell-list based neighbor list
 *
 * Efficiently finds all pairs of atoms within a cutoff distance.
 * Uses spatial binning (cell lists) for O(N) scaling.
 *
 * This replaces the Fortran neighbors_t type.
 */
class NeighborList {
public:
    NeighborList() = default;

    // Set the interaction cutoff
    void set_cutoff(Scalar cutoff);
    Scalar cutoff() const { return cutoff_; }

    // Set Verlet shell for delayed rebuilding
    void set_verlet_shell(Scalar shell);
    Scalar verlet_shell() const { return verlet_shell_; }

    // Update neighbor list (rebuilds if needed)
    void update(const AtomicSystem& system);

    // Force rebuild on next update
    void invalidate();

    // Number of atoms in the list
    std::size_t num_atoms() const { return seed_.empty() ? 0 : seed_.size() - 1; }

    // Get neighbors of atom i (returns pair of iterators)
    using NeighborIterator = std::vector<Neighbor>::const_iterator;
    std::pair<NeighborIterator, NeighborIterator> neighbors(std::size_t i) const;

    // Number of neighbors for atom i
    std::size_t num_neighbors(std::size_t i) const;

    // Total number of pairs
    std::size_t num_pairs() const { return neighbors_.size(); }

    // Check if update is needed
    bool needs_update(const AtomicSystem& system) const;

private:
    // Build the neighbor list from scratch
    void build(const AtomicSystem& system);

    // Cell list data structures
    void build_cell_list(const AtomicSystem& system);

    Scalar cutoff_ = 0.0;
    Scalar verlet_shell_ = 0.0;

    // Cached system state for checking if rebuild needed
    int cached_position_revision_ = -1;
    int cached_cell_revision_ = -1;

    // Per-atom neighbor storage
    // neighbors_[seed_[i]] to neighbors_[seed_[i+1]-1] are neighbors of atom i
    std::vector<std::size_t> seed_;  // Start index for each atom (size: num_atoms + 1)
    std::vector<Neighbor> neighbors_;

    // Cell list internals
    std::array<int, 3> num_cells_;
    std::vector<int> cell_list_;        // Linked list: cell_list_[i] = next atom in cell, -1 if end
    std::vector<int> cell_head_;        // Head of linked list for each cell
    Vec3 cell_size_;
    Mat3 inverse_cell_;

    // Saved positions for Verlet shell checking
    Array3X saved_positions_;
};

} // namespace atomistica
