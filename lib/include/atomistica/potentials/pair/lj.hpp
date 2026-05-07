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

#include <cmath>
#include <map>
#include <tuple>
#include <utility>

#include "../potential_base.hpp"
#include "../../config.hpp"

namespace atomistica {

/**
 * @brief Lennard-Jones 12-6 pair potential
 *
 * V(r) = 4 * epsilon * [(sigma/r)^12 - (sigma/r)^6]
 *
 * with optional energy shift at cutoff.
 *
 * @tparam Shift If true, shift potential to zero at cutoff
 */
template<bool Shift = false>
class LJPotential : public PotentialBase<LJPotential<Shift>> {
public:
    /**
     * @brief Parameters for a pair of element types
     */
    struct PairParams {
        Scalar epsilon;  // Well depth
        Scalar sigma;    // Zero-crossing distance
        Scalar cutoff;   // Interaction cutoff

        // Precomputed values for efficiency
        Scalar sigma6;
        Scalar sigma12;
        Scalar shift_energy;  // Energy at cutoff (for shifting)
    };

    LJPotential() = default;

    /**
     * @brief Construct with single element parameters
     */
    LJPotential(int Z, Scalar epsilon, Scalar sigma, Scalar cutoff) {
        set_params(Z, Z, epsilon, sigma, cutoff);
    }

    /**
     * @brief Set parameters for a pair of elements
     *
     * @param Z1 First element (atomic number)
     * @param Z2 Second element (atomic number)
     * @param epsilon Well depth
     * @param sigma Zero-crossing distance
     * @param cutoff Interaction cutoff
     */
    void set_params(int Z1, int Z2, Scalar epsilon, Scalar sigma, Scalar cutoff) {
        PairParams params;
        params.epsilon = epsilon;
        params.sigma = sigma;
        params.cutoff = cutoff;
        params.sigma6 = std::pow(sigma, 6);
        params.sigma12 = params.sigma6 * params.sigma6;

        // Compute shift energy
        if constexpr (Shift) {
            Scalar r6 = std::pow(cutoff, 6);
            Scalar r12 = r6 * r6;
            params.shift_energy = 4.0 * epsilon * (params.sigma12 / r12 - params.sigma6 / r6);
        } else {
            params.shift_energy = 0.0;
        }

        // Store for both orderings
        auto key1 = std::make_pair(std::min(Z1, Z2), std::max(Z1, Z2));
        params_[key1] = params;

        // Update max cutoff
        if (cutoff > max_cutoff_) {
            max_cutoff_ = cutoff;
        }
    }

    /**
     * @brief Get parameters for a pair
     */
    const PairParams* get_params(int Z1, int Z2) const {
        auto key = std::make_pair(std::min(Z1, Z2), std::max(Z1, Z2));
        auto it = params_.find(key);
        if (it != params_.end()) {
            return &it->second;
        }
        return nullptr;
    }

    // CRTP implementation
    Scalar cutoff_impl() const {
        return max_cutoff_;
    }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial) {
        PotentialResults results;
        const std::size_t num_atoms = system.num_atoms();
        const Mat3& cell = system.cell();

        for (std::size_t i = 0; i < num_atoms; ++i) {
            int Zi = system.atomic_numbers()(i);
            Vec3 ri = system.position(i).matrix();

            auto [begin, end] = neighbors.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                const auto& neigh = *it;
                std::size_t j = neigh.index;
                int Zj = system.atomic_numbers()(j);

                // Get parameters for this pair
                const PairParams* params = get_params(Zi, Zj);
                if (!params) continue;

                // Compute distance vector with periodic shift
                Vec3 rj = system.position(j).matrix();
                Vec3 dr = rj - ri;
                dr += cell.col(0) * neigh.cell_shift[0];
                dr += cell.col(1) * neigh.cell_shift[1];
                dr += cell.col(2) * neigh.cell_shift[2];

                Scalar r_sq = dr.squaredNorm();
                Scalar cutoff_sq = params->cutoff * params->cutoff;

                if (r_sq >= cutoff_sq) continue;

                Scalar r6 = r_sq * r_sq * r_sq;
                Scalar r12 = r6 * r6;

                // LJ energy: 4*eps*[(sigma/r)^12 - (sigma/r)^6]
                Scalar energy = 4.0 * params->epsilon *
                    (params->sigma12 / r12 - params->sigma6 / r6);

                if constexpr (Shift) {
                    energy -= params->shift_energy;
                }

                // Full neighbor list: each pair counted twice, so halve energy
                results.energy += 0.5 * energy;

                if (compute_forces || compute_virial) {
                    // Force on i due to j: F_i = -dV/dr * (ri - rj)/|ri - rj| = dV/dr * dr/r
                    // where dr = rj - ri
                    // dV/dr = 4*eps*[-12*sigma^12/r^13 + 6*sigma^6/r^7]
                    //       = 24*eps/r * [sigma^6/r^6 - 2*sigma^12/r^12]
                    // F_i = -F_ij = 24*eps/r^2 * [2*sigma^12/r^12 - sigma^6/r^6] * dr
                    Scalar force_over_r = 24.0 * params->epsilon / r_sq *
                        (2.0 * params->sigma12 / r12 - params->sigma6 / r6);

                    Vec3 force = force_over_r * dr;

                    if (compute_forces) {
                        // Full neighbor list: only add force to atom i
                        // The reverse pair j->i will add force to j
                        system.forces().col(i) -= force.array();
                    }

                    if (compute_virial) {
                        // Virial: W_ab = sum_{i,j} r_ij,a * f_i,b / 2
                        // Full neighbor list: halve contribution
                        results.virial -= 0.5 * dr * force.transpose();
                    }
                }
            }
        }

        return results;
    }

private:
    std::map<std::pair<int, int>, PairParams> params_;
    Scalar max_cutoff_ = 0.0;
};

// Common type aliases
using LJCut = LJPotential<false>;
using LJCutShift = LJPotential<true>;

} // namespace atomistica
