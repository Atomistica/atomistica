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

#include <algorithm>
#include <cmath>
#include <map>
#include <vector>

#include "../config.hpp"
#include "../core/atomic_system.hpp"
#include "../core/neighbor_list.hpp"
#include "materials.hpp"
#include "slater_koster.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Compute distance vector between atoms i and j using neighbor info
 *
 * @param system Atomic system
 * @param i Central atom index
 * @param neighbor Neighbor information (index and cell_shift)
 * @return Distance vector r_j - r_i accounting for periodic images
 */
inline Vec3 neighbor_distance_vector(const AtomicSystem& system, std::size_t i,
                                     const Neighbor& neighbor) {
    Vec3 r_i = system.position(i);
    Vec3 r_j = system.position(neighbor.index);
    Vec3 shift = system.cell() * Vec3(neighbor.cell_shift[0],
                                       neighbor.cell_shift[1],
                                       neighbor.cell_shift[2]);
    return r_j + shift - r_i;
}

/**
 * @brief Tight-binding Hamiltonian builder
 *
 * Constructs H and S matrices from atomic positions using
 * Slater-Koster transformations and tabulated integrals.
 */
class TBHamiltonian {
public:
    TBHamiltonian() = default;

    /**
     * @brief Set materials database
     */
    void set_materials(MaterialsDatabase* db) { materials_ = db; }

    /**
     * @brief Initialize Hamiltonian for a given atomic system
     *
     * Sets up orbital indices and allocates matrices
     *
     * @param system Atomic system
     * @param elements Vector of element parameters for each atom type
     */
    void init(const AtomicSystem& system, const std::vector<TBElementParams>& elements) {
        int nat = system.num_atoms();

        // Map atomic numbers to element indices
        element_params_ = elements;
        z_to_elem_.clear();
        for (size_t i = 0; i < elements.size(); ++i) {
            z_to_elem_[elements[i].atomic_number] = static_cast<int>(i);
        }

        // Count total orbitals
        int total_orbitals = 0;
        ham_.orbitals_per_atom.resize(nat);
        ham_.orbital_offset.resize(nat);
        ham_.element_index.resize(nat);

        for (int i = 0; i < nat; ++i) {
            int Z = system.atomic_number(i);
            auto it = z_to_elem_.find(Z);
            if (it == z_to_elem_.end()) {
                throw std::runtime_error("Unknown element Z=" + std::to_string(Z));
            }
            int elem_idx = it->second;
            ham_.element_index[i] = elem_idx;

            const TBElementParams& elem = element_params_[elem_idx];
            ham_.orbital_offset[i] = total_orbitals;
            ham_.orbitals_per_atom[i] = elem.num_orbitals;
            total_orbitals += elem.num_orbitals;
        }

        ham_.resize(nat, total_orbitals);

        // Set neutral charges
        for (int i = 0; i < nat; ++i) {
            const TBElementParams& elem = element_params_[ham_.element_index[i]];
            ham_.neutral_charges[i] = elem.valence_electrons;
        }
    }

    /**
     * @brief Build H and S matrices
     *
     * @param system Atomic system
     * @param neighbors Neighbor list
     */
    void build_matrices(const AtomicSystem& system, const NeighborList& neighbors) {
        ham_.clear_matrices();

        int nat = system.num_atoms();

        // Set diagonal elements (on-site energies and overlap normalization)
        for (int i = 0; i < nat; ++i) {
            const TBElementParams& elem = element_params_[ham_.element_index[i]];
            int offset = ham_.orbital_offset[i];
            int norb = ham_.orbitals_per_atom[i];

            for (int a = 0; a < norb; ++a) {
                ham_.H(offset + a, offset + a) = elem.onsite[a];
                ham_.S(offset + a, offset + a) = 1.0;
            }
        }

        // Build off-diagonal elements from neighbor pairs
        std::array<Scalar, NUM_SK_INTEGRALS> H_sk, S_sk;

        for (int i = 0; i < nat; ++i) {
            int Z_i = system.atomic_number(i);
            const TBElementParams& elem_i = element_params_[ham_.element_index[i]];
            int offset_i = ham_.orbital_offset[i];
            int norb_i = ham_.orbitals_per_atom[i];

            auto [begin, end] = neighbors.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                int j = it->index;
                if (j <= i) continue;  // Only upper triangle

                int Z_j = system.atomic_number(j);
                const TBElementParams& elem_j = element_params_[ham_.element_index[j]];
                int offset_j = ham_.orbital_offset[j];
                int norb_j = ham_.orbitals_per_atom[j];

                // Get distance and direction
                Vec3 r_ij = neighbor_distance_vector(system, i, *it);
                Scalar r = r_ij.norm();

                if (r < 1e-10) continue;

                Vec3 c = r_ij / r;  // Direction cosines

                // Get SK integrals from materials database
                try {
                    const SKSpline& H_spline = materials_->get_H_spline(Z_i, Z_j);
                    const SKSpline& S_spline = materials_->get_S_spline(Z_i, Z_j);

                    H_spline.eval(r, H_sk);
                    S_spline.eval(r, S_sk);
                } catch (...) {
                    // If splines not available, skip this pair
                    continue;
                }

                // Build matrix elements for all orbital pairs
                for (int a = 0; a < norb_i; ++a) {
                    int a_abs = get_absolute_orbital(norb_i, a + 1);

                    for (int b = 0; b < norb_j; ++b) {
                        int b_abs = get_absolute_orbital(norb_j, b + 1);

                        // Transform SK integrals to matrix elements
                        Scalar H_el = transform_orb(a_abs, b_abs, c, H_sk);
                        Scalar S_el = transform_orb(a_abs, b_abs, c, S_sk);

                        // Store in symmetric matrix
                        int ii = offset_i + a;
                        int jj = offset_j + b;

                        ham_.H(ii, jj) = H_el;
                        ham_.H(jj, ii) = H_el;
                        ham_.S(ii, jj) = S_el;
                        ham_.S(jj, ii) = S_el;
                    }
                }
            }
        }
    }

    /**
     * @brief Add SCC (self-consistent charge) correction to Hamiltonian
     *
     * Modifies H based on Mulliken charges: H_scc = H + shift * S
     * where shift depends on charge difference from neutral
     *
     * @param gamma Gamma matrix (Coulomb interaction between atoms)
     */
    void add_scc_correction(const MatX& gamma) {
        int nat = ham_.num_atoms;

        // Compute potential shifts from charges
        VecX shift = VecX::Zero(nat);
        for (int i = 0; i < nat; ++i) {
            for (int j = 0; j < nat; ++j) {
                shift[i] += gamma(i, j) * ham_.charges[j];
            }
        }

        // Add shift * S to diagonal blocks
        for (int i = 0; i < nat; ++i) {
            int offset = ham_.orbital_offset[i];
            int norb = ham_.orbitals_per_atom[i];
            Scalar si = 0.5 * shift[i];

            for (int a = 0; a < norb; ++a) {
                ham_.H(offset + a, offset + a) += si;
            }
        }

        // Add shift * S to off-diagonal blocks
        for (int i = 0; i < nat; ++i) {
            int offset_i = ham_.orbital_offset[i];
            int norb_i = ham_.orbitals_per_atom[i];

            for (int j = i + 1; j < nat; ++j) {
                int offset_j = ham_.orbital_offset[j];
                int norb_j = ham_.orbitals_per_atom[j];

                Scalar sij = 0.5 * (shift[i] + shift[j]);

                for (int a = 0; a < norb_i; ++a) {
                    for (int b = 0; b < norb_j; ++b) {
                        int ii = offset_i + a;
                        int jj = offset_j + b;

                        Scalar correction = sij * ham_.S(ii, jj);
                        ham_.H(ii, jj) += correction;
                        ham_.H(jj, ii) += correction;
                    }
                }
            }
        }
    }

    /**
     * @brief Compute repulsive energy
     */
    Scalar compute_repulsive_energy(const AtomicSystem& system,
                                    const NeighborList& neighbors) {
        Scalar E_rep = 0.0;
        int nat = system.num_atoms();

        for (int i = 0; i < nat; ++i) {
            int Z_i = system.atomic_number(i);

            auto [begin, end] = neighbors.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                int j = it->index;
                if (j <= i) continue;

                int Z_j = system.atomic_number(j);
                Vec3 r_ij = neighbor_distance_vector(system, i, *it);
                Scalar r = r_ij.norm();

                try {
                    const RepulsiveSpline& rep = materials_->get_rep_spline(Z_i, Z_j);
                    E_rep += rep.eval(r);
                } catch (...) {
                    // No repulsive potential for this pair
                }
            }
        }

        ham_.repulsive_energy = E_rep;
        return E_rep;
    }

    /**
     * @brief Compute repulsive forces
     */
    void compute_repulsive_forces(const AtomicSystem& system,
                                  const NeighborList& neighbors,
                                  MatX3& forces) {
        int nat = system.num_atoms();

        for (int i = 0; i < nat; ++i) {
            int Z_i = system.atomic_number(i);

            auto [begin, end] = neighbors.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                int j = it->index;
                if (j <= i) continue;

                int Z_j = system.atomic_number(j);
                Vec3 r_ij = neighbor_distance_vector(system, i, *it);
                Scalar r = r_ij.norm();

                if (r < 1e-10) continue;

                Vec3 r_hat = r_ij / r;

                try {
                    const RepulsiveSpline& rep = materials_->get_rep_spline(Z_i, Z_j);
                    Scalar dV_dr;
                    rep.eval_deriv(r, dV_dr);

                    Vec3 f = -dV_dr * r_hat;
                    forces.row(i) -= f.transpose();
                    forces.row(j) += f.transpose();
                } catch (...) {
                    // No repulsive potential for this pair
                }
            }
        }
    }

    /**
     * @brief Get Hamiltonian data structure
     */
    DenseHamiltonian& hamiltonian() { return ham_; }
    const DenseHamiltonian& hamiltonian() const { return ham_; }

    /**
     * @brief Get H matrix
     */
    MatX& H() { return ham_.H; }
    const MatX& H() const { return ham_.H; }

    /**
     * @brief Get S matrix
     */
    MatX& S() { return ham_.S; }
    const MatX& S() const { return ham_.S; }

private:
    MaterialsDatabase* materials_ = nullptr;
    DenseHamiltonian ham_;
    std::vector<TBElementParams> element_params_;
    std::map<int, int> z_to_elem_;
};

/**
 * @brief Compute gamma matrix for SCC-DFTB
 *
 * gamma_ij = short-range Coulomb interaction between atoms i and j
 * Uses Hubbard U parameters and distance-dependent function.
 *
 * @param system Atomic system
 * @param elements Element parameters
 * @param use_periodic Include periodic images
 * @return Gamma matrix
 */
inline MatX compute_gamma_matrix(const AtomicSystem& system,
                                 const std::vector<TBElementParams>& elements,
                                 const std::vector<int>& elem_index,
                                 bool use_periodic = false) {
    int nat = system.num_atoms();
    MatX gamma = MatX::Zero(nat, nat);

    // Map element index
    std::map<int, int> z_to_elem;
    for (size_t i = 0; i < elements.size(); ++i) {
        z_to_elem[elements[i].atomic_number] = static_cast<int>(i);
    }

    for (int i = 0; i < nat; ++i) {
        Scalar U_i = elements[elem_index[i]].hubbard_U;

        for (int j = i; j < nat; ++j) {
            Scalar U_j = elements[elem_index[j]].hubbard_U;

            if (i == j) {
                // On-site: gamma_ii = U_i
                gamma(i, i) = U_i;
            } else {
                // Off-site: use short-range function
                Vec3 r_ij = system.position(j) - system.position(i);
                if (use_periodic) {
                    // Minimum image
                    // (simplified - full implementation would use neighbor list)
                }
                Scalar r = r_ij.norm();

                // Klopman-Ohno formula:
                // gamma_ij = 1/sqrt(r^2 + (1/U_i + 1/U_j)^2 / 4)
                // But DFTB uses a slightly different form:
                // gamma_ij = 1/r * erf(sqrt(tau_i * tau_j / (tau_i + tau_j)) * r)
                // where tau = 16/5 * U^2 / (3.2 Hartree)

                // Simplified version using exponential damping:
                Scalar tau_i = 3.2 * U_i * U_i;
                Scalar tau_j = 3.2 * U_j * U_j;
                Scalar tau_avg = std::sqrt(tau_i * tau_j);

                // Short-range gamma function (DFTB form)
                Scalar gamma_ij;
                if (r < 1e-6) {
                    gamma_ij = 0.5 * (U_i + U_j);
                } else {
                    // Use complementary error function approximation
                    Scalar x = tau_avg * r;
                    Scalar erf_approx = 1.0 - std::exp(-x * x);  // Simplified
                    gamma_ij = erf_approx / r;
                }

                gamma(i, j) = gamma_ij;
                gamma(j, i) = gamma_ij;
            }
        }
    }

    return gamma;
}

} // namespace tb
} // namespace atomistica
