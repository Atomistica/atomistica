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

#include <atomistica/tightbinding/dftb.hpp>

#include <set>

namespace atomistica {
namespace tb {

void DFTB::init(const AtomicSystem& system) {
    // Collect unique elements
    elements_.clear();
    std::set<int> unique_Z;
    for (std::size_t i = 0; i < system.num_atoms(); ++i) {
        unique_Z.insert(system.atomic_number(i));
    }

    // Load SKF files for all required pairs (diagonal files contain element params)
    for (int Z1 : unique_Z) {
        for (int Z2 : unique_Z) {
            if (Z2 >= Z1) {
                materials_.load_pair(Z1, Z2);
            }
        }
    }

    // Now collect element parameters
    for (int Z : unique_Z) {
        if (materials_.has_element(Z)) {
            elements_.push_back(materials_.get_element(Z));
        }
    }

    // Initialize Hamiltonian
    hamiltonian_.init(system, elements_);

    // Count electrons
    n_electrons_ = 0.0;
    for (std::size_t i = 0; i < system.num_atoms(); ++i) {
        int Z = system.atomic_number(i);
        n_electrons_ += materials_.get_element(Z).valence_electrons;
    }

    // Reset mixer for new system
    mixer_.reset();

    initialized_ = true;
}

Scalar DFTB::compute_energy(const AtomicSystem& system,
                            const NeighborList& neighbors) {
    if (!initialized_) init(system);

    // Build H and S matrices
    hamiltonian_.build_matrices(system, neighbors);

    if (enable_scc_) {
        // SCC iteration
        run_scc(system, neighbors);
    } else {
        // Single diagonalization
        solve_electronic();
    }

    // Compute total energy
    DenseHamiltonian& ham = hamiltonian_.hamiltonian();
    ham.band_energy = solver_.compute_band_energy(ham);

    // Add repulsive energy
    Scalar E_rep = hamiltonian_.compute_repulsive_energy(system, neighbors);

    total_energy_ = ham.band_energy + E_rep;

    if (enable_scc_) {
        // Add SCC energy correction
        total_energy_ += compute_scc_energy();
    }

    return total_energy_;
}

Scalar DFTB::compute(const AtomicSystem& system, const NeighborList& neighbors,
                     MatX3& forces) {
    // First compute energy
    Scalar energy = compute_energy(system, neighbors);

    // Initialize forces
    int nat = system.num_atoms();
    forces = MatX3::Zero(nat, 3);

    // Compute band structure forces (Hellmann-Feynman)
    compute_band_forces(system, neighbors, forces);

    // Add repulsive forces
    hamiltonian_.compute_repulsive_forces(system, neighbors, forces);

    if (enable_scc_) {
        // Add SCC force corrections
        compute_scc_forces(system, neighbors, forces);
    }

    return energy;
}

Scalar DFTB::compute_with_stress(const AtomicSystem& system, const NeighborList& neighbors,
                                  MatX3& forces, Mat3& stress) {
    // First compute energy
    Scalar energy = compute_energy(system, neighbors);

    // Initialize forces and stress
    int nat = system.num_atoms();
    forces = MatX3::Zero(nat, 3);
    stress = Mat3::Zero();

    // Compute band structure forces and stress
    compute_band_forces_and_stress(system, neighbors, forces, stress);

    // Add repulsive forces and stress
    compute_repulsive_forces_and_stress(system, neighbors, forces, stress);

    if (enable_scc_) {
        // Add SCC force and stress corrections
        compute_scc_forces_and_stress(system, neighbors, forces, stress);
    }

    // Convert stress to per-volume (divide by cell volume)
    Scalar volume = std::abs(system.cell().determinant());
    if (volume > 1e-10) {
        stress /= volume;
    }

    return energy;
}

void DFTB::solve_electronic() {
    DenseHamiltonian& ham = hamiltonian_.hamiltonian();

    // Solve generalized eigenvalue problem
    solver_.solve(ham);

    // Compute occupation numbers
    solver_.compute_occupation(ham, n_electrons_, 2);

    // Build density matrix
    solver_.build_density_matrix(ham);

    // Compute Mulliken charges
    solver_.compute_mulliken_charges(ham);

    // Build energy-weighted density for forces
    solver_.build_energy_weighted_density(ham);
}

void DFTB::compute_gamma(const AtomicSystem& system, const NeighborList& neighbors) {
    const DenseHamiltonian& ham = hamiltonian_.hamiltonian();
    int nat = ham.num_atoms;

    gamma_ = MatX::Zero(nat, nat);

    // On-site terms (diagonal)
    for (int i = 0; i < nat; ++i) {
        Scalar U_i = elements_[ham.element_index[i]].hubbard_U;
        gamma_(i, i) = U_i;
    }

    // Off-site terms from neighbor list
    for (std::size_t i = 0; i < static_cast<std::size_t>(nat); ++i) {
        Scalar U_i = elements_[ham.element_index[i]].hubbard_U;

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;
            if (j <= i) continue;

            Scalar U_j = elements_[ham.element_index[j]].hubbard_U;

            Vec3 r_ij = neighbor_distance_vector(system, i, *it);
            Scalar r = r_ij.norm();

            Scalar gamma_ij = gamma_function(r, U_i, U_j);
            gamma_(i, j) = gamma_ij;
            gamma_(j, i) = gamma_ij;
        }
    }
}

void DFTB::run_scc(const AtomicSystem& system, const NeighborList& neighbors) {
    DenseHamiltonian& ham = hamiltonian_.hamiltonian();
    int nat = ham.num_atoms;

    // Compute gamma matrix
    compute_gamma(system, neighbors);

    // Store original H matrix (without SCC correction)
    H0_ = ham.H;

    // Initialize charges from neutral (dq = 0)
    VecX charges_input = VecX::Zero(nat);

    // Reset mixer for new SCC cycle
    mixer_.reset();

    last_scc_iterations_ = 0;

    for (int iter = 0; iter < scc_params_.max_iterations; ++iter) {
        last_scc_iterations_ = iter + 1;

        // Update charges in Hamiltonian for correction
        ham.charges = charges_input;

        // Restore original H and add SCC correction
        ham.H = H0_;
        hamiltonian_.add_scc_correction(gamma_);

        // Solve electronic structure
        solve_electronic();

        // Get new charges from Mulliken analysis (these are delta charges)
        VecX charges_output = ham.charges;

        // Mix charges using Anderson mixer
        bool converged = mixer_.mix(iter + 1, charges_input, charges_output,
                                    scc_params_.mixing_parameter,
                                    scc_params_.convergence_threshold);

        if (converged) {
            break;
        }
    }

    // Final charges
    ham.charges = charges_input;
}

Scalar DFTB::compute_scc_energy() const {
    const DenseHamiltonian& ham = hamiltonian_.hamiltonian();
    int nat = ham.num_atoms;

    Scalar E_scc = 0.0;
    for (int i = 0; i < nat; ++i) {
        for (int j = 0; j < nat; ++j) {
            E_scc += 0.5 * gamma_(i, j) * ham.charges[i] * ham.charges[j];
        }
    }

    return E_scc;
}

void DFTB::compute_band_forces(const AtomicSystem& system,
                               const NeighborList& neighbors,
                               MatX3& forces) {
    Mat3 dummy_stress;
    compute_band_forces_impl(system, neighbors, forces, dummy_stress, false);
}

void DFTB::compute_band_forces_and_stress(const AtomicSystem& system,
                                          const NeighborList& neighbors,
                                          MatX3& forces, Mat3& stress) {
    compute_band_forces_impl(system, neighbors, forces, stress, true);
}

void DFTB::compute_band_forces_impl(const AtomicSystem& system,
                                    const NeighborList& neighbors,
                                    MatX3& forces, Mat3& stress,
                                    bool compute_stress) {
    const DenseHamiltonian& ham = hamiltonian_.hamiltonian();
    int nat = system.num_atoms();

    std::array<Scalar, NUM_SK_INTEGRALS> H_sk, S_sk, dH_sk, dS_sk;

    for (std::size_t i = 0; i < static_cast<std::size_t>(nat); ++i) {
        int Z_i = system.atomic_number(i);
        int offset_i = ham.orbital_offset[i];
        int norb_i = ham.orbitals_per_atom[i];

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;
            if (j <= i) continue;

            int Z_j = system.atomic_number(j);
            int offset_j = ham.orbital_offset[j];
            int norb_j = ham.orbitals_per_atom[j];

            Vec3 r_ij = neighbor_distance_vector(system, i, *it);
            Scalar r = r_ij.norm();

            if (r < 1e-10) continue;

            Vec3 r_hat = r_ij / r;

            // Get SK integrals and derivatives
            try {
                const SKSpline& H_spline = materials_.get_H_spline(Z_i, Z_j);
                const SKSpline& S_spline = materials_.get_S_spline(Z_i, Z_j);

                H_spline.eval_deriv(r, H_sk, dH_sk);
                S_spline.eval_deriv(r, S_sk, dS_sk);
            } catch (...) {
                continue;
            }

            // Compute force contribution from each orbital pair
            Vec3 force_ij = Vec3::Zero();

            for (int a = 0; a < norb_i; ++a) {
                int a_abs = get_absolute_orbital(norb_i, a + 1);
                int ii = offset_i + a;

                for (int b = 0; b < norb_j; ++b) {
                    int b_abs = get_absolute_orbital(norb_j, b + 1);
                    int jj = offset_j + b;

                    // Get matrix element derivatives
                    Vec3 dH = transform_orb_derivative(a_abs, b_abs, r_hat, r, H_sk, dH_sk);
                    Vec3 dS = transform_orb_derivative(a_abs, b_abs, r_hat, r, S_sk, dS_sk);

                    // Hellmann-Feynman contribution
                    // F = -Tr(rho * dH) + Tr(E * dS)
                    Scalar rho_ij = ham.rho(ii, jj);
                    Scalar e_ij = ham.e_matrix(ii, jj);

                    force_ij -= 2.0 * rho_ij * dH;  // Factor of 2 for symmetric matrix
                    force_ij += 2.0 * e_ij * dS;
                }
            }

            // Add to forces (Newton's 3rd law)
            forces.row(i) -= force_ij.transpose();
            forces.row(j) += force_ij.transpose();

            // Add to stress tensor (virial contribution)
            if (compute_stress) {
                // stress_ab = -sum_ij f_ij_a * r_ij_b
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        stress(a, b) -= force_ij[a] * r_ij[b];
                    }
                }
            }
        }
    }
}

void DFTB::compute_repulsive_forces_and_stress(const AtomicSystem& system,
                                                const NeighborList& neighbors,
                                                MatX3& forces, Mat3& stress) {
    int nat = system.num_atoms();

    for (std::size_t i = 0; i < static_cast<std::size_t>(nat); ++i) {
        int Z_i = system.atomic_number(i);

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;
            if (j <= i) continue;

            int Z_j = system.atomic_number(j);
            Vec3 r_ij = neighbor_distance_vector(system, i, *it);
            Scalar r = r_ij.norm();

            if (r < 1e-10) continue;

            Vec3 r_hat = r_ij / r;

            try {
                const RepulsiveSpline& rep = materials_.get_rep_spline(Z_i, Z_j);
                Scalar dV_dr;
                rep.eval_deriv(r, dV_dr);

                Vec3 f = -dV_dr * r_hat;
                forces.row(i) -= f.transpose();
                forces.row(j) += f.transpose();

                // Stress contribution
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        stress(a, b) -= f[a] * r_ij[b];
                    }
                }
            } catch (...) {
                // No repulsive potential for this pair
            }
        }
    }
}

void DFTB::compute_scc_forces(const AtomicSystem& system,
                              const NeighborList& neighbors,
                              MatX3& forces) {
    Mat3 dummy_stress;
    compute_scc_forces_impl(system, neighbors, forces, dummy_stress, false);
}

void DFTB::compute_scc_forces_and_stress(const AtomicSystem& system,
                                         const NeighborList& neighbors,
                                         MatX3& forces, Mat3& stress) {
    compute_scc_forces_impl(system, neighbors, forces, stress, true);
}

void DFTB::compute_scc_forces_impl(const AtomicSystem& system,
                                   const NeighborList& neighbors,
                                   MatX3& forces, Mat3& stress,
                                   bool compute_stress) {
    const DenseHamiltonian& ham = hamiltonian_.hamiltonian();
    int nat = system.num_atoms();

    // Force from gamma derivative: F_I = -sum_J dq_I * dq_J * d(gamma_IJ)/dR_I
    for (std::size_t i = 0; i < static_cast<std::size_t>(nat); ++i) {
        Scalar U_i = elements_[ham.element_index[i]].hubbard_U;

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;
            if (j <= i) continue;

            Scalar U_j = elements_[ham.element_index[j]].hubbard_U;

            Vec3 r_ij = neighbor_distance_vector(system, i, *it);
            Scalar r = r_ij.norm();

            if (r < 1e-10) continue;

            Vec3 r_hat = r_ij / r;

            // Derivative of gamma
            Scalar dgamma_dr = gamma_derivative(r, U_i, U_j);

            // Force contribution
            // Note: d(gamma)/dR_I = d(gamma)/dr * r_hat (for r = R_J - R_I)
            Vec3 f = -ham.charges[i] * ham.charges[j] * dgamma_dr * r_hat;

            forces.row(i) -= f.transpose();
            forces.row(j) += f.transpose();

            // Stress contribution
            if (compute_stress) {
                for (int a = 0; a < 3; ++a) {
                    for (int b = 0; b < 3; ++b) {
                        stress(a, b) -= f[a] * r_ij[b];
                    }
                }
            }
        }
    }
}

} // namespace tb
} // namespace atomistica
