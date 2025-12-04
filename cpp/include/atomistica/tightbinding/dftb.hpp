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
#include <memory>
#include <set>
#include <vector>

#include "../config.hpp"
#include "../core/atomic_system.hpp"
#include "../core/neighbor_list.hpp"
#include "hamiltonian.hpp"
#include "materials.hpp"
#include "slater_koster.hpp"
#include "solver.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief DFTB (Density Functional Tight Binding) potential
 *
 * Implements non-orthogonal tight-binding with optional SCC
 * (Self-Consistent Charge) corrections.
 */
class DFTB {
public:
    /**
     * @brief Constructor
     *
     * @param skf_path Path to directory containing SKF files
     * @param enable_scc Enable self-consistent charges
     */
    explicit DFTB(const std::string& skf_path = "", bool enable_scc = false)
        : enable_scc_(enable_scc) {
        if (!skf_path.empty()) {
            materials_.load_skf_directory(skf_path);
        }
        hamiltonian_.set_materials(&materials_);
    }

    /**
     * @brief Get potential name
     */
    std::string name() const { return "DFTB"; }

    /**
     * @brief Get cutoff distance
     */
    Scalar cutoff() const { return materials_.get_max_cutoff(); }

    /**
     * @brief Add element to materials database
     */
    void add_element(const TBElementParams& elem) {
        materials_.add_element(elem);
        update_elements();
    }

    /**
     * @brief Load pair parameters from SKF file
     */
    void load_pair(int Z1, int Z2) {
        materials_.load_pair(Z1, Z2);
    }

    /**
     * @brief Set SKF directory
     */
    void set_skf_path(const std::string& path) {
        materials_.load_skf_directory(path);
    }

    /**
     * @brief Enable/disable SCC
     */
    void set_scc(bool enable) { enable_scc_ = enable; }

    /**
     * @brief Set SCC parameters
     */
    void set_scc_params(const SCCParams& params) { scc_params_ = params; }

    /**
     * @brief Set solver parameters
     */
    void set_solver_params(const SolverParams& params) {
        solver_params_ = params;
        solver_.set_params(params);
    }

    /**
     * @brief Initialize potential for atomic system
     */
    void init(const AtomicSystem& system) {
        // Collect unique elements
        elements_.clear();
        std::set<int> unique_Z;
        for (int i = 0; i < system.num_atoms(); ++i) {
            unique_Z.insert(system.atomic_number(i));
        }

        for (int Z : unique_Z) {
            if (materials_.has_element(Z)) {
                elements_.push_back(materials_.get_element(Z));
            }
        }

        // Initialize Hamiltonian
        hamiltonian_.init(system, elements_);

        // Count electrons
        n_electrons_ = 0.0;
        for (int i = 0; i < system.num_atoms(); ++i) {
            int Z = system.atomic_number(i);
            n_electrons_ += materials_.get_element(Z).valence_electrons;
        }

        initialized_ = true;
    }

    /**
     * @brief Compute energy
     */
    Scalar compute_energy(const AtomicSystem& system,
                         const NeighborList& neighbors) {
        if (!initialized_) init(system);

        // Build H and S matrices
        hamiltonian_.build_matrices(system, neighbors);

        if (enable_scc_) {
            // SCC iteration
            run_scc(system);
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

    /**
     * @brief Compute energy and forces
     */
    Scalar compute(const AtomicSystem& system, const NeighborList& neighbors,
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

    /**
     * @brief Get the Hamiltonian structure
     */
    DenseHamiltonian& hamiltonian() { return hamiltonian_.hamiltonian(); }
    const DenseHamiltonian& hamiltonian() const { return hamiltonian_.hamiltonian(); }

    /**
     * @brief Get eigenvalues
     */
    const VecX& eigenvalues() const { return hamiltonian_.hamiltonian().eigenvalues; }

    /**
     * @brief Get Fermi level
     */
    Scalar fermi_level() const { return hamiltonian_.hamiltonian().fermi_level; }

    /**
     * @brief Get Mulliken charges
     */
    const VecX& charges() const { return hamiltonian_.hamiltonian().charges; }

    /**
     * @brief Get band energy
     */
    Scalar band_energy() const { return hamiltonian_.hamiltonian().band_energy; }

    /**
     * @brief Get repulsive energy
     */
    Scalar repulsive_energy() const { return hamiltonian_.hamiltonian().repulsive_energy; }

    /**
     * @brief Get materials database
     */
    MaterialsDatabase& materials() { return materials_; }
    const MaterialsDatabase& materials() const { return materials_; }

private:
    MaterialsDatabase materials_;
    TBHamiltonian hamiltonian_;
    TBSolver solver_;

    std::vector<TBElementParams> elements_;
    Scalar n_electrons_ = 0.0;
    Scalar total_energy_ = 0.0;

    bool enable_scc_ = false;
    bool initialized_ = false;
    SCCParams scc_params_;
    SolverParams solver_params_;

    MatX gamma_;  // SCC gamma matrix

    /**
     * @brief Update elements list from materials database
     */
    void update_elements() {
        // Called when new elements are added
    }

    /**
     * @brief Solve electronic structure (single iteration)
     */
    void solve_electronic() {
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

    /**
     * @brief Run SCC iteration to self-consistency
     */
    void run_scc(const AtomicSystem& system) {
        DenseHamiltonian& ham = hamiltonian_.hamiltonian();
        int nat = ham.num_atoms;

        // Compute gamma matrix
        gamma_ = compute_gamma_matrix(system, elements_, ham.element_index);

        // Store original H matrix (without SCC correction)
        MatX H0 = ham.H;

        // Initialize charges
        VecX charges_old = VecX::Zero(nat);

        // Anderson mixing history
        std::vector<VecX> F_history, X_history;

        for (int iter = 0; iter < scc_params_.max_iterations; ++iter) {
            // Restore original H
            ham.H = H0;

            // Add SCC correction
            hamiltonian_.add_scc_correction(gamma_);

            // Solve
            solve_electronic();

            // Check convergence
            Scalar max_diff = (ham.charges - charges_old).cwiseAbs().maxCoeff();

            if (max_diff < scc_params_.convergence_threshold) {
                break;  // Converged
            }

            // Update charges with mixing
            if (scc_params_.anderson_memory > 0 && iter > 0) {
                // Anderson mixing
                VecX F = ham.charges - charges_old;

                F_history.push_back(F);
                X_history.push_back(charges_old);

                if (static_cast<int>(F_history.size()) > scc_params_.anderson_memory) {
                    F_history.erase(F_history.begin());
                    X_history.erase(X_history.begin());
                }

                // Simple Anderson mixing with one history point
                if (F_history.size() >= 2) {
                    VecX dF = F_history.back() - F_history[F_history.size()-2];
                    VecX dX = X_history.back() - X_history[X_history.size()-2];

                    Scalar beta = F.dot(dF) / dF.squaredNorm();
                    beta = std::min(std::max(beta, 0.0), 1.0);

                    charges_old = (1.0 - beta) * (charges_old + scc_params_.mixing_parameter * F)
                                + beta * (X_history.back() + scc_params_.mixing_parameter * F_history.back());
                } else {
                    charges_old = charges_old + scc_params_.mixing_parameter * F;
                }
            } else {
                // Simple linear mixing
                charges_old = (1.0 - scc_params_.mixing_parameter) * charges_old
                            + scc_params_.mixing_parameter * ham.charges;
            }

            ham.charges = charges_old;
        }
    }

    /**
     * @brief Compute SCC energy correction
     */
    Scalar compute_scc_energy() const {
        const DenseHamiltonian& ham = hamiltonian_.hamiltonian();
        int nat = ham.num_atoms;

        // E_scc = 0.5 * sum_ij gamma_ij * dq_i * dq_j
        Scalar E_scc = 0.0;
        for (int i = 0; i < nat; ++i) {
            for (int j = 0; j < nat; ++j) {
                E_scc += 0.5 * gamma_(i, j) * ham.charges[i] * ham.charges[j];
            }
        }

        return E_scc;
    }

    /**
     * @brief Compute band structure forces using Hellmann-Feynman theorem
     *
     * F_I = -Tr(rho * dH/dR_I) + Tr(E * dS/dR_I)
     *
     * where E is the energy-weighted density matrix
     */
    void compute_band_forces(const AtomicSystem& system,
                            const NeighborList& neighbors,
                            MatX3& forces) {
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

                forces.row(i) -= force_ij.transpose();
                forces.row(j) += force_ij.transpose();
            }
        }
    }

    /**
     * @brief Compute SCC force corrections
     */
    void compute_scc_forces(const AtomicSystem& system,
                           const NeighborList& /*neighbors*/,
                           MatX3& forces) {
        const DenseHamiltonian& ham = hamiltonian_.hamiltonian();
        int nat = system.num_atoms();

        // SCC forces come from:
        // 1. Derivative of gamma matrix: F_I += -0.5 * sum_J dq_I * dq_J * dgamma_IJ/dR_I
        // 2. Derivative of Mulliken charges (implicit through rho)

        // For now, only include explicit gamma derivative contribution
        for (int i = 0; i < nat; ++i) {
            for (int j = i + 1; j < nat; ++j) {
                Vec3 r_ij = system.position(j) - system.position(i);
                Scalar r = r_ij.norm();

                if (r < 1e-10) continue;

                Vec3 r_hat = r_ij / r;

                // Derivative of gamma (simplified - full version would use DFTB formula)
                Scalar U_i = elements_[ham.element_index[i]].hubbard_U;
                Scalar U_j = elements_[ham.element_index[j]].hubbard_U;

                Scalar tau_i = 3.2 * U_i * U_i;
                Scalar tau_j = 3.2 * U_j * U_j;
                Scalar tau = std::sqrt(tau_i * tau_j);

                // d(gamma)/dr = d(erf(tau*r)/r)/dr
                Scalar dgamma_dr;
                if (r < 0.1) {
                    dgamma_dr = -2.0 * tau * tau * tau / (3.0 * std::sqrt(M_PI));
                } else {
                    Scalar x = tau * r;
                    Scalar erf_val = std::erf(x);
                    Scalar exp_val = std::exp(-x * x);
                    dgamma_dr = (2.0 * tau * exp_val / std::sqrt(M_PI) - erf_val / r) / r;
                }

                // Force contribution
                Vec3 f = -0.5 * ham.charges[i] * ham.charges[j] * dgamma_dr * r_hat;
                forces.row(i) -= f.transpose();
                forces.row(j) += f.transpose();
            }
        }
    }
};

} // namespace tb
} // namespace atomistica
