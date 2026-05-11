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
#include <memory>
#include <string>
#include <vector>

#include "../config.hpp"
#include "../core/atomic_system.hpp"
#include "../core/neighbor_list.hpp"
#include "anderson_mixer.hpp"
#include "hamiltonian.hpp"
#include "materials.hpp"
#include "slater_koster.hpp"
#include "solver.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Compute short-range gamma function for SCC-DFTB
 *
 * The gamma function describes the Coulomb interaction between
 * atomic charge distributions. For DFTB2:
 *
 *   gamma_ij = 1/r * erf(C_ij * r)
 *
 * where C_ij depends on the Hubbard U parameters.
 *
 * @param r Distance between atoms
 * @param U_i Hubbard U of atom i
 * @param U_j Hubbard U of atom j
 * @return gamma_ij value
 */
inline Scalar gamma_function(Scalar r, Scalar U_i, Scalar U_j) {
    // tau = 16/5 * U (in atomic units, approximately 3.2 * U in eV)
    // Using the DFTB convention
    Scalar tau_i = 3.2 * U_i * U_i;  // tau^2 actually
    Scalar tau_j = 3.2 * U_j * U_j;
    Scalar tau = std::sqrt(tau_i * tau_j);

    if (r < 1e-10) {
        // On-site: gamma_ii = U_i
        return 0.5 * (U_i + U_j);
    }

    Scalar C = std::sqrt(tau);
    Scalar x = C * r;
    Scalar erf_val = std::erf(x);
    return erf_val / r;
}

/**
 * @brief Compute derivative of gamma function with respect to r
 *
 * d(gamma)/dr = (2*C*exp(-C^2*r^2)/sqrt(pi) - erf(C*r)/r) / r
 *
 * @param r Distance between atoms
 * @param U_i Hubbard U of atom i
 * @param U_j Hubbard U of atom j
 * @return d(gamma_ij)/dr
 */
inline Scalar gamma_derivative(Scalar r, Scalar U_i, Scalar U_j) {
    Scalar tau_i = 3.2 * U_i * U_i;
    Scalar tau_j = 3.2 * U_j * U_j;
    Scalar tau = std::sqrt(tau_i * tau_j);
    Scalar C = std::sqrt(tau);

    if (r < 1e-6) {
        // Small r expansion: gamma ≈ 2C/sqrt(pi) - 2C^3*r^2/(3*sqrt(pi)) + ...
        // d(gamma)/dr ≈ -4C^3*r/(3*sqrt(pi)) + ...
        return -4.0 * C * C * C * r / (3.0 * std::sqrt(M_PI));
    }

    Scalar x = C * r;
    Scalar erf_val = std::erf(x);
    Scalar exp_val = std::exp(-x * x);

    return (2.0 * C * exp_val / std::sqrt(M_PI) - erf_val / r) / r;
}

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
        : enable_scc_(enable_scc), mixer_(3) {
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
    void set_scc_params(const SCCParams& params) {
        scc_params_ = params;
        mixer_ = AndersonMixer(params.anderson_memory);
    }

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
    void init(const AtomicSystem& system);

    /**
     * @brief Compute energy
     */
    Scalar compute_energy(const AtomicSystem& system,
                         const NeighborList& neighbors);

    /**
     * @brief Compute energy and forces
     */
    Scalar compute(const AtomicSystem& system, const NeighborList& neighbors,
                  MatX3& forces);

    /**
     * @brief Compute energy, forces, and stress tensor
     */
    Scalar compute_with_stress(const AtomicSystem& system, const NeighborList& neighbors,
                               MatX3& forces, Mat3& stress);

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

    /**
     * @brief Get gamma matrix (for SCC)
     */
    const MatX& gamma_matrix() const { return gamma_; }

    /**
     * @brief Get number of SCC iterations from last computation
     */
    int scc_iterations() const { return last_scc_iterations_; }

private:
    MaterialsDatabase materials_;
    TBHamiltonian hamiltonian_;
    TBSolver solver_;

    std::vector<TBElementParams> elements_;
    Scalar n_electrons_ = 0.0;
    Scalar total_energy_ = 0.0;

    bool enable_scc_ = false;
    AndersonMixer mixer_;
    bool initialized_ = false;
    SCCParams scc_params_;
    SolverParams solver_params_;

    MatX gamma_;                    // SCC gamma matrix
    MatX H0_;                       // Original H matrix (without SCC)
    int last_scc_iterations_ = 0;   // Number of SCC iterations in last computation

    /**
     * @brief Update elements list from materials database
     */
    void update_elements() {
        // Called when new elements are added
    }

    /**
     * @brief Solve electronic structure (single iteration)
     */
    void solve_electronic();

    /**
     * @brief Compute gamma matrix for SCC-DFTB using neighbor list
     */
    void compute_gamma(const AtomicSystem& system, const NeighborList& neighbors);

    /**
     * @brief Run SCC iteration to self-consistency using Anderson mixing
     */
    void run_scc(const AtomicSystem& system, const NeighborList& neighbors);

    /**
     * @brief Compute SCC energy correction
     *
     * E_scc = 0.5 * sum_ij gamma_ij * dq_i * dq_j
     */
    Scalar compute_scc_energy() const;

    /**
     * @brief Compute band structure forces using Hellmann-Feynman theorem
     *
     * F_I = -Tr(rho * dH/dR_I) + Tr(E * dS/dR_I)
     *
     * where E is the energy-weighted density matrix
     */
    void compute_band_forces(const AtomicSystem& system,
                            const NeighborList& neighbors,
                            MatX3& forces);

    /**
     * @brief Compute band structure forces and stress tensor
     */
    void compute_band_forces_and_stress(const AtomicSystem& system,
                                        const NeighborList& neighbors,
                                        MatX3& forces, Mat3& stress);

    /**
     * @brief Implementation of band force/stress calculation
     */
    void compute_band_forces_impl(const AtomicSystem& system,
                                  const NeighborList& neighbors,
                                  MatX3& forces, Mat3& stress,
                                  bool compute_stress);

    /**
     * @brief Compute repulsive forces and stress
     */
    void compute_repulsive_forces_and_stress(const AtomicSystem& system,
                                             const NeighborList& neighbors,
                                             MatX3& forces, Mat3& stress);

    /**
     * @brief Compute SCC force corrections
     *
     * The SCC force has two contributions:
     * 1. From d(gamma)/dR: F_I = -sum_J dq_I * dq_J * d(gamma_IJ)/dR_I
     * 2. From d(shift*S)/dR in Hamiltonian (already included in band forces via rho * dH)
     */
    void compute_scc_forces(const AtomicSystem& system,
                           const NeighborList& neighbors,
                           MatX3& forces);

    /**
     * @brief Compute SCC force and stress corrections
     */
    void compute_scc_forces_and_stress(const AtomicSystem& system,
                                       const NeighborList& neighbors,
                                       MatX3& forces, Mat3& stress);

    /**
     * @brief Implementation of SCC force/stress calculation
     */
    void compute_scc_forces_impl(const AtomicSystem& system,
                                 const NeighborList& neighbors,
                                 MatX3& forces, Mat3& stress,
                                 bool compute_stress);
};

} // namespace tb
} // namespace atomistica
