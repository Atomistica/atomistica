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
#include <stdexcept>
#include <vector>

#include <Eigen/Eigenvalues>

#include "../config.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Fermi-Dirac distribution
 *
 * @param e Energy
 * @param mu Chemical potential (Fermi level)
 * @param kT Temperature in energy units
 * @return Occupation number (0 to 1)
 */
inline Scalar fermi_dirac(Scalar e, Scalar mu, Scalar kT) {
    if (kT < 1e-10) {
        // Zero temperature: step function
        return (e < mu) ? 1.0 : ((e > mu) ? 0.0 : 0.5);
    }

    Scalar x = (e - mu) / kT;

    // Avoid overflow
    if (x > 40.0) return 0.0;
    if (x < -40.0) return 1.0;

    return 1.0 / (1.0 + std::exp(x));
}

/**
 * @brief Derivative of Fermi-Dirac distribution
 */
inline Scalar fermi_dirac_derivative(Scalar e, Scalar mu, Scalar kT) {
    if (kT < 1e-10) return 0.0;

    Scalar x = (e - mu) / kT;
    if (std::abs(x) > 40.0) return 0.0;

    Scalar f = fermi_dirac(e, mu, kT);
    return -f * (1.0 - f) / kT;
}

/**
 * @brief Electronic entropy contribution
 *
 * S_el = -k_B * sum_i [f_i * ln(f_i) + (1-f_i) * ln(1-f_i)]
 */
inline Scalar electronic_entropy(Scalar f) {
    if (f < 1e-15 || f > 1.0 - 1e-15) return 0.0;
    return -f * std::log(f) - (1.0 - f) * std::log(1.0 - f);
}

/**
 * @brief Tight-binding eigenvalue solver using LAPACK
 */
class TBSolver {
public:
    TBSolver() = default;

    /**
     * @brief Set solver parameters
     */
    void set_params(const SolverParams& params) { params_ = params; }

    /**
     * @brief Solve generalized eigenvalue problem H*C = S*C*E
     *
     * Finds eigenvalues and eigenvectors of the tight-binding
     * Hamiltonian with overlap using Eigen's GeneralizedSelfAdjointEigenSolver.
     *
     * @param ham Hamiltonian structure (H, S matrices modified on output)
     */
    void solve(DenseHamiltonian& ham) {
        int n = ham.num_orbitals;
        if (n == 0) return;

        // Use Eigen's generalized eigenvalue solver
        // Solves H*C = S*C*E where H and S are symmetric/self-adjoint
        Eigen::GeneralizedSelfAdjointEigenSolver<MatX> solver(ham.H, ham.S);

        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Generalized eigenvalue solver failed");
        }

        // Store results
        ham.eigenvalues = solver.eigenvalues();
        ham.eigenvectors = solver.eigenvectors();
    }

    /**
     * @brief Compute occupation numbers using Fermi-Dirac distribution
     *
     * @param ham Hamiltonian with eigenvalues
     * @param n_electrons Total number of electrons
     * @param spin_degeneracy Spin degeneracy (1 or 2)
     */
    void compute_occupation(DenseHamiltonian& ham, Scalar n_electrons,
                           int spin_degeneracy = 2) {
        int n = ham.num_orbitals;
        Scalar kT = params_.electronic_temperature;

        // Find Fermi level using bisection
        ham.fermi_level = find_fermi_level(ham.eigenvalues, n_electrons,
                                           spin_degeneracy, kT);

        // Compute occupation numbers
        ham.occupation.resize(n);
        for (int i = 0; i < n; ++i) {
            ham.occupation[i] = spin_degeneracy *
                               fermi_dirac(ham.eigenvalues[i], ham.fermi_level, kT);
        }
    }

    /**
     * @brief Build density matrix from eigenvectors and occupations
     *
     * rho = C * diag(f) * C^T where C are eigenvectors
     */
    void build_density_matrix(DenseHamiltonian& ham) {
        int n = ham.num_orbitals;
        ham.rho = MatX::Zero(n, n);

        for (int k = 0; k < n; ++k) {
            Scalar f_k = ham.occupation[k];
            if (f_k < 1e-15) continue;  // Skip unoccupied states

            for (int i = 0; i < n; ++i) {
                for (int j = i; j < n; ++j) {
                    Scalar contrib = f_k * ham.eigenvectors(i, k) *
                                    ham.eigenvectors(j, k);
                    ham.rho(i, j) += contrib;
                    if (i != j) ham.rho(j, i) += contrib;
                }
            }
        }
    }

    /**
     * @brief Compute band energy from eigenvalues and occupations
     *
     * E_band = sum_i f_i * epsilon_i
     */
    Scalar compute_band_energy(const DenseHamiltonian& ham) const {
        int n = ham.num_orbitals;
        Scalar E_band = 0.0;

        for (int i = 0; i < n; ++i) {
            E_band += ham.occupation[i] * ham.eigenvalues[i];
        }

        return E_band;
    }

    /**
     * @brief Compute electronic free energy (including entropy)
     *
     * A = E_band - T*S_el
     */
    Scalar compute_free_energy(const DenseHamiltonian& ham) const {
        Scalar E_band = compute_band_energy(ham);
        Scalar kT = params_.electronic_temperature;

        if (kT < 1e-10) return E_band;

        // Electronic entropy
        Scalar S_el = 0.0;
        int n = ham.num_orbitals;
        for (int i = 0; i < n; ++i) {
            Scalar f = ham.occupation[i] / 2.0;  // Per spin channel
            S_el += 2.0 * electronic_entropy(f);  // Both spin channels
        }

        return E_band - kT * S_el;
    }

    /**
     * @brief Compute Mulliken charges
     *
     * q_i = sum_a (rho * S)_{a,a} for orbitals a on atom i
     */
    void compute_mulliken_charges(DenseHamiltonian& ham) {
        int nat = ham.num_atoms;
        ham.charges = VecX::Zero(nat);

        // Compute rho * S
        MatX rhoS = ham.rho * ham.S;

        for (int i = 0; i < nat; ++i) {
            int offset = ham.orbital_offset[i];
            int norb = ham.orbitals_per_atom[i];

            Scalar q = 0.0;
            for (int a = 0; a < norb; ++a) {
                q += rhoS(offset + a, offset + a);
            }

            // Net charge = q0 - q (positive = electron deficient)
            ham.charges[i] = ham.neutral_charges[i] - q;
        }
    }

    /**
     * @brief Build E matrix for force calculation
     *
     * E_ij = sum_k f_k * epsilon_k * C_ik * C_jk
     */
    void build_energy_weighted_density(DenseHamiltonian& ham) {
        int n = ham.num_orbitals;
        ham.e_matrix = MatX::Zero(n, n);

        for (int k = 0; k < n; ++k) {
            Scalar fe = ham.occupation[k] * ham.eigenvalues[k];
            if (std::abs(fe) < 1e-15) continue;

            for (int i = 0; i < n; ++i) {
                for (int j = i; j < n; ++j) {
                    Scalar contrib = fe * ham.eigenvectors(i, k) *
                                    ham.eigenvectors(j, k);
                    ham.e_matrix(i, j) += contrib;
                    if (i != j) ham.e_matrix(j, i) += contrib;
                }
            }
        }
    }

private:
    SolverParams params_;

    /**
     * @brief Find Fermi level using bisection
     */
    Scalar find_fermi_level(const VecX& eigenvalues, Scalar n_electrons,
                            int spin_deg, Scalar kT, Scalar tol = 1e-12) {
        int n = eigenvalues.size();
        if (n == 0) return 0.0;

        // Initial bounds
        Scalar mu_lo = eigenvalues[0] - 10.0 * kT;
        Scalar mu_hi = eigenvalues[n-1] + 10.0 * kT;

        // Bisection
        const int max_iter = 100;
        for (int iter = 0; iter < max_iter; ++iter) {
            Scalar mu = 0.5 * (mu_lo + mu_hi);

            // Count electrons at this chemical potential
            Scalar n_el = 0.0;
            for (int i = 0; i < n; ++i) {
                n_el += spin_deg * fermi_dirac(eigenvalues[i], mu, kT);
            }

            if (n_el < n_electrons) {
                mu_lo = mu;
            } else {
                mu_hi = mu;
            }

            if (mu_hi - mu_lo < tol) break;
        }

        return 0.5 * (mu_lo + mu_hi);
    }
};

/**
 * @brief Canonical purification solver (O(N) scaling alternative)
 *
 * Builds density matrix directly without diagonalization.
 * Uses iterative purification: rho_{n+1} = 3*rho_n^2 - 2*rho_n^3
 */
class PurificationSolver {
public:
    PurificationSolver() = default;

    void set_params(const SolverParams& params) { params_ = params; }

    /**
     * @brief Solve using canonical purification
     *
     * @param ham Hamiltonian structure
     * @param n_electrons Target number of electrons (currently unused, for future extension)
     * @param max_iter Maximum iterations
     * @param tol Convergence tolerance
     */
    void solve(DenseHamiltonian& ham, Scalar /*n_electrons*/,
               int max_iter = 100, Scalar tol = 1e-8) {
        int n = ham.num_orbitals;

        // Compute S^(-1/2) for orthogonalization
        // Use Cholesky: S = L * L^T, then S^(-1/2) = L^(-T)
        MatX S_inv_sqrt = compute_s_inv_sqrt(ham.S);

        // Transform H to orthogonal basis: H' = S^(-1/2)^T * H * S^(-1/2)
        MatX H_orth = S_inv_sqrt.transpose() * ham.H * S_inv_sqrt;

        // Estimate spectral bounds using SelfAdjointEigenSolver
        Eigen::SelfAdjointEigenSolver<MatX> es(H_orth, Eigen::EigenvaluesOnly);
        Scalar e_min = es.eigenvalues().minCoeff();
        Scalar e_max = es.eigenvalues().maxCoeff();

        // Scale H to [0, 1] interval
        Scalar scale = 1.0 / (e_max - e_min);
        MatX rho = MatX::Identity(n, n) - scale * (H_orth - e_min * MatX::Identity(n, n));

        // Purification iterations
        for (int iter = 0; iter < max_iter; ++iter) {
            // McWeeny purification: rho = 3*rho^2 - 2*rho^3
            MatX rho2 = rho * rho;
            MatX rho_new = 3.0 * rho2 - 2.0 * rho2 * rho;

            // Check convergence
            Scalar diff = (rho_new - rho).norm();
            rho = rho_new;

            if (diff < tol) break;
        }

        // Transform back to non-orthogonal basis
        ham.rho = S_inv_sqrt * rho * S_inv_sqrt.transpose();

        // Compute band energy: E = Tr(rho * H)
        ham.band_energy = (ham.rho * ham.H).trace();
    }

private:
    SolverParams params_;

    /**
     * @brief Compute S^(-1/2) using Cholesky decomposition
     */
    MatX compute_s_inv_sqrt(const MatX& S) {
        // Cholesky decomposition: S = L * L^T
        Eigen::LLT<MatX> llt(S);
        if (llt.info() != Eigen::Success) {
            throw std::runtime_error("Overlap matrix not positive definite");
        }

        MatX L = llt.matrixL();

        // S^(-1/2) = L^(-T)
        MatX L_inv = L.inverse();
        return L_inv.transpose();
    }
};

} // namespace tb
} // namespace atomistica
