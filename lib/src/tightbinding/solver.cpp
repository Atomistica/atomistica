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

#include <atomistica/tightbinding/solver.hpp>

#include <cmath>
#include <stdexcept>

#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>

namespace atomistica {
namespace tb {

// ============================================================
// TBSolver implementation
// ============================================================

void TBSolver::solve(DenseHamiltonian& ham) {
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

void TBSolver::compute_occupation(DenseHamiltonian& ham, Scalar n_electrons,
                                  int spin_degeneracy) {
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

void TBSolver::build_density_matrix(DenseHamiltonian& ham) {
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

void TBSolver::compute_mulliken_charges(DenseHamiltonian& ham) {
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

void TBSolver::build_energy_weighted_density(DenseHamiltonian& ham) {
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

Scalar TBSolver::find_fermi_level(const VecX& eigenvalues, Scalar n_electrons,
                                  int spin_deg, Scalar kT, Scalar tol) {
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

// ============================================================
// PurificationSolver implementation
// ============================================================

void PurificationSolver::solve(DenseHamiltonian& ham, Scalar /*n_electrons*/,
                               int max_iter, Scalar tol) {
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

MatX PurificationSolver::compute_s_inv_sqrt(const MatX& S) {
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

} // namespace tb
} // namespace atomistica
