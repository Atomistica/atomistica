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

#include <atomistica/tightbinding/anderson_mixer.hpp>

#include <algorithm>
#include <cmath>

#include <Eigen/Cholesky>

namespace atomistica {
namespace tb {

bool AndersonMixer::mix(int it, VecX& xi, const VecX& yi, Scalar beta, Scalar limit) {
    int n = xi.size();

    // Handle empty vectors - trivially converged
    if (n == 0) {
        return true;
    }

    // Initialize if needed
    if (n_ < n) {
        set_dimension(n);
    }

    // Current residual F = y - x
    VecX F = yi - xi;

    // Check convergence
    Scalar max_residual = F.cwiseAbs().maxCoeff();
    bool converged = (limit > 0.0) && (max_residual < limit);

    // Determine effective history size M
    int M = std::min(it - 1, memory_);
    M = std::min(M, static_cast<int>(F_hist_.size()));

    VecX x_new;

    if (M > 0) {
        // Build and solve the linear system A * z = b
        // A_ij = <F0 - Fi, F0 - Fj>
        // b_i  = <F0 - Fi, F0>
        MatX A = MatX::Zero(M, M);
        VecX b = VecX::Zero(M);

        for (int i = 0; i < M; ++i) {
            VecX dF_i = F - F_hist_[i];
            b(i) = dF_i.dot(F);

            for (int j = 0; j < M; ++j) {
                VecX dF_j = F - F_hist_[j];
                A(i, j) = dF_i.dot(dF_j);
            }
        }

        // Solve A * z = b using Eigen's built-in solver
        // Use LDLT for symmetric positive semi-definite (fallback to simple mixing if singular)
        Eigen::LDLT<MatX> ldlt(A);

        if (ldlt.info() == Eigen::Success && ldlt.isPositive()) {
            VecX z = ldlt.solve(b);

            // Compute optimal linear combination
            // xbar = x0 + sum_j z_j * (x_j - x0)
            // Fbar = F0 + sum_j z_j * (F_j - F0)
            VecX xbar = xi;
            VecX Fbar = F;

            for (int j = 0; j < M; ++j) {
                xbar += z(j) * (x_hist_[j] - xi);
                Fbar += z(j) * (F_hist_[j] - F);
            }

            // New input: xbar + beta * Fbar
            x_new = xbar + beta * Fbar;
        } else {
            // Matrix singular - fall back to simple mixing
            x_new = (1.0 - beta) * xi + beta * yi;
        }
    } else {
        // No history yet - use simple mixing
        x_new = (1.0 - beta) * xi + beta * yi;
    }

    // Shift history: insert current at front, remove oldest if necessary
    F_hist_.insert(F_hist_.begin(), F);
    x_hist_.insert(x_hist_.begin(), xi);

    if (static_cast<int>(F_hist_.size()) > memory_) {
        F_hist_.pop_back();
        x_hist_.pop_back();
    }

    // Update xi with new mixed value
    xi = x_new;

    return converged;
}

} // namespace tb
} // namespace atomistica
