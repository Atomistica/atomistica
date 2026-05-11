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

#include <vector>

#include "../config.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Anderson mixer for self-consistent charge iterations
 *
 * Implements Anderson mixing to accelerate convergence in iterative processes.
 * See: V. Eyert, J. Comp. Phys. 124, 271 (1996)
 *
 * The mixer solves for optimal linear combination coefficients to minimize
 * the residual in the least-squares sense.
 */
class AndersonMixer {
public:
    /**
     * @brief Constructor
     *
     * @param memory Maximum number of history vectors to keep (M)
     */
    explicit AndersonMixer(int memory = 3) : memory_(memory) {}

    /**
     * @brief Reset the mixer (clear history)
     */
    void reset() {
        x_hist_.clear();
        F_hist_.clear();
        n_ = -1;
    }

    /**
     * @brief Set the dimension and reset history
     */
    void set_dimension(int n) {
        reset();
        n_ = n;
        x_hist_.reserve(memory_ + 1);
        F_hist_.reserve(memory_ + 1);
    }

    /**
     * @brief Perform one mixing iteration
     *
     * @param it     Current iteration number (1-based)
     * @param xi     In/Out: Current input vector (modified to next input)
     * @param yi     New output vector from solver
     * @param beta   Mixing parameter (0 < beta <= 1)
     * @param limit  Convergence threshold (optional)
     * @return true if converged (max|F| < limit), false otherwise
     */
    bool mix(int it, VecX& xi, const VecX& yi, Scalar beta, Scalar limit = 0.0);

    /**
     * @brief Get the current history size
     */
    int history_size() const { return static_cast<int>(F_hist_.size()); }

    /**
     * @brief Get the memory parameter
     */
    int memory() const { return memory_; }

private:
    int memory_;                    // Maximum history length (M)
    int n_ = -1;                    // Vector dimension
    std::vector<VecX> x_hist_;      // History of input vectors
    std::vector<VecX> F_hist_;      // History of residuals
};

} // namespace tb
} // namespace atomistica
