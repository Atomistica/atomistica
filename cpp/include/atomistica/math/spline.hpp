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

#include <cstddef>
#include <stdexcept>
#include <vector>

#include "../config.hpp"

namespace atomistica {

/**
 * @brief Result of spline evaluation
 */
struct SplineResult {
    Scalar value;       // f(x)
    Scalar derivative;  // f'(x)
};

/**
 * @brief Cubic spline interpolation on uniform grid
 *
 * Natural cubic spline with continuous second derivatives.
 * Efficient evaluation using uniform spacing.
 */
class CubicSpline {
public:
    CubicSpline() = default;

    /**
     * @brief Construct spline from data points
     *
     * @param x_min Lower bound of interpolation range
     * @param x_max Upper bound of interpolation range
     * @param y Values at uniformly spaced points (including endpoints)
     */
    CubicSpline(Scalar x_min, Scalar x_max, const std::vector<Scalar>& y);

    /**
     * @brief Initialize spline from data points
     */
    void init(Scalar x_min, Scalar x_max, const std::vector<Scalar>& y);

    /**
     * @brief Check if spline is initialized
     */
    bool is_valid() const { return !coeffs_.empty(); }

    /**
     * @brief Evaluate spline at point x
     */
    SplineResult eval(Scalar x) const;

    /**
     * @brief Evaluate just the value (faster if derivative not needed)
     */
    Scalar value(Scalar x) const;

    /**
     * @brief Get interpolation range
     */
    Scalar x_min() const { return x_min_; }
    Scalar x_max() const { return x_max_; }

    /**
     * @brief Number of data points
     */
    std::size_t size() const { return n_; }

private:
    void compute_coefficients(const std::vector<Scalar>& y);

    Scalar x_min_ = 0.0;
    Scalar x_max_ = 0.0;
    Scalar dx_ = 0.0;      // Grid spacing
    Scalar inv_dx_ = 0.0;  // 1/dx for efficiency
    std::size_t n_ = 0;    // Number of points

    // Spline coefficients: f(x) = a + b*t + c*t^2 + d*t^3, where t is local parameter
    // Stored as [a0, b0, c0, d0, a1, b1, c1, d1, ...]
    std::vector<Scalar> coeffs_;
};

/**
 * @brief Cubic spline on non-uniform grid
 *
 * For tabulated data with irregular spacing (e.g., Slater-Koster tables)
 */
class NonUniformSpline {
public:
    NonUniformSpline() = default;

    /**
     * @brief Construct spline from data points
     *
     * @param x X coordinates (must be monotonically increasing)
     * @param y Y values at each x
     */
    NonUniformSpline(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

    /**
     * @brief Initialize spline from data points
     */
    void init(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

    /**
     * @brief Check if spline is initialized
     */
    bool is_valid() const { return !x_.empty(); }

    /**
     * @brief Evaluate spline at point x
     */
    SplineResult eval(Scalar x) const;

    /**
     * @brief Evaluate just the value
     */
    Scalar value(Scalar x) const;

    /**
     * @brief Get interpolation range
     */
    Scalar x_min() const { return x_.empty() ? 0.0 : x_.front(); }
    Scalar x_max() const { return x_.empty() ? 0.0 : x_.back(); }

private:
    std::size_t find_interval(Scalar x) const;
    void compute_coefficients(const std::vector<Scalar>& y);

    std::vector<Scalar> x_;
    std::vector<Scalar> y_;
    std::vector<Scalar> y2_;  // Second derivatives for natural spline
};

} // namespace atomistica
