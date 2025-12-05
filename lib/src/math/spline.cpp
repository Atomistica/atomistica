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

#include "atomistica/math/spline.hpp"

#include <algorithm>
#include <cmath>

namespace atomistica {

// CubicSpline implementation

CubicSpline::CubicSpline(Scalar x_min, Scalar x_max, const std::vector<Scalar>& y) {
    init(x_min, x_max, y);
}

void CubicSpline::init(Scalar x_min, Scalar x_max, const std::vector<Scalar>& y) {
    if (y.size() < 2) {
        throw std::invalid_argument("Spline requires at least 2 points");
    }

    x_min_ = x_min;
    x_max_ = x_max;
    n_ = y.size();
    dx_ = (x_max - x_min) / static_cast<Scalar>(n_ - 1);
    inv_dx_ = 1.0 / dx_;

    compute_coefficients(y);
}

void CubicSpline::compute_coefficients(const std::vector<Scalar>& y) {
    // Compute second derivatives using natural spline boundary conditions
    // y''(x_min) = y''(x_max) = 0

    const std::size_t n = y.size();
    std::vector<Scalar> y2(n, 0.0);  // Second derivatives

    // Tridiagonal system for natural cubic spline with uniform spacing h:
    // y2[i-1] + 4*y2[i] + y2[i+1] = 6/h^2 * (y[i+1] - 2*y[i] + y[i-1])
    // with y2[0] = y2[n-1] = 0 (natural boundary conditions)

    // Thomas algorithm for tridiagonal system
    // a[i] * x[i-1] + b[i] * x[i] + c[i] * x[i+1] = d[i]
    // Here: a=1, b=4, c=1 (all uniform)

    std::vector<Scalar> c_star(n, 0.0);  // Modified upper diagonal
    std::vector<Scalar> d_star(n, 0.0);  // Modified RHS

    // Forward sweep (natural BC: y2[0] = 0 means we start from i=1)
    // First interior point (i=1): 4*y2[1] + y2[2] = rhs[1] (since y2[0]=0)
    Scalar rhs1 = (y[2] - 2.0 * y[1] + y[0]) * 6.0 * inv_dx_ * inv_dx_;
    c_star[1] = 1.0 / 4.0;
    d_star[1] = rhs1 / 4.0;

    for (std::size_t i = 2; i < n - 1; ++i) {
        Scalar rhs = (y[i + 1] - 2.0 * y[i] + y[i - 1]) * 6.0 * inv_dx_ * inv_dx_;
        Scalar denom = 4.0 - c_star[i - 1];  // b - a * c_star[i-1]
        c_star[i] = 1.0 / denom;
        d_star[i] = (rhs - d_star[i - 1]) / denom;
    }

    // Backsubstitution
    y2[n - 1] = 0.0;  // Natural boundary condition
    for (std::size_t k = n - 2; k >= 1; --k) {
        y2[k] = d_star[k] - c_star[k] * y2[k + 1];
    }
    y2[0] = 0.0;  // Natural boundary condition

    // Compute polynomial coefficients for each interval
    coeffs_.resize(4 * (n - 1));

    for (std::size_t i = 0; i < n - 1; ++i) {
        // Cubic: f(t) = a + b*t + c*t^2 + d*t^3, where t in [0, h]
        Scalar a = y[i];
        Scalar b = (y[i + 1] - y[i]) / dx_ - dx_ * (2.0 * y2[i] + y2[i + 1]) / 6.0;
        Scalar c = y2[i] / 2.0;
        Scalar d = (y2[i + 1] - y2[i]) / (6.0 * dx_);

        coeffs_[4 * i + 0] = a;
        coeffs_[4 * i + 1] = b;
        coeffs_[4 * i + 2] = c;
        coeffs_[4 * i + 3] = d;
    }
}

SplineResult CubicSpline::eval(Scalar x) const {
    if (!is_valid()) {
        return {0.0, 0.0};
    }

    // Clamp to valid range
    x = std::clamp(x, x_min_, x_max_);

    // Find interval
    Scalar t = (x - x_min_) * inv_dx_;
    std::size_t i = static_cast<std::size_t>(t);
    if (i >= n_ - 1) {
        i = n_ - 2;
    }

    // Local parameter within interval
    Scalar dt = (x - x_min_) - static_cast<Scalar>(i) * dx_;

    // Evaluate cubic polynomial and its derivative
    Scalar a = coeffs_[4 * i + 0];
    Scalar b = coeffs_[4 * i + 1];
    Scalar c = coeffs_[4 * i + 2];
    Scalar d = coeffs_[4 * i + 3];

    Scalar value = a + dt * (b + dt * (c + dt * d));
    Scalar deriv = b + dt * (2.0 * c + 3.0 * dt * d);

    return {value, deriv};
}

Scalar CubicSpline::value(Scalar x) const {
    return eval(x).value;
}

// NonUniformSpline implementation

NonUniformSpline::NonUniformSpline(const std::vector<Scalar>& x,
                                     const std::vector<Scalar>& y) {
    init(x, y);
}

void NonUniformSpline::init(const std::vector<Scalar>& x,
                            const std::vector<Scalar>& y) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("x and y must have same size");
    }
    if (x.size() < 2) {
        throw std::invalid_argument("Spline requires at least 2 points");
    }

    x_ = x;
    y_ = y;
    compute_coefficients(y);
}

void NonUniformSpline::compute_coefficients(const std::vector<Scalar>& y) {
    const std::size_t n = y.size();
    y2_.resize(n, 0.0);

    // Natural spline: y2[0] = y2[n-1] = 0
    // Solve tridiagonal system using Thomas algorithm
    std::vector<Scalar> c_star(n, 0.0);  // Modified upper diagonal ratio
    std::vector<Scalar> d_star(n, 0.0);  // Modified RHS

    // Forward sweep starting at i=1 (y2[0]=0 by natural BC)
    {
        Scalar h_prev = x_[1] - x_[0];
        Scalar h_next = x_[2] - x_[1];
        Scalar mu = h_prev / (h_prev + h_next);  // lower diagonal coefficient
        Scalar lambda = 1.0 - mu;                 // upper diagonal coefficient (= h_next/(h+h))
        // Diagonal is always 2 for natural spline
        // RHS: 6*f[x_{i-1}, x_i, x_{i+1}]
        Scalar rhs = 6.0 * ((y[2] - y[1]) / h_next - (y[1] - y[0]) / h_prev) / (h_prev + h_next);
        c_star[1] = lambda / 2.0;
        d_star[1] = rhs / 2.0;
    }

    for (std::size_t i = 2; i < n - 1; ++i) {
        Scalar h_prev = x_[i] - x_[i - 1];
        Scalar h_next = x_[i + 1] - x_[i];
        Scalar mu = h_prev / (h_prev + h_next);
        Scalar lambda = 1.0 - mu;
        Scalar rhs = 6.0 * ((y[i + 1] - y[i]) / h_next - (y[i] - y[i - 1]) / h_prev) / (h_prev + h_next);

        Scalar denom = 2.0 - mu * c_star[i - 1];
        c_star[i] = lambda / denom;
        d_star[i] = (rhs - mu * d_star[i - 1]) / denom;
    }

    // Backsubstitution
    y2_[n - 1] = 0.0;  // Natural BC
    for (std::size_t k = n - 2; k >= 1; --k) {
        y2_[k] = d_star[k] - c_star[k] * y2_[k + 1];
    }
    y2_[0] = 0.0;  // Natural BC
}

std::size_t NonUniformSpline::find_interval(Scalar x) const {
    // Binary search for interval
    auto it = std::lower_bound(x_.begin(), x_.end(), x);
    if (it == x_.begin()) {
        return 0;
    }
    if (it == x_.end()) {
        return x_.size() - 2;
    }
    return static_cast<std::size_t>(std::distance(x_.begin(), it)) - 1;
}

SplineResult NonUniformSpline::eval(Scalar x) const {
    if (!is_valid()) {
        return {0.0, 0.0};
    }

    // Clamp to range
    x = std::clamp(x, x_.front(), x_.back());

    std::size_t i = find_interval(x);

    Scalar h = x_[i + 1] - x_[i];
    Scalar a = (x_[i + 1] - x) / h;
    Scalar b = (x - x_[i]) / h;

    Scalar value = a * y_[i] + b * y_[i + 1] +
                   ((a * a * a - a) * y2_[i] + (b * b * b - b) * y2_[i + 1]) *
                   (h * h) / 6.0;

    Scalar deriv = (y_[i + 1] - y_[i]) / h -
                   (3.0 * a * a - 1.0) * h * y2_[i] / 6.0 +
                   (3.0 * b * b - 1.0) * h * y2_[i + 1] / 6.0;

    return {value, deriv};
}

Scalar NonUniformSpline::value(Scalar x) const {
    return eval(x).value;
}

} // namespace atomistica
