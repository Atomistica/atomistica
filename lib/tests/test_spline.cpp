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

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <atomistica/math/spline.hpp>

#include <cmath>
#include <vector>

using namespace atomistica;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("CubicSpline interpolation", "[Spline]") {
    SECTION("Linear function") {
        // f(x) = 2x + 1
        std::vector<Scalar> y = {1.0, 3.0, 5.0, 7.0, 9.0};
        CubicSpline spline(0.0, 4.0, y);

        REQUIRE(spline.is_valid());

        // Test at data points
        REQUIRE_THAT(spline.value(0.0), WithinRel(1.0, 1e-6));
        REQUIRE_THAT(spline.value(2.0), WithinRel(5.0, 1e-6));
        REQUIRE_THAT(spline.value(4.0), WithinRel(9.0, 1e-6));

        // Test interpolation
        REQUIRE_THAT(spline.value(1.5), WithinRel(4.0, 1e-6));
        REQUIRE_THAT(spline.value(3.5), WithinRel(8.0, 1e-6));

        // Test derivative
        auto result = spline.eval(2.0);
        REQUIRE_THAT(result.derivative, WithinRel(2.0, 1e-4));
    }

    SECTION("Quadratic function") {
        // f(x) = x^2, x in [0, 4]
        std::vector<Scalar> y;
        for (int i = 0; i <= 10; ++i) {
            Scalar x = i * 0.4;
            y.push_back(x * x);
        }
        CubicSpline spline(0.0, 4.0, y);

        // Test interpolation
        REQUIRE_THAT(spline.value(1.5), WithinRel(2.25, 1e-3));
        REQUIRE_THAT(spline.value(2.5), WithinRel(6.25, 1e-3));

        // Test derivative: f'(x) = 2x
        auto result = spline.eval(2.0);
        REQUIRE_THAT(result.derivative, WithinRel(4.0, 1e-2));
    }

    SECTION("Sine function") {
        // f(x) = sin(x), x in [0, 2*pi]
        const int n = 21;
        std::vector<Scalar> y(n);
        for (int i = 0; i < n; ++i) {
            Scalar x = i * 2.0 * PI / (n - 1);
            y[i] = std::sin(x);
        }
        CubicSpline spline(0.0, 2.0 * PI, y);

        // Test at pi/2
        REQUIRE_THAT(spline.value(PI / 2), WithinRel(1.0, 1e-3));

        // Test derivative at 0: cos(0) = 1
        auto result = spline.eval(0.0);
        REQUIRE_THAT(result.derivative, WithinRel(1.0, 1e-2));

        // Test derivative at pi/2: cos(pi/2) = 0
        result = spline.eval(PI / 2);
        REQUIRE_THAT(result.derivative, WithinAbs(0.0, 0.05));
    }
}

TEST_CASE("CubicSpline edge cases", "[Spline]") {
    std::vector<Scalar> y = {1.0, 4.0, 9.0, 16.0, 25.0};  // x^2 for x = 1..5
    CubicSpline spline(1.0, 5.0, y);

    SECTION("Below range") {
        // Should clamp to x_min
        REQUIRE_THAT(spline.value(0.0), WithinRel(1.0, 1e-6));
        REQUIRE_THAT(spline.value(-10.0), WithinRel(1.0, 1e-6));
    }

    SECTION("Above range") {
        // Should clamp to x_max
        REQUIRE_THAT(spline.value(6.0), WithinRel(25.0, 1e-6));
        REQUIRE_THAT(spline.value(100.0), WithinRel(25.0, 1e-6));
    }
}

TEST_CASE("NonUniformSpline interpolation", "[Spline]") {
    SECTION("Linear function") {
        std::vector<Scalar> x = {0.0, 1.0, 3.0, 6.0, 10.0};
        std::vector<Scalar> y = {0.0, 2.0, 6.0, 12.0, 20.0};  // f(x) = 2x

        NonUniformSpline spline(x, y);

        REQUIRE(spline.is_valid());

        // Test at data points
        REQUIRE_THAT(spline.value(0.0), WithinRel(0.0, 1e-6));
        REQUIRE_THAT(spline.value(3.0), WithinRel(6.0, 1e-6));

        // Test interpolation
        REQUIRE_THAT(spline.value(2.0), WithinRel(4.0, 1e-3));
        REQUIRE_THAT(spline.value(5.0), WithinRel(10.0, 1e-3));

        // Test derivative
        auto result = spline.eval(5.0);
        REQUIRE_THAT(result.derivative, WithinRel(2.0, 1e-2));
    }

    SECTION("Non-uniform spacing") {
        // Dense near x=0, sparse for larger x
        std::vector<Scalar> x = {0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0};
        std::vector<Scalar> y;
        for (Scalar xi : x) {
            y.push_back(std::exp(-xi));
        }

        NonUniformSpline spline(x, y);

        // Test interpolation in well-sampled region
        REQUIRE_THAT(spline.value(0.05), WithinRel(std::exp(-0.05), 1e-2));
        REQUIRE_THAT(spline.value(0.15), WithinRel(std::exp(-0.15), 1e-2));
        // x=3.0 is between widely-spaced points (2.0, 5.0) - relax tolerance
        // Cubic spline interpolation on sparse exp decay has significant error
        REQUIRE_THAT(spline.value(3.0), WithinRel(std::exp(-3.0), 0.25));
    }
}

TEST_CASE("Spline numerical derivative", "[Spline]") {
    // Test that analytical derivative matches numerical
    std::vector<Scalar> y;
    for (int i = 0; i <= 20; ++i) {
        Scalar x = i * 0.5;
        y.push_back(std::sin(x));
    }
    CubicSpline spline(0.0, 10.0, y);

    const Scalar dx = 1e-6;

    for (Scalar x = 0.5; x < 9.5; x += 0.7) {
        auto result = spline.eval(x);

        // Numerical derivative
        Scalar f_plus = spline.value(x + dx);
        Scalar f_minus = spline.value(x - dx);
        Scalar numerical_deriv = (f_plus - f_minus) / (2 * dx);

        REQUIRE_THAT(result.derivative, WithinRel(numerical_deriv, 1e-4));
    }
}
