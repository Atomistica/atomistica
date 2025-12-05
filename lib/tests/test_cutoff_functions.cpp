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

#include <atomistica/math/cutoff_functions.hpp>

#include <cmath>

using namespace atomistica;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("TrigOffCutoff boundary values", "[Cutoff]") {
    TrigOffCutoff fc(2.0, 3.0);

    SECTION("Below r1") {
        auto result = fc(1.5);
        REQUIRE_THAT(result.fc, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("At r1") {
        auto result = fc(2.0);
        REQUIRE_THAT(result.fc, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("At midpoint") {
        auto result = fc(2.5);
        REQUIRE_THAT(result.fc, WithinRel(0.5, 1e-10));
    }

    SECTION("At r2") {
        auto result = fc(3.0);
        REQUIRE_THAT(result.fc, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("Above r2") {
        auto result = fc(4.0);
        REQUIRE_THAT(result.fc, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("TrigOnCutoff boundary values", "[Cutoff]") {
    TrigOnCutoff fc(2.0, 3.0);

    SECTION("Below r1") {
        auto result = fc(1.5);
        REQUIRE_THAT(result.fc, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("At r1") {
        auto result = fc(2.0);
        REQUIRE_THAT(result.fc, WithinAbs(0.0, 1e-10));
    }

    SECTION("At midpoint") {
        auto result = fc(2.5);
        REQUIRE_THAT(result.fc, WithinRel(0.5, 1e-10));
    }

    SECTION("At r2") {
        auto result = fc(3.0);
        REQUIRE_THAT(result.fc, WithinRel(1.0, 1e-10));
    }

    SECTION("Above r2") {
        auto result = fc(4.0);
        REQUIRE_THAT(result.fc, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("ExpCutoff boundary values", "[Cutoff]") {
    ExpCutoff fc(2.0, 3.0);

    SECTION("Below r1") {
        auto result = fc(1.5);
        REQUIRE_THAT(result.fc, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("At r1") {
        auto result = fc(2.0);
        REQUIRE_THAT(result.fc, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("At r2") {
        auto result = fc(3.0);
        REQUIRE_THAT(result.fc, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("Above r2") {
        auto result = fc(4.0);
        REQUIRE_THAT(result.fc, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(result.dfc, WithinAbs(0.0, 1e-10));
    }

    SECTION("Monotonically decreasing") {
        Scalar prev = 1.0;
        for (Scalar r = 2.05; r < 3.0; r += 0.1) {
            auto result = fc(r);
            REQUIRE(result.fc < prev);
            REQUIRE(result.fc > 0.0);
            REQUIRE(result.dfc < 0.0);  // Negative derivative
            prev = result.fc;
        }
    }
}

TEST_CASE("Cutoff numerical derivative", "[Cutoff]") {
    const Scalar dx = 1e-6;

    SECTION("TrigOffCutoff") {
        TrigOffCutoff fc(2.0, 3.0);

        for (Scalar r = 2.1; r < 2.9; r += 0.2) {
            auto result = fc(r);
            auto fp = fc(r + dx);
            auto fm = fc(r - dx);
            Scalar numerical_deriv = (fp.fc - fm.fc) / (2 * dx);
            REQUIRE_THAT(result.dfc, WithinRel(numerical_deriv, 1e-5));
        }
    }

    SECTION("TrigOnCutoff") {
        TrigOnCutoff fc(2.0, 3.0);

        for (Scalar r = 2.1; r < 2.9; r += 0.2) {
            auto result = fc(r);
            auto fp = fc(r + dx);
            auto fm = fc(r - dx);
            Scalar numerical_deriv = (fp.fc - fm.fc) / (2 * dx);
            REQUIRE_THAT(result.dfc, WithinRel(numerical_deriv, 1e-5));
        }
    }

    SECTION("ExpCutoff") {
        ExpCutoff fc(2.0, 3.0);

        for (Scalar r = 2.1; r < 2.9; r += 0.2) {
            auto result = fc(r);
            auto fp = fc(r + dx);
            auto fm = fc(r - dx);
            Scalar numerical_deriv = (fp.fc - fm.fc) / (2 * dx);
            REQUIRE_THAT(result.dfc, WithinRel(numerical_deriv, 1e-4));
        }
    }
}

TEST_CASE("BOPCutoff wrapper", "[Cutoff]") {
    SECTION("TrigOff type") {
        BOPCutoff fc(BOPCutoff::Type::TrigOff, 2.0, 3.0);

        auto result = fc(2.5);
        REQUIRE_THAT(result.fc, WithinRel(0.5, 1e-10));
        REQUIRE_THAT(fc.r1(), WithinRel(2.0, 1e-10));
        REQUIRE_THAT(fc.r2(), WithinRel(3.0, 1e-10));
    }

    SECTION("Exp type") {
        BOPCutoff fc(BOPCutoff::Type::Exp, 2.0, 3.0);

        auto result_low = fc(1.5);
        auto result_high = fc(3.5);
        REQUIRE_THAT(result_low.fc, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(result_high.fc, WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("ExpCutoff C2 continuity", "[Cutoff]") {
    // Verify the exponential cutoff has continuous second derivative
    ExpCutoff fc(2.0, 3.0);
    const Scalar dx = 1e-5;

    SECTION("Second derivative at r2") {
        // At r2, the second derivative should also be zero
        Scalar r = 3.0 - dx;
        auto fp = fc(r + dx);
        auto f0 = fc(r);
        auto fm = fc(r - dx);

        // Second derivative approximation
        Scalar d2f = (fp.fc - 2*f0.fc + fm.fc) / (dx * dx);

        // Should be close to zero at boundary
        REQUIRE_THAT(d2f, WithinAbs(0.0, 1e-2));
    }
}
