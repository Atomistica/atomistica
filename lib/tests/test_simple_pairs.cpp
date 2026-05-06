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

#include <atomistica/atomistica.hpp>

#include <cmath>

using namespace atomistica;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// Helper: build two-atom system
static std::pair<AtomicSystem, NeighborList> make_dimer(int Z, Scalar r = 1.0,
                                                         Scalar cutoff = 3.0) {
    AtomicSystem sys(2);
    Mat3 cell = Mat3::Identity() * 20.0;
    sys.set_cell(cell);
    sys.pbc() = {false, false, false};
    sys.positions().col(0) << 0.0, 0.0, 0.0;
    sys.positions().col(1) << r, 0.0, 0.0;
    sys.atomic_numbers()(0) = Z;
    sys.atomic_numbers()(1) = Z;

    NeighborList nl;
    nl.set_cutoff(cutoff);
    nl.update(sys);
    return {sys, nl};
}

// Numerical force via central differences
template<typename PotType>
static Scalar numerical_fx(PotType& pot, AtomicSystem sys, NeighborList nl) {
    const Scalar dx = 1e-5;
    sys.positions()(0, 0) += dx;
    sys.positions_changed();
    nl.update(sys);
    Scalar Ep = pot.compute(sys, nl, false, false).energy;

    sys.positions()(0, 0) -= 2 * dx;
    sys.positions_changed();
    nl.update(sys);
    Scalar Em = pot.compute(sys, nl, false, false).energy;
    return -(Ep - Em) / (2 * dx);
}

// ============================================================================
// BornMayer tests
// ============================================================================

TEST_CASE("BornMayer potential", "[BornMayer]") {
    BornMayer pot(2.0, 0.5, 3.0);  // A=2, rho=0.5, cutoff=3.0

    SECTION("Energy at r < cutoff is positive (repulsive)") {
        auto [sys, nl] = make_dimer(2, 1.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE(res.energy > 0.0);
    }

    SECTION("Energy is zero at cutoff") {
        // shift = A*exp(-cutoff/rho), so V(cutoff)=A*exp(-c/rho)-shift=0
        auto [sys, nl] = make_dimer(2, 3.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE_THAT(res.energy, WithinAbs(0.0, 1e-10));
    }

    SECTION("Force-energy consistency") {
        auto [sys, nl] = make_dimer(2, 1.5);
        auto res = pot.compute(sys, nl, true, false);
        Scalar fx_ana = static_cast<Scalar>(sys.forces()(0, 0));
        Scalar fx_num = numerical_fx(pot, sys, nl);
        REQUIRE_THAT(fx_ana, WithinRel(fx_num, 1e-3));
    }

    SECTION("Element filtering: Z=0 matches all") {
        BornMayer pot_any(2.0, 0.5, 3.0, 0, 0);  // wildcard
        auto [sys1, nl1] = make_dimer(6, 1.0);  // Carbon
        auto [sys2, nl2] = make_dimer(14, 1.0); // Silicon
        Scalar E1 = pot_any.compute(sys1, nl1, false, false).energy;
        Scalar E2 = pot_any.compute(sys2, nl2, false, false).energy;
        REQUIRE(std::abs(E1) > 0.0);
        REQUIRE(std::abs(E2) > 0.0);
    }

    SECTION("Element filtering: wrong Z gives zero energy") {
        BornMayer pot_he(2.0, 0.5, 3.0, 2, 2);  // He only
        auto [sys, nl] = make_dimer(6, 1.0);     // Carbon, not He
        Scalar E = pot_he.compute(sys, nl, false, false).energy;
        REQUIRE_THAT(E, WithinAbs(0.0, 1e-15));
    }
}

// ============================================================================
// Harmonic tests
// ============================================================================

TEST_CASE("Harmonic potential", "[Harmonic]") {
    Harmonic pot(1.0, 1.0, 2.0);  // k=1, r0=1, cutoff=2

    SECTION("Energy is zero at equilibrium") {
        auto [sys, nl] = make_dimer(2, 1.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE_THAT(res.energy, WithinAbs(0.0, 1e-10));
    }

    SECTION("Energy is positive away from equilibrium") {
        auto [sys, nl] = make_dimer(2, 1.5);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE(res.energy > 0.0);
        REQUIRE_THAT(res.energy, WithinAbs(0.5 * 1.0 * 0.25, 1e-10));  // 0.5*k*(dr)^2
    }

    SECTION("Force points toward equilibrium") {
        auto [sys, nl] = make_dimer(2, 1.5);  // r > r0, atom 0 pulled in +x
        pot.compute(sys, nl, true, false);
        REQUIRE(static_cast<Scalar>(sys.forces()(0, 0)) > 0.0);
        REQUIRE(static_cast<Scalar>(sys.forces()(0, 1)) < 0.0);
    }

    SECTION("Force-energy consistency") {
        auto [sys, nl] = make_dimer(2, 1.3);
        auto res = pot.compute(sys, nl, true, false);
        Scalar fx_ana = static_cast<Scalar>(sys.forces()(0, 0));
        Scalar fx_num = numerical_fx(pot, sys, nl);
        REQUIRE_THAT(fx_ana, WithinRel(fx_num, 1e-3));
    }

    SECTION("Shift makes V(cutoff) = 0") {
        Harmonic pot_shift(1.0, 1.0, 1.5, true);
        auto [sys, nl] = make_dimer(2, 1.5);
        Scalar E = pot_shift.compute(sys, nl, false, false).energy;
        REQUIRE_THAT(E, WithinAbs(0.0, 1e-10));
    }
}

// ============================================================================
// DoubleHarmonic tests
// ============================================================================

TEST_CASE("DoubleHarmonic potential", "[DoubleHarmonic]") {
    const Scalar sqrt2 = std::sqrt(2.0);
    DoubleHarmonic pot(1.0, 1.0, 1.0, sqrt2, 1.6);

    SECTION("Energy at r1 is zero") {
        auto [sys, nl] = make_dimer(2, 1.0, 2.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE_THAT(res.energy, WithinAbs(0.0, 1e-10));
    }

    SECTION("Energy at r2 is zero") {
        auto [sys, nl] = make_dimer(2, sqrt2, 2.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE_THAT(res.energy, WithinAbs(0.0, 1e-10));
    }

    SECTION("Energy is positive between r1 and r2") {
        Scalar rm = 0.5 * (1.0 + sqrt2);
        auto [sys, nl] = make_dimer(2, rm, 2.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE(res.energy >= 0.0);
    }

    SECTION("Force-energy consistency, r < rm") {
        auto [sys, nl] = make_dimer(2, 1.1, 2.0);
        auto res = pot.compute(sys, nl, true, false);
        Scalar fx_ana = static_cast<Scalar>(sys.forces()(0, 0));
        Scalar fx_num = numerical_fx(pot, sys, nl);
        REQUIRE_THAT(fx_ana, WithinRel(fx_num, 1e-3));
    }

    SECTION("Force-energy consistency, r > rm") {
        auto [sys, nl] = make_dimer(2, 1.3, 2.0);
        auto res = pot.compute(sys, nl, true, false);
        Scalar fx_ana = static_cast<Scalar>(sys.forces()(0, 0));
        Scalar fx_num = numerical_fx(pot, sys, nl);
        REQUIRE_THAT(fx_ana, WithinRel(fx_num, 1e-3));
    }
}

// ============================================================================
// R6 tests
// ============================================================================

TEST_CASE("R6 potential", "[R6]") {
    R6 pot(-1.0, 0.0, 5.0);  // A=-1 (attractive), r0=0, cutoff=5

    SECTION("Energy is negative (attractive)") {
        auto [sys, nl] = make_dimer(2, 2.0, 6.0);
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE(res.energy < 0.0);
        REQUIRE_THAT(res.energy, WithinAbs(-1.0 / std::pow(2.0, 6), 1e-10));
    }

    SECTION("Energy decays with distance") {
        auto [sys1, nl1] = make_dimer(2, 1.5, 6.0);
        auto [sys2, nl2] = make_dimer(2, 3.0, 6.0);
        Scalar E1 = pot.compute(sys1, nl1, false, false).energy;
        Scalar E2 = pot.compute(sys2, nl2, false, false).energy;
        REQUIRE(std::abs(E1) > std::abs(E2));
    }

    SECTION("Force-energy consistency") {
        R6 pot_rep(1.0, 0.5, 5.0);  // repulsive
        auto [sys, nl] = make_dimer(2, 2.0, 6.0);
        auto res = pot_rep.compute(sys, nl, true, false);
        Scalar fx_ana = static_cast<Scalar>(sys.forces()(0, 0));
        Scalar fx_num = numerical_fx(pot_rep, sys, nl);
        REQUIRE_THAT(fx_ana, WithinRel(fx_num, 1e-3));
    }

    SECTION("r0 offset shifts the singularity") {
        R6 pot_off(1.0, 1.0, 5.0);  // V(r) = 1/(r0+r)^6 = 1/(1+r)^6
        // At r=0.5: V = 1/(1+0.5)^6 = 1/1.5^6
        auto [sys, nl] = make_dimer(2, 0.5, 6.0);
        auto res = pot_off.compute(sys, nl, false, false);
        REQUIRE(std::isfinite(res.energy));
        REQUIRE_THAT(res.energy, WithinRel(1.0 / std::pow(1.5, 6), 1e-6));
    }
}
