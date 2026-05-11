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

TEST_CASE("Juslin parameter loading - WCH", "[Juslin]") {
    Juslin<false> pot;

    SECTION("Load W-C-H parameters") {
        pot.load_parameters("Juslin_JAP_98_123520_WCH");

        REQUIRE(pot.element_index(74) == 0);  // W
        REQUIRE(pot.element_index(6)  == 1);  // C
        REQUIRE(pot.element_index(1)  == 2);  // H
        REQUIRE(pot.element_index(14) == -1); // Si not defined
        REQUIRE(pot.num_elements() == 3);
        REQUIRE(pot.cutoff() > 1.0);
    }

    SECTION("Unknown parameter set throws") {
        REQUIRE_THROWS(pot.load_parameters("NonExistent"));
    }
}

TEST_CASE("Juslin pair functions - WCH", "[Juslin]") {
    Juslin<false> pot;
    pot.load_parameters("Juslin_JAP_98_123520_WCH");

    SECTION("W-W pair functions") {
        int ptype_ww = pot.pair_type(0, 0);

        auto [VR, dVR] = pot.repulsive(ptype_ww, 2.5);
        REQUIRE(VR > 0.0);
        REQUIRE(dVR < 0.0);

        auto [VA, dVA] = pot.attractive(ptype_ww, 2.5);
        REQUIRE(VA < 0.0);
        REQUIRE(dVA > 0.0);
    }

    SECTION("C-H pair functions") {
        int ptype_ch = pot.pair_type(1, 2);  // C=1, H=2

        auto [VR, dVR] = pot.repulsive(ptype_ch, 1.2);
        REQUIRE(VR > 0.0);

        auto [VA, dVA] = pot.attractive(ptype_ch, 1.2);
        REQUIRE(VA < 0.0);
    }

    SECTION("W-C is inactive (r2=0)") {
        int ptype_wc = pot.pair_type(0, 1);  // W=0, C=1
        REQUIRE_THAT(pot.pair_cutoff(ptype_wc), WithinAbs(0.0, 1e-15));
    }
}

TEST_CASE("Juslin triplet distance function", "[Juslin]") {
    Juslin<false> pot;
    pot.load_parameters("Juslin_JAP_98_123520_WCH");

    SECTION("W central, W-W triplet: alpha=0.45876") {
        // eli=W=0, elj=W=0, elk=W=0
        auto [h, dh_drik, dh_drij] = pot.distance_function(0, 0, 0, 0, 0, 2.5, 2.4);
        // dr = 2.5 - 2.4 = 0.1, h = omega*exp(alpha*dr) = 1*exp(0.45876*0.1) > 1
        REQUIRE(h > 1.0);
        REQUIRE(dh_drik < 0.0);  // dh/dr_ik = -alpha*h < 0
        REQUIRE(dh_drij > 0.0);  // dh/dr_ij = +alpha*h > 0
    }

    SECTION("W central, W-C triplet: alpha=0, h=omega=1 constant") {
        // eli=W=0, elj=W=0, elk=C=1
        auto [h, dh_drik, dh_drij] = pot.distance_function(0, 0, 1, 0, 1, 2.5, 2.4);
        REQUIRE_THAT(h, WithinRel(1.0, 1e-10));
        REQUIRE_THAT(dh_drik, WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(dh_drij, WithinAbs(0.0, 1e-10));
    }

    SECTION("C central, H-C triplet: alpha=4, omega=2.94586") {
        // eli=C=1, elj=H=2, elk=C=1
        auto [h, dh_drik, dh_drij] = pot.distance_function(1, 2, 1, 0, 0, 1.5, 1.4);
        // dr = 0.1, h = 2.94586 * exp(4 * 0.1) = 2.94586 * exp(0.4)
        Scalar expected = 2.94586 * std::exp(4.0 * 0.1);
        REQUIRE_THAT(h, WithinRel(expected, 1e-6));
    }
}

// Helper: build a simple W dimer
static std::pair<AtomicSystem, NeighborList> make_w_dimer(Scalar r = 2.5) {
    AtomicSystem sys(2);
    Mat3 cell = Mat3::Identity() * 20.0;
    sys.set_cell(cell);
    sys.pbc() = {false, false, false};
    sys.positions().col(0) << 0.0, 0.0, 0.0;
    sys.positions().col(1) << r, 0.0, 0.0;
    sys.atomic_numbers()(0) = 74;  // W
    sys.atomic_numbers()(1) = 74;

    NeighborList nl;
    nl.set_cutoff(5.0);
    nl.update(sys);
    return {sys, nl};
}

TEST_CASE("Juslin W-dimer energy", "[Juslin]") {
    Juslin<false> pot;
    pot.load_parameters("Juslin_JAP_98_123520_WCH");

    SECTION("Energy is finite for W dimer at equilibrium distance") {
        auto [sys, nl] = make_w_dimer(2.34095);  // r0 for W-W
        auto res = pot.compute(sys, nl, false, false);
        REQUIRE(std::isfinite(res.energy));
        REQUIRE(res.energy < 0.0);  // Should be negative (bound)
    }

    SECTION("Energy increases at shorter distance") {
        auto [sys1, nl1] = make_w_dimer(2.1);
        auto [sys2, nl2] = make_w_dimer(2.34095);
        auto r1 = pot.compute(sys1, nl1, false, false);
        auto r2 = pot.compute(sys2, nl2, false, false);
        REQUIRE(r1.energy > r2.energy);
    }
}

TEST_CASE("Juslin force-energy consistency", "[Juslin]") {
    Juslin<false> pot;
    pot.load_parameters("Juslin_JAP_98_123520_WCH");

    const Scalar dx = 1e-5;
    const Scalar tol = 1e-3;

    auto [sys, nl] = make_w_dimer(2.5);

    // Compute analytical forces
    auto res = pot.compute(sys, nl, true, false);
    Scalar fx_analytical = static_cast<Scalar>(sys.forces()(0, 0));

    // Numerical force via finite differences
    sys.positions()(0, 0) += dx;
    sys.positions_changed();
    nl.update(sys);
    Scalar E_plus = pot.compute(sys, nl, false, false).energy;

    sys.positions()(0, 0) -= 2 * dx;
    sys.positions_changed();
    nl.update(sys);
    Scalar E_minus = pot.compute(sys, nl, false, false).energy;

    sys.positions()(0, 0) += dx;  // restore

    Scalar fx_numerical = -(E_plus - E_minus) / (2 * dx);

    REQUIRE_THAT(fx_analytical, WithinRel(fx_numerical, tol));
}

TEST_CASE("JuslinScr parameter loading", "[Juslin]") {
    Juslin<true> pot;
    pot.load_parameters("Juslin_JAP_98_123520_WCH");

    REQUIRE(pot.element_index(74) == 0);
    REQUIRE(pot.element_index(6)  == 1);
    REQUIRE(pot.element_index(1)  == 2);
    REQUIRE(pot.cutoff() > 1.0);
}
