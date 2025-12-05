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

TEST_CASE("Tersoff parameter loading", "[Tersoff]") {
    Tersoff<false> pot;

    SECTION("Load Si-C parameters") {
        pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

        REQUIRE(pot.element_index(14) == 0);  // Si
        REQUIRE(pot.element_index(6) == 1);   // C
        REQUIRE(pot.element_index(1) == -1);  // H not defined
        REQUIRE(pot.num_elements() == 2);
        REQUIRE(pot.cutoff() > 2.0);
    }

    SECTION("Unknown parameter set throws") {
        REQUIRE_THROWS(pot.load_parameters("NonExistent"));
    }
}

TEST_CASE("Tersoff pair functions", "[Tersoff]") {
    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    // Si-Si pair type
    int ptype_si_si = pot.pair_type(0, 0);

    SECTION("Repulsive potential") {
        auto [VR, dVR] = pot.repulsive(ptype_si_si, 2.35);

        REQUIRE(VR > 0.0);  // Repulsive is positive
        REQUIRE(dVR < 0.0); // Derivative is negative (decays with distance)
    }

    SECTION("Attractive potential") {
        auto [VA, dVA] = pot.attractive(ptype_si_si, 2.35);

        REQUIRE(VA < 0.0);  // Attractive is negative
        REQUIRE(dVA > 0.0); // Derivative is positive (becomes less negative with distance)
    }

    SECTION("Angular function") {
        // At cos_theta = h, angular function has minimum
        auto [g, dg] = pot.angular_function(0, 0, 0, ptype_si_si, ptype_si_si, 0.0);
        REQUIRE(g > 0.0);
    }

    SECTION("Bond order at z=0") {
        auto [b, db] = pot.bond_order(0, ptype_si_si, 0.0);
        REQUIRE_THAT(b, WithinRel(1.0, 1e-6));
    }

    SECTION("Bond order decreases with z") {
        auto [b1, db1] = pot.bond_order(0, ptype_si_si, 1.0);
        auto [b2, db2] = pot.bond_order(0, ptype_si_si, 2.0);

        REQUIRE(b1 < 1.0);
        REQUIRE(b2 < b1);  // Higher coordination -> lower bond order
        REQUIRE(db1 < 0.0); // Derivative is negative
    }
}

TEST_CASE("Tersoff Si dimer", "[Tersoff]") {
    // Two Si atoms
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;  // Si
    system.atomic_numbers()(1) = 14;  // Si

    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());

    // Test at typical Si-Si bond length
    Scalar r_bond = 2.35;  // Angstrom

    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r_bond, 10.0, 10.0;

    nl.update(system);
    system.zero_forces();

    auto result = pot.compute(system, nl, true, true);

    SECTION("Energy is negative") {
        // For a dimer within cutoff, energy should be negative (bound)
        REQUIRE(result.energy < 0.0);
    }

    SECTION("Newton's third law") {
        Vec3 total_force = system.forces().col(0).matrix() + system.forces().col(1).matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }

    SECTION("Forces along bond axis") {
        // Forces should be along the x-axis (bond direction)
        REQUIRE_THAT(system.forces()(1, 0), WithinAbs(0.0, 1e-10));  // y component
        REQUIRE_THAT(system.forces()(2, 0), WithinAbs(0.0, 1e-10));  // z component
    }
}

TEST_CASE("Tersoff Si3 trimer", "[Tersoff]") {
    // Three Si atoms in a triangle to test angular forces
    AtomicSystem system(3);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;
    system.atomic_numbers()(1) = 14;
    system.atomic_numbers()(2) = 14;

    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());

    // Equilateral triangle
    Scalar r = 2.35;
    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r, 10.0, 10.0;
    system.positions().col(2) << 10.0 + 0.5*r, 10.0 + r*std::sqrt(3.0)/2.0, 10.0;

    nl.update(system);
    system.zero_forces();

    auto result = pot.compute(system, nl, true, true);

    SECTION("Energy is negative") {
        REQUIRE(result.energy < 0.0);
    }

    SECTION("Total force is zero") {
        Vec3 total_force = system.forces().col(0).matrix() +
                          system.forces().col(1).matrix() +
                          system.forces().col(2).matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }

    SECTION("Symmetric forces for equilateral triangle") {
        // By symmetry, force magnitudes should be equal
        Scalar f0 = system.forces().col(0).matrix().norm();
        Scalar f1 = system.forces().col(1).matrix().norm();
        Scalar f2 = system.forces().col(2).matrix().norm();

        REQUIRE_THAT(f0, WithinRel(f1, 1e-6));
        REQUIRE_THAT(f1, WithinRel(f2, 1e-6));
    }
}

TEST_CASE("Tersoff numerical force test", "[Tersoff]") {
    AtomicSystem system(3);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;
    system.atomic_numbers()(1) = 14;
    system.atomic_numbers()(2) = 14;

    // Asymmetric configuration
    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 12.3, 10.1, 10.0;
    system.positions().col(2) << 10.5, 12.2, 10.2;

    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    // Analytical forces
    system.zero_forces();
    auto result = pot.compute(system, nl, true, false);
    Array3X analytical_forces = system.forces();

    // Numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 3);

    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            system.positions()(dir, atom) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) += dx;
            system.positions_changed();

            numerical_forces(dir, atom) = -(r_plus.energy - r_minus.energy) / (2 * dx);
        }
    }

    // Compare
    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            if (std::abs(numerical_forces(dir, atom)) > 1e-8) {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinRel(numerical_forces(dir, atom), 1e-4));
            } else {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinAbs(numerical_forces(dir, atom), 1e-8));
            }
        }
    }
}

TEST_CASE("Tersoff SiC heteroatomic", "[Tersoff]") {
    // Si-C bond to test mixed parameters
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;  // Si
    system.atomic_numbers()(1) = 6;   // C

    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());

    // Typical Si-C bond length
    Scalar r_bond = 1.89;  // Angstrom

    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r_bond, 10.0, 10.0;

    nl.update(system);
    system.zero_forces();

    auto result = pot.compute(system, nl, true, true);

    SECTION("Energy is negative") {
        REQUIRE(result.energy < 0.0);
    }

    SECTION("Forces obey Newton's third law") {
        Vec3 total_force = system.forces().col(0).matrix() + system.forces().col(1).matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }
}

// ============================================================================
// Screened Tersoff Tests
// ============================================================================

TEST_CASE("Screened Tersoff parameter loading", "[TersoffScr]") {
    Tersoff<true> pot;

    SECTION("Load Si-C parameters") {
        pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

        REQUIRE(pot.element_index(14) == 0);  // Si
        REQUIRE(pot.element_index(6) == 1);   // C
        REQUIRE(pot.num_elements() == 2);
        // Screened cutoff should be larger than non-screened
        REQUIRE(pot.cutoff() > 3.0);
    }
}

TEST_CASE("Screened Tersoff Si dimer (unscreened region)", "[TersoffScr]") {
    // Two Si atoms at short distance - should be unscreened
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;  // Si
    system.atomic_numbers()(1) = 14;  // Si

    Tersoff<true> pot_scr;
    pot_scr.load_parameters("Tersoff_PRB_39_5566_Si_C");

    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot_scr.cutoff());

    // At short distance (well within inner cutoff), screening should not apply
    Scalar r_bond = 2.35;  // Angstrom

    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r_bond, 10.0, 10.0;

    nl.update(system);

    // Compute with non-screened
    system.zero_forces();
    auto result = pot.compute(system, nl, true, true);

    // Compute with screened (should give same result for isolated dimer)
    system.zero_forces();
    auto result_scr = pot_scr.compute(system, nl, true, true);

    SECTION("Energy is negative") {
        REQUIRE(result_scr.energy < 0.0);
    }

    SECTION("Dimer energy matches non-screened") {
        // For an isolated dimer with no screening atoms, energies should match
        REQUIRE_THAT(result_scr.energy, WithinRel(result.energy, 1e-6));
    }

    SECTION("Newton's third law") {
        Vec3 total_force = system.forces().col(0).matrix() + system.forces().col(1).matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("Screened Tersoff trimer screening", "[TersoffScr]") {
    // Three Si atoms arranged linearly - middle atom should screen the bond
    AtomicSystem system(3);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;
    system.atomic_numbers()(1) = 14;
    system.atomic_numbers()(2) = 14;

    // Linear arrangement: 0 -- 1 -- 2
    // Atom 1 should screen the 0-2 bond
    Scalar r = 2.5;
    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r, 10.0, 10.0;  // Middle atom
    system.positions().col(2) << 10.0 + 2*r, 10.0, 10.0;  // r_02 = 2r = 5.0

    Tersoff<true> pot_scr;
    pot_scr.load_parameters("Tersoff_PRB_39_5566_Si_C");

    Tersoff<false> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot_scr.cutoff());
    nl.update(system);

    system.zero_forces();
    auto result_scr = pot_scr.compute(system, nl, true, true);

    SECTION("Energy is finite") {
        REQUIRE(std::isfinite(result_scr.energy));
    }

    SECTION("Total force is zero") {
        Vec3 total_force = system.forces().col(0).matrix() +
                          system.forces().col(1).matrix() +
                          system.forces().col(2).matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("Screened Tersoff numerical force test (unscreened config)", "[TersoffScr]") {
    // Test numerical forces in a configuration where screening is minimal
    // Full screening derivative forces are not yet implemented
    AtomicSystem system(3);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;
    system.atomic_numbers()(1) = 14;
    system.atomic_numbers()(2) = 14;

    // Equilateral triangle at short distance - minimal screening effect
    Scalar r = 2.35;
    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r, 10.0, 10.0;
    system.positions().col(2) << 10.0 + 0.5*r, 10.0 + r*std::sqrt(3.0)/2.0, 10.0;

    Tersoff<true> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    // Analytical forces
    system.zero_forces();
    auto result = pot.compute(system, nl, true, false);
    Array3X analytical_forces = system.forces();

    // Numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 3);

    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            system.positions()(dir, atom) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) += dx;
            system.positions_changed();

            numerical_forces(dir, atom) = -(r_plus.energy - r_minus.energy) / (2 * dx);
        }
    }

    // Compare - this should match well for unscreened configurations
    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            if (std::abs(numerical_forces(dir, atom)) > 1e-8) {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinRel(numerical_forces(dir, atom), 1e-4));
            } else {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinAbs(numerical_forces(dir, atom), 1e-8));
            }
        }
    }
}

TEST_CASE("Screened Tersoff numerical force test (linear config with screening)", "[TersoffScr]") {
    // Test numerical forces in a linear configuration where screening is active
    // The middle atom (1) screens the 0-2 bond
    AtomicSystem system(3);

    Mat3 cell;
    cell << 30.0, 0.0, 0.0,
            0.0, 30.0, 0.0,
            0.0, 0.0, 30.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;
    system.atomic_numbers()(1) = 14;
    system.atomic_numbers()(2) = 14;

    // Linear arrangement: 0 -- 1 -- 2
    // Distances: r_01 = 2.5, r_12 = 2.5, r_02 = 5.0
    // Atom 1 is exactly between 0 and 2, so it should provide maximum screening
    Scalar r = 2.5;
    system.positions().col(0) << 15.0, 15.0, 15.0;
    system.positions().col(1) << 15.0 + r, 15.0, 15.0;  // Middle atom
    system.positions().col(2) << 15.0 + 2*r, 15.0, 15.0;

    Tersoff<true> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    // Analytical forces
    system.zero_forces();
    auto result = pot.compute(system, nl, true, false);
    Array3X analytical_forces = system.forces();

    // Numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 3);

    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            system.positions()(dir, atom) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) += dx;
            system.positions_changed();

            numerical_forces(dir, atom) = -(r_plus.energy - r_minus.energy) / (2 * dx);
        }
    }

    // Compare analytical and numerical forces
    // This tests that screening force derivatives are correctly implemented
    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            if (std::abs(numerical_forces(dir, atom)) > 1e-8) {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinRel(numerical_forces(dir, atom), 1e-4));
            } else {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinAbs(numerical_forces(dir, atom), 1e-8));
            }
        }
    }
}

TEST_CASE("Screened Tersoff numerical force test (off-axis screener)", "[TersoffScr]") {
    // Test with a screener that's not exactly on the bond axis
    // This provides a more thorough test of the screening derivatives
    AtomicSystem system(3);

    Mat3 cell;
    cell << 30.0, 0.0, 0.0,
            0.0, 30.0, 0.0,
            0.0, 0.0, 30.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;
    system.atomic_numbers()(1) = 14;
    system.atomic_numbers()(2) = 14;

    // Asymmetric configuration: atom 2 slightly off the 0-1 axis
    // r_01 = 4.0, atom 2 is near the midpoint but offset in y
    system.positions().col(0) << 15.0, 15.0, 15.0;
    system.positions().col(1) << 19.0, 15.0, 15.0;  // r_01 = 4.0
    system.positions().col(2) << 17.0, 15.5, 15.0;  // Near midpoint, offset by 0.5 in y

    Tersoff<true> pot;
    pot.load_parameters("Tersoff_PRB_39_5566_Si_C");

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    // Analytical forces
    system.zero_forces();
    auto result = pot.compute(system, nl, true, false);
    Array3X analytical_forces = system.forces();

    // Numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 3);

    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            system.positions()(dir, atom) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = pot.compute(system, nl, false, false);

            system.positions()(dir, atom) += dx;
            system.positions_changed();

            numerical_forces(dir, atom) = -(r_plus.energy - r_minus.energy) / (2 * dx);
        }
    }

    // Compare analytical and numerical forces
    for (int atom = 0; atom < 3; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            if (std::abs(numerical_forces(dir, atom)) > 1e-8) {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinRel(numerical_forces(dir, atom), 1e-4));
            } else {
                REQUIRE_THAT(analytical_forces(dir, atom),
                            WithinAbs(numerical_forces(dir, atom), 1e-8));
            }
        }
    }
}
