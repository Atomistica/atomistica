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

TEST_CASE("LJ dimer energy", "[LJ]") {
    // Two-atom system
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    // Argon parameters
    const Scalar epsilon = 0.0103;  // eV
    const Scalar sigma = 3.40;      // Angstrom
    const Scalar cutoff = 10.0;

    system.atomic_numbers()(0) = 18;  // Argon
    system.atomic_numbers()(1) = 18;

    LJCut lj(18, epsilon, sigma, cutoff);

    NeighborList nl;
    nl.set_cutoff(lj.cutoff());

    SECTION("At equilibrium distance") {
        // Equilibrium distance: r0 = 2^(1/6) * sigma
        Scalar r0 = std::pow(2.0, 1.0/6.0) * sigma;

        system.positions().col(0) << 10.0, 10.0, 10.0;
        system.positions().col(1) << 10.0 + r0, 10.0, 10.0;

        nl.update(system);
        system.zero_forces();

        auto result = lj.compute(system, nl);

        // Energy at equilibrium: -epsilon
        REQUIRE_THAT(result.energy, WithinRel(-epsilon, 1e-6));

        // Forces should be zero at equilibrium
        REQUIRE_THAT(system.forces().col(0).matrix().norm(), WithinAbs(0.0, 1e-10));
        REQUIRE_THAT(system.forces().col(1).matrix().norm(), WithinAbs(0.0, 1e-10));
    }

    SECTION("At sigma distance") {
        // At r = sigma: V = 0
        system.positions().col(0) << 10.0, 10.0, 10.0;
        system.positions().col(1) << 10.0 + sigma, 10.0, 10.0;

        nl.update(system);
        system.zero_forces();

        auto result = lj.compute(system, nl);

        REQUIRE_THAT(result.energy, WithinAbs(0.0, 1e-10));
    }

    SECTION("Force direction") {
        // At distance less than r0: repulsive
        Scalar r = sigma * 0.9;

        system.positions().col(0) << 10.0, 10.0, 10.0;
        system.positions().col(1) << 10.0 + r, 10.0, 10.0;

        nl.update(system);
        system.zero_forces();

        lj.compute(system, nl);

        // Atom 0 should be pushed in -x direction
        REQUIRE(system.forces()(0, 0) < 0);
        // Atom 1 should be pushed in +x direction
        REQUIRE(system.forces()(0, 1) > 0);
    }

    SECTION("Newton's third law") {
        system.positions().col(0) << 10.0, 10.0, 10.0;
        system.positions().col(1) << 14.0, 10.0, 10.0;

        nl.update(system);
        system.zero_forces();

        lj.compute(system, nl);

        // F_0 + F_1 = 0
        Vec3 total_force = system.forces().col(0).matrix() + system.forces().col(1).matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("LJ numerical force test", "[LJ]") {
    // Verify forces match numerical derivative
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    const Scalar epsilon = 0.0103;
    const Scalar sigma = 3.40;
    const Scalar cutoff = 10.0;

    system.atomic_numbers()(0) = 18;
    system.atomic_numbers()(1) = 18;

    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 14.5, 10.3, 10.7;

    LJCut lj(18, epsilon, sigma, cutoff);

    NeighborList nl;
    nl.set_cutoff(lj.cutoff());
    nl.update(system);

    // Compute analytical forces
    system.zero_forces();
    auto result = lj.compute(system, nl);

    Array3X analytical_forces = system.forces();

    // Compute numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 2);

    for (int atom = 0; atom < 2; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            // Forward
            system.positions()(dir, atom) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = lj.compute(system, nl, false, false);

            // Backward
            system.positions()(dir, atom) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = lj.compute(system, nl, false, false);

            // Restore
            system.positions()(dir, atom) += dx;
            system.positions_changed();

            // F = -dE/dr
            numerical_forces(dir, atom) = -(r_plus.energy - r_minus.energy) / (2 * dx);
        }
    }

    // Compare
    for (int atom = 0; atom < 2; ++atom) {
        for (int dir = 0; dir < 3; ++dir) {
            REQUIRE_THAT(analytical_forces(dir, atom),
                        WithinRel(numerical_forces(dir, atom), 1e-4));
        }
    }
}

TEST_CASE("LJ shifted vs unshifted", "[LJ]") {
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    const Scalar epsilon = 0.0103;
    const Scalar sigma = 3.40;
    const Scalar cutoff = 10.0;

    system.atomic_numbers()(0) = 18;
    system.atomic_numbers()(1) = 18;

    LJCut lj_unshifted(18, epsilon, sigma, cutoff);
    LJCutShift lj_shifted(18, epsilon, sigma, cutoff);

    NeighborList nl;
    nl.set_cutoff(cutoff);

    SECTION("Energy difference is constant") {
        // The difference between shifted and unshifted should be constant
        // for all distances < cutoff

        std::vector<Scalar> distances = {3.5, 4.0, 5.0, 6.0, 8.0, 9.0};
        Scalar diff = 0.0;
        bool first = true;

        for (Scalar r : distances) {
            system.positions().col(0) << 10.0, 10.0, 10.0;
            system.positions().col(1) << 10.0 + r, 10.0, 10.0;

            nl.update(system);

            auto result_unshifted = lj_unshifted.compute(system, nl, false, false);
            auto result_shifted = lj_shifted.compute(system, nl, false, false);

            Scalar current_diff = result_unshifted.energy - result_shifted.energy;

            if (first) {
                diff = current_diff;
                first = false;
            } else {
                REQUIRE_THAT(current_diff, WithinRel(diff, 1e-10));
            }
        }
    }

    SECTION("Shifted energy goes to zero at cutoff") {
        // Just below cutoff
        system.positions().col(0) << 10.0, 10.0, 10.0;
        system.positions().col(1) << 10.0 + cutoff - 0.01, 10.0, 10.0;

        nl.update(system);
        auto result = lj_shifted.compute(system, nl, false, false);

        // Energy should be very small
        REQUIRE_THAT(result.energy, WithinAbs(0.0, 1e-6));
    }

    SECTION("Forces are identical") {
        // Forces should be identical for shifted and unshifted
        system.positions().col(0) << 10.0, 10.0, 10.0;
        system.positions().col(1) << 14.0, 10.0, 10.0;

        nl.update(system);

        system.zero_forces();
        lj_unshifted.compute(system, nl);
        Array3X forces_unshifted = system.forces();

        system.zero_forces();
        lj_shifted.compute(system, nl);
        Array3X forces_shifted = system.forces();

        REQUIRE(forces_unshifted.isApprox(forces_shifted));
    }
}

TEST_CASE("LJ FCC bulk", "[LJ]") {
    // Create 2x2x2 FCC Argon
    const Scalar a = 5.26;  // Argon lattice constant at 0K
    const Scalar epsilon = 0.0103;
    const Scalar sigma = 3.40;
    const Scalar cutoff = 2.5 * sigma;

    AtomicSystem system(32);

    Mat3 cell;
    cell << 2*a, 0.0, 0.0,
            0.0, 2*a, 0.0,
            0.0, 0.0, 2*a;
    system.set_cell(cell);

    Vec3 basis[4] = {
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
        {0.5, 0.0, 0.5},
        {0.0, 0.5, 0.5}
    };

    int idx = 0;
    for (int iz = 0; iz < 2; ++iz) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int ix = 0; ix < 2; ++ix) {
                for (int b = 0; b < 4; ++b) {
                    system.positions().col(idx) << (ix + basis[b](0)) * a,
                                                   (iy + basis[b](1)) * a,
                                                   (iz + basis[b](2)) * a;
                    system.atomic_numbers()(idx) = 18;
                    ++idx;
                }
            }
        }
    }

    LJCut lj(18, epsilon, sigma, cutoff);

    NeighborList nl;
    nl.set_cutoff(lj.cutoff());
    nl.update(system);

    system.zero_forces();
    auto result = lj.compute(system, nl);

    SECTION("Energy is negative") {
        // Bound system should have negative energy
        REQUIRE(result.energy < 0);
    }

    SECTION("Energy per atom is reasonable") {
        Scalar energy_per_atom = result.energy / 32.0;
        // For FCC LJ, cohesive energy is about -8.6 * epsilon
        // (with infinite cutoff). With our cutoff, it's less negative.
        REQUIRE(energy_per_atom > -10 * epsilon);
        REQUIRE(energy_per_atom < 0);
    }

    SECTION("Total force is zero") {
        // In a perfect crystal, total force should be zero
        Vec3 total_force = system.forces().rowwise().sum().matrix();
        REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-10));
    }

    SECTION("Virial is symmetric") {
        REQUIRE(result.virial.isApprox(result.virial.transpose()));
    }
}
