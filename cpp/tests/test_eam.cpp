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

#include <cmath>
#include <filesystem>

#include "atomistica/atomistica.hpp"
#include "atomistica/potentials/eam/eam.hpp"

using namespace atomistica;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// Helper to find test data files
std::string find_eam_test_file(const std::string& filename) {
    // Try various paths relative to build directory
    std::vector<std::string> paths = {
        filename,
        "../tests/" + filename,
        "../../tests/" + filename,
        "../../../tests/" + filename,
        "../../../../tests/" + filename,
    };

    for (const auto& path : paths) {
        if (std::filesystem::exists(path)) {
            return path;
        }
    }

    throw std::runtime_error("Cannot find test file: " + filename);
}

TEST_CASE("TabulatedEAM file loading", "[eam]") {
    TabulatedEAM eam;

    SECTION("Load Au_u3.eam") {
        REQUIRE_NOTHROW(eam.load(find_eam_test_file("Au_u3.eam")));

        REQUIRE(eam.is_valid());

        // Check element info
        CHECK(eam.element_info().atomic_number == 79);
        CHECK_THAT(eam.element_info().mass, WithinRel(196.97, 0.01));
        CHECK_THAT(eam.element_info().lattice_constant, WithinRel(4.08, 0.01));

        // Check cutoff
        CHECK_THAT(eam.cutoff(), WithinRel(5.55, 0.01));
    }
}

TEST_CASE("TabulatedEAM spline functions", "[eam]") {
    TabulatedEAM eam;
    eam.load(find_eam_test_file("Au_u3.eam"));

    SECTION("Embedding function F(rho)") {
        // F(0) should be 0 (from file)
        auto F0 = eam.embedding(0.0);
        CHECK_THAT(F0.value, WithinAbs(0.0, 0.01));

        // F(rho) should be negative for positive rho (typical EAM)
        auto F1 = eam.embedding(0.1);
        CHECK(F1.value < 0.0);

        // Check that derivative exists
        CHECK(std::isfinite(F1.derivative));
    }

    SECTION("Density function rho(r)") {
        // rho(0) typically large
        auto rho0 = eam.density(0.1);
        CHECK(rho0.value > 0.0);

        // rho decreases with distance
        auto rho1 = eam.density(2.0);
        auto rho2 = eam.density(3.0);
        CHECK(rho1.value > rho2.value);

        // rho near cutoff should be small
        auto rho_cut = eam.density(5.5);
        CHECK_THAT(rho_cut.value, WithinAbs(0.0, 0.1));
    }

    SECTION("Pair potential phi(r)") {
        // phi should be positive (repulsive) at short distances
        auto phi1 = eam.pair_potential(2.0);
        CHECK(phi1.value > 0.0);

        // phi should decrease with distance
        auto phi2 = eam.pair_potential(3.0);
        CHECK(phi1.value > phi2.value);
    }
}

TEST_CASE("TabulatedEAM dimer energy", "[eam]") {
    TabulatedEAM eam;
    eam.load(find_eam_test_file("Au_u3.eam"));

    // Create a Au dimer at various separations
    std::vector<Scalar> distances = {2.5, 2.8, 3.0, 3.5, 4.0};

    for (Scalar d : distances) {
        AtomicSystem system(2);

        Mat3 cell;
        cell << 20.0, 0.0, 0.0,
                0.0, 20.0, 0.0,
                0.0, 0.0, 20.0;
        system.set_cell(cell);
        system.pbc() = {false, false, false};

        system.position(0) << 10.0, 10.0, 10.0;
        system.position(1) << 10.0 + d, 10.0, 10.0;
        system.atomic_numbers()(0) = 79;  // Au
        system.atomic_numbers()(1) = 79;

        NeighborList nl;
        nl.set_cutoff(eam.cutoff());
        nl.update(system);

        system.zero_forces();
        auto result = eam.compute(system, nl);

        // Energy should be finite
        CHECK(std::isfinite(result.energy));

        // For a dimer, forces should be equal and opposite
        Vec3 f0 = system.forces().col(0).matrix();
        Vec3 f1 = system.forces().col(1).matrix();
        CHECK_THAT(f0.x(), WithinAbs(-f1.x(), 1e-10));
        CHECK_THAT(f0.y(), WithinAbs(-f1.y(), 1e-10));
        CHECK_THAT(f0.z(), WithinAbs(-f1.z(), 1e-10));

        // Forces should be along the bond direction (x)
        CHECK_THAT(f0.y(), WithinAbs(0.0, 1e-10));
        CHECK_THAT(f0.z(), WithinAbs(0.0, 1e-10));
    }
}

TEST_CASE("TabulatedEAM numerical force test", "[eam]") {
    TabulatedEAM eam;
    eam.load(find_eam_test_file("Au_u3.eam"));

    // Create a small cluster
    AtomicSystem system(3);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.position(0) << 10.0, 10.0, 10.0;
    system.position(1) << 12.8, 10.0, 10.0;
    system.position(2) << 11.4, 12.4, 10.0;
    for (int i = 0; i < 3; ++i) {
        system.atomic_numbers()(i) = 79;
    }

    NeighborList nl;
    nl.set_cutoff(eam.cutoff());
    nl.update(system);

    system.zero_forces();
    eam.compute(system, nl);

    Array3X analytical_forces = system.forces();

    // Numerical derivative
    const Scalar dx = 1e-6;

    for (std::size_t i = 0; i < 3; ++i) {
        for (int d = 0; d < 3; ++d) {
            // Forward
            system.position(i)(d) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = eam.compute(system, nl, false, false);

            // Backward
            system.position(i)(d) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = eam.compute(system, nl, false, false);

            // Restore
            system.position(i)(d) += dx;
            system.positions_changed();

            // F = -dE/dr
            Scalar force_num = -(r_plus.energy - r_minus.energy) / (2 * dx);

            // Compare with analytical
            CHECK_THAT(analytical_forces(d, i), WithinRel(force_num, 1e-4));
        }
    }
}

TEST_CASE("TabulatedAlloyEAM file loading", "[eam][alloy]") {
    TabulatedAlloyEAM eam;

    SECTION("Load Cu_mishin1.eam.alloy") {
        REQUIRE_NOTHROW(eam.load(find_eam_test_file("Cu_mishin1.eam.alloy")));

        REQUIRE(eam.is_valid());

        // Check element count
        CHECK(eam.num_elements() == 1);

        // Check element info
        CHECK(eam.element_index("Cu") == 0);
        CHECK(eam.element_info(0).atomic_number == 29);

        // Check cutoff
        CHECK_THAT(eam.cutoff(), WithinRel(5.5, 0.02));
    }

    SECTION("Load Au-Grochola-JCP05.eam.alloy") {
        REQUIRE_NOTHROW(eam.load(find_eam_test_file("Au-Grochola-JCP05.eam.alloy")));

        REQUIRE(eam.is_valid());

        // Check element count
        CHECK(eam.num_elements() == 1);
        CHECK(eam.element_index("Au") == 0);
    }
}

TEST_CASE("TabulatedAlloyEAM dimer energy", "[eam][alloy]") {
    TabulatedAlloyEAM eam;
    eam.load(find_eam_test_file("Cu_mishin1.eam.alloy"));

    // Create a Cu dimer
    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.position(0) << 10.0, 10.0, 10.0;
    system.position(1) << 12.5, 10.0, 10.0;
    system.atomic_numbers()(0) = 29;  // Cu
    system.atomic_numbers()(1) = 29;

    NeighborList nl;
    nl.set_cutoff(eam.cutoff());
    nl.update(system);

    system.zero_forces();
    auto result = eam.compute(system, nl);

    // Energy should be finite
    CHECK(std::isfinite(result.energy));

    // Forces should be equal and opposite
    Vec3 f0 = system.forces().col(0).matrix();
    Vec3 f1 = system.forces().col(1).matrix();
    CHECK_THAT(f0.x(), WithinAbs(-f1.x(), 1e-10));
    CHECK_THAT(f0.y(), WithinAbs(-f1.y(), 1e-10));
    CHECK_THAT(f0.z(), WithinAbs(-f1.z(), 1e-10));
}

TEST_CASE("TabulatedAlloyEAM numerical force test", "[eam][alloy]") {
    TabulatedAlloyEAM eam;
    eam.load(find_eam_test_file("Cu_mishin1.eam.alloy"));

    // Create a small cluster
    AtomicSystem system(3);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.position(0) << 10.0, 10.0, 10.0;
    system.position(1) << 12.5, 10.0, 10.0;
    system.position(2) << 11.25, 12.2, 10.0;
    for (int i = 0; i < 3; ++i) {
        system.atomic_numbers()(i) = 29;  // Cu
    }

    NeighborList nl;
    nl.set_cutoff(eam.cutoff());
    nl.update(system);

    system.zero_forces();
    eam.compute(system, nl);

    Array3X analytical_forces = system.forces();

    // Numerical derivative
    const Scalar dx = 1e-6;

    for (std::size_t i = 0; i < 3; ++i) {
        for (int d = 0; d < 3; ++d) {
            // Forward
            system.position(i)(d) += dx;
            system.positions_changed();
            nl.update(system);
            auto r_plus = eam.compute(system, nl, false, false);

            // Backward
            system.position(i)(d) -= 2 * dx;
            system.positions_changed();
            nl.update(system);
            auto r_minus = eam.compute(system, nl, false, false);

            // Restore
            system.position(i)(d) += dx;
            system.positions_changed();

            // F = -dE/dr
            Scalar force_num = -(r_plus.energy - r_minus.energy) / (2 * dx);

            // Compare with analytical
            CHECK_THAT(analytical_forces(d, i), WithinRel(force_num, 1e-4));
        }
    }
}

TEST_CASE("TabulatedEAM FCC bulk", "[eam]") {
    TabulatedEAM eam;
    eam.load(find_eam_test_file("Au_u3.eam"));

    // Create 2x2x2 FCC Au
    const Scalar a = 4.08;  // Au lattice constant

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
                    system.position(idx) << (ix + basis[b](0)) * a,
                                            (iy + basis[b](1)) * a,
                                            (iz + basis[b](2)) * a;
                    system.atomic_numbers()(idx) = 79;
                    ++idx;
                }
            }
        }
    }

    NeighborList nl;
    nl.set_cutoff(eam.cutoff());
    nl.update(system);

    system.zero_forces();
    auto result = eam.compute(system, nl);

    SECTION("Energy is negative") {
        // Bound system should have negative energy
        REQUIRE(result.energy < 0);
    }

    SECTION("Energy per atom is reasonable") {
        Scalar energy_per_atom = result.energy / 32.0;
        // For FCC Au, cohesive energy is about -3.8 eV per atom
        // With our cutoff and periodic boundaries, it should be close
        REQUIRE(energy_per_atom < 0);
        REQUIRE(energy_per_atom > -10.0);  // Sanity check
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
