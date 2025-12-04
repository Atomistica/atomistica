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
#include <atomistica/atomistica.hpp>

using namespace atomistica;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("Kumagai parameter loading", "[Kumagai]") {
    Kumagai<false> pot;

    SECTION("Load Kumagai Si parameters") {
        pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

        REQUIRE(pot.num_elements() == 1);
        REQUIRE(pot.element_index(14) == 0);  // Si
        REQUIRE(pot.element_index(6) == -1);  // C not defined
        REQUIRE(pot.cutoff() > 0.0);
    }

    SECTION("Unknown parameter set throws") {
        REQUIRE_THROWS(pot.load_parameters("NonExistent"));
    }
}

TEST_CASE("Kumagai pair functions", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

    // Si-Si pair type
    int ptype = pot.pair_type(0, 0);

    SECTION("Repulsive function") {
        Scalar r = 2.35;  // Near Si equilibrium distance
        auto [V_R, dV_R] = pot.repulsive(ptype, r);

        // V_R should be positive (repulsive)
        REQUIRE(V_R > 0.0);
        // dV_R should be negative (energy decreases with distance)
        REQUIRE(dV_R < 0.0);

        // Test numerical derivative
        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.repulsive(ptype, r + dr);
        auto [V_minus, _2] = pot.repulsive(ptype, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV_R, WithinRel(numerical_deriv, 1e-5));
    }

    SECTION("Attractive function") {
        Scalar r = 2.35;
        auto [V_A, dV_A] = pot.attractive(ptype, r);

        // V_A should be negative (attractive)
        REQUIRE(V_A < 0.0);
        // dV_A should be positive (magnitude decreases with distance)
        REQUIRE(dV_A > 0.0);

        // Test numerical derivative
        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.attractive(ptype, r + dr);
        auto [V_minus, _2] = pot.attractive(ptype, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV_A, WithinRel(numerical_deriv, 1e-5));
    }
}

TEST_CASE("Kumagai angular function", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

    int eli = 0;  // Si
    int ptype = 0;

    SECTION("Angular function values") {
        // Test at various angles
        for (Scalar cos_theta : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
            auto [g, dg] = pot.angular_function(eli, 0, 0, ptype, ptype, cos_theta);

            // g should be positive for these parameter sets
            REQUIRE(std::isfinite(g));
            REQUIRE(std::isfinite(dg));
        }
    }

    SECTION("Angular function derivative") {
        Scalar cos_theta = 0.3;
        auto [g, dg] = pot.angular_function(eli, 0, 0, ptype, ptype, cos_theta);

        // Test numerical derivative
        const Scalar dcos = 1e-6;
        auto [g_plus, _1] = pot.angular_function(eli, 0, 0, ptype, ptype, cos_theta + dcos);
        auto [g_minus, _2] = pot.angular_function(eli, 0, 0, ptype, ptype, cos_theta - dcos);
        Scalar numerical_deriv = (g_plus - g_minus) / (2 * dcos);

        REQUIRE_THAT(dg, WithinRel(numerical_deriv, 1e-4));
    }
}

TEST_CASE("Kumagai distance function", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

    int ptype = 0;

    SECTION("Distance function values") {
        Scalar r_ij = 2.35;
        Scalar r_ik = 2.40;
        auto [h, dh_drik, dh_drij] = pot.distance_function(0, 0, 0, ptype, ptype, r_ij, r_ik);

        // h = exp(alpha * dr), with dr = r_ij - r_ik = -0.05
        // For alpha > 0 and dr < 0, h < 1
        REQUIRE(h > 0.0);
        REQUIRE(std::isfinite(dh_drik));
        REQUIRE(std::isfinite(dh_drij));

        // dh/dr_ik should be opposite sign to dh/dr_ij
        REQUIRE_THAT(dh_drik, WithinAbs(-dh_drij, 1e-10));
    }

    SECTION("Distance function numerical derivatives") {
        Scalar r_ij = 2.35;
        Scalar r_ik = 2.40;
        auto [h, dh_drik, dh_drij] = pot.distance_function(0, 0, 0, ptype, ptype, r_ij, r_ik);

        const Scalar dr = 1e-6;

        // Test dh/dr_ik
        auto [h_plus_ik, _1, _2] = pot.distance_function(0, 0, 0, ptype, ptype, r_ij, r_ik + dr);
        auto [h_minus_ik, _3, _4] = pot.distance_function(0, 0, 0, ptype, ptype, r_ij, r_ik - dr);
        Scalar numerical_dh_drik = (h_plus_ik - h_minus_ik) / (2 * dr);
        REQUIRE_THAT(dh_drik, WithinRel(numerical_dh_drik, 1e-4));

        // Test dh/dr_ij
        auto [h_plus_ij, _5, _6] = pot.distance_function(0, 0, 0, ptype, ptype, r_ij + dr, r_ik);
        auto [h_minus_ij, _7, _8] = pot.distance_function(0, 0, 0, ptype, ptype, r_ij - dr, r_ik);
        Scalar numerical_dh_drij = (h_plus_ij - h_minus_ij) / (2 * dr);
        REQUIRE_THAT(dh_drij, WithinRel(numerical_dh_drij, 1e-4));
    }
}

TEST_CASE("Kumagai bond order function", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

    int eli = 0;
    int ptype = 0;

    SECTION("Bond order at z=0") {
        auto [b, db] = pot.bond_order(eli, ptype, 0.0);
        REQUIRE_THAT(b, WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(db, WithinAbs(0.0, 1e-10));
    }

    SECTION("Bond order decreases with z") {
        auto [b1, _1] = pot.bond_order(eli, ptype, 1.0);
        auto [b2, _2] = pot.bond_order(eli, ptype, 2.0);
        auto [b3, _3] = pot.bond_order(eli, ptype, 3.0);

        REQUIRE(b1 < 1.0);
        REQUIRE(b2 < b1);
        REQUIRE(b3 < b2);
    }

    SECTION("Bond order derivative") {
        Scalar z = 2.0;
        auto [b, db] = pot.bond_order(eli, ptype, z);

        // Numerical derivative
        const Scalar dz = 1e-6;
        auto [b_plus, _1] = pot.bond_order(eli, ptype, z + dz);
        auto [b_minus, _2] = pot.bond_order(eli, ptype, z - dz);
        Scalar numerical_db = (b_plus - b_minus) / (2 * dz);

        REQUIRE_THAT(db, WithinRel(numerical_db, 1e-4));
    }
}

TEST_CASE("Kumagai Si dimer", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 14;  // Si
    system.atomic_numbers()(1) = 14;  // Si

    SECTION("Dimer energy varies with distance") {
        NeighborList nl;
        nl.set_cutoff(pot.cutoff());

        std::vector<Scalar> distances = {2.0, 2.2, 2.35, 2.5, 2.7};
        std::vector<Scalar> energies;

        for (Scalar r : distances) {
            system.positions().col(0) << 10.0, 10.0, 10.0;
            system.positions().col(1) << 10.0 + r, 10.0, 10.0;
            system.positions_changed();

            nl.update(system);
            auto result = pot.compute(system, nl, false, false);
            energies.push_back(result.energy);
        }

        // Energy should have a minimum (not at the smallest distance)
        // Find minimum index
        auto min_it = std::min_element(energies.begin(), energies.end());
        size_t min_idx = std::distance(energies.begin(), min_it);

        // Minimum shouldn't be at the endpoints
        REQUIRE(min_idx > 0);
        REQUIRE(min_idx < energies.size() - 1);
    }
}

TEST_CASE("Kumagai Si trimer", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

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

    SECTION("Trimer energy depends on angle") {
        Scalar r = 2.35;  // Si-Si distance

        NeighborList nl;
        nl.set_cutoff(pot.cutoff());

        std::vector<Scalar> angles = {60.0, 90.0, 109.47, 120.0, 180.0};
        std::vector<Scalar> energies;

        for (Scalar angle_deg : angles) {
            Scalar angle_rad = angle_deg * M_PI / 180.0;

            system.positions().col(0) << 15.0, 15.0, 15.0;
            system.positions().col(1) << 15.0 + r, 15.0, 15.0;
            system.positions().col(2) << 15.0 + r * std::cos(angle_rad),
                                         15.0 + r * std::sin(angle_rad), 15.0;
            system.positions_changed();

            nl.update(system);
            auto result = pot.compute(system, nl, false, false);
            energies.push_back(result.energy);
        }

        // Energies should vary with angle
        Scalar min_e = *std::min_element(energies.begin(), energies.end());
        Scalar max_e = *std::max_element(energies.begin(), energies.end());
        REQUIRE(max_e > min_e);
    }
}

TEST_CASE("Kumagai numerical force test", "[Kumagai]") {
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

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

    // Asymmetric configuration for thorough testing
    system.positions().col(0) << 15.0, 15.0, 15.0;
    system.positions().col(1) << 17.3, 15.2, 15.1;
    system.positions().col(2) << 15.8, 17.1, 15.3;

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

TEST_CASE("Kumagai screened numerical force test", "[KumagaiScr]") {
    Kumagai<true> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

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

    // Linear arrangement for screening test
    Scalar r = 2.5;
    system.positions().col(0) << 15.0, 15.0, 15.0;
    system.positions().col(1) << 15.0 + r, 15.0, 15.0;
    system.positions().col(2) << 15.0 + 2*r, 15.0, 15.0;

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

TEST_CASE("Kumagai larger system numerical forces", "[Kumagai]") {
    // Test with a 4-atom configuration
    Kumagai<false> pot;
    pot.load_parameters("Kumagai_CompMaterSci_39_457_Si");

    AtomicSystem system(4);

    Mat3 cell;
    cell << 30.0, 0.0, 0.0,
            0.0, 30.0, 0.0,
            0.0, 0.0, 30.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    // All Si
    for (int i = 0; i < 4; ++i) {
        system.atomic_numbers()(i) = 14;
    }

    // Tetrahedral-like configuration
    Scalar a = 2.35;  // Si-Si distance
    system.positions().col(0) << 15.0, 15.0, 15.0;
    system.positions().col(1) << 15.0 + a, 15.0, 15.0;
    system.positions().col(2) << 15.0 + a/2, 15.0 + a*0.866, 15.0;
    system.positions().col(3) << 15.0 + a/2, 15.0 + a*0.289, 15.0 + a*0.816;

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    // Analytical forces
    system.zero_forces();
    auto result = pot.compute(system, nl, true, false);
    Array3X analytical_forces = system.forces();

    // Numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 4);

    for (int atom = 0; atom < 4; ++atom) {
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
    for (int atom = 0; atom < 4; ++atom) {
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
