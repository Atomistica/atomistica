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

TEST_CASE("REBO2 initialization", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    REQUIRE(pot.cutoff() > 0.0);
    REQUIRE(pot.element_type(6) == REBO2_C);   // Carbon
    REQUIRE(pot.element_type(1) == REBO2_H);   // Hydrogen
    REQUIRE(pot.element_type(14) == -1);       // Silicon not supported
}

TEST_CASE("REBO2 pair type", "[REBO2]") {
    REQUIRE(REBO2::pair_type(REBO2_C, REBO2_C) == REBO2_C_C);
    REQUIRE(REBO2::pair_type(REBO2_C, REBO2_H) == REBO2_C_H);
    REQUIRE(REBO2::pair_type(REBO2_H, REBO2_C) == REBO2_C_H);
    REQUIRE(REBO2::pair_type(REBO2_H, REBO2_H) == REBO2_H_H);
}

TEST_CASE("REBO2 repulsive function", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    SECTION("C-C repulsive") {
        Scalar r = 1.5;
        auto [V, dV] = pot.repulsive(REBO2_C_C, r);

        // V_R should be positive
        REQUIRE(V > 0.0);
        // dV_R should be negative (energy decreases with distance)
        REQUIRE(dV < 0.0);

        // Test numerical derivative
        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.repulsive(REBO2_C_C, r + dr);
        auto [V_minus, _2] = pot.repulsive(REBO2_C_C, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV, WithinRel(numerical_deriv, 1e-5));
    }

    SECTION("C-H repulsive") {
        Scalar r = 1.2;
        auto [V, dV] = pot.repulsive(REBO2_C_H, r);

        REQUIRE(V > 0.0);
        REQUIRE(dV < 0.0);

        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.repulsive(REBO2_C_H, r + dr);
        auto [V_minus, _2] = pot.repulsive(REBO2_C_H, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV, WithinRel(numerical_deriv, 1e-5));
    }

    SECTION("H-H repulsive") {
        Scalar r = 1.0;
        auto [V, dV] = pot.repulsive(REBO2_H_H, r);

        REQUIRE(V > 0.0);
        REQUIRE(dV < 0.0);

        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.repulsive(REBO2_H_H, r + dr);
        auto [V_minus, _2] = pot.repulsive(REBO2_H_H, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV, WithinRel(numerical_deriv, 1e-5));
    }
}

TEST_CASE("REBO2 attractive function", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    SECTION("C-C attractive") {
        Scalar r = 1.5;
        auto [V, dV] = pot.attractive(REBO2_C_C, r);

        // V_A should be negative
        REQUIRE(V < 0.0);
        // dV_A should be positive (magnitude decreases with distance)
        REQUIRE(dV > 0.0);

        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.attractive(REBO2_C_C, r + dr);
        auto [V_minus, _2] = pot.attractive(REBO2_C_C, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV, WithinRel(numerical_deriv, 1e-5));
    }

    SECTION("C-H attractive") {
        Scalar r = 1.2;
        auto [V, dV] = pot.attractive(REBO2_C_H, r);

        REQUIRE(V < 0.0);
        REQUIRE(dV > 0.0);

        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.attractive(REBO2_C_H, r + dr);
        auto [V_minus, _2] = pot.attractive(REBO2_C_H, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV, WithinRel(numerical_deriv, 1e-5));
    }

    SECTION("H-H attractive") {
        Scalar r = 1.0;
        auto [V, dV] = pot.attractive(REBO2_H_H, r);

        REQUIRE(V < 0.0);
        REQUIRE(dV > 0.0);

        const Scalar dr = 1e-6;
        auto [V_plus, _1] = pot.attractive(REBO2_H_H, r + dr);
        auto [V_minus, _2] = pot.attractive(REBO2_H_H, r - dr);
        Scalar numerical_deriv = (V_plus - V_minus) / (2 * dr);
        REQUIRE_THAT(dV, WithinRel(numerical_deriv, 1e-5));
    }
}

TEST_CASE("REBO2 angular function for H", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    SECTION("H angular function values") {
        // Test at various angles
        for (Scalar cos_theta : {-1.0, -0.5, 0.0, 0.5, 1.0}) {
            auto [g, dg, dg_dN] = pot.angular_function(REBO2_H, cos_theta, 0.0);
            REQUIRE(std::isfinite(g));
            REQUIRE(std::isfinite(dg));
            REQUIRE(dg_dN == 0.0);  // H has no N dependence
        }
    }

    SECTION("H angular function derivative") {
        Scalar cos_theta = 0.3;
        auto [g, dg, _] = pot.angular_function(REBO2_H, cos_theta, 0.0);

        const Scalar dcos = 1e-6;
        auto [g_plus, _1, _2] = pot.angular_function(REBO2_H, cos_theta + dcos, 0.0);
        auto [g_minus, _3, _4] = pot.angular_function(REBO2_H, cos_theta - dcos, 0.0);
        Scalar numerical_deriv = (g_plus - g_minus) / (2 * dcos);

        REQUIRE_THAT(dg, WithinRel(numerical_deriv, 1e-4));
    }
}

TEST_CASE("REBO2 distance weighting", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    SECTION("C-C to C-C (no weighting)") {
        auto [h, dh] = pot.distance_weight(REBO2_C_C, REBO2_C_C, 0.1);
        // C-C to C-C should have weight 1 (no hydrogen involved)
        REQUIRE_THAT(h, WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(dh, WithinAbs(0.0, 1e-10));
    }

    SECTION("C-H to C-H") {
        Scalar dr = 0.1;
        auto [h, dh] = pot.distance_weight(REBO2_C_H, REBO2_C_H, dr);

        // Should apply exponential weighting
        REQUIRE(h > 0.0);
        REQUIRE(std::isfinite(dh));

        // Test numerical derivative
        const Scalar ddr = 1e-6;
        auto [h_plus, _1] = pot.distance_weight(REBO2_C_H, REBO2_C_H, dr + ddr);
        auto [h_minus, _2] = pot.distance_weight(REBO2_C_H, REBO2_C_H, dr - ddr);
        Scalar numerical_deriv = (h_plus - h_minus) / (2 * ddr);
        REQUIRE_THAT(dh, WithinRel(numerical_deriv, 1e-4));
    }
}

TEST_CASE("REBO2 C-C dimer", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 6;  // C
    system.atomic_numbers()(1) = 6;  // C

    SECTION("Dimer energy varies with distance") {
        NeighborList nl;
        nl.set_cutoff(pot.cutoff());

        // For a C-C dimer with b=1 (no angular contributions),
        // the equilibrium distance is around 1.2-1.3 A (similar to double bond)
        std::vector<Scalar> distances = {1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
        std::vector<Scalar> energies;

        for (Scalar r : distances) {
            system.positions().col(0) << 10.0, 10.0, 10.0;
            system.positions().col(1) << 10.0 + r, 10.0, 10.0;
            system.positions_changed();

            nl.update(system);
            auto result = pot.compute(system, nl, false, false);
            energies.push_back(result.energy);
        }

        // Energy should have a minimum somewhere in this range
        auto min_it = std::min_element(energies.begin(), energies.end());
        size_t min_idx = std::distance(energies.begin(), min_it);

        // Minimum shouldn't be at the endpoints
        REQUIRE(min_idx > 0);
        REQUIRE(min_idx < energies.size() - 1);
    }
}

TEST_CASE("REBO2 C-H dimer", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 6;  // C
    system.atomic_numbers()(1) = 1;  // H

    Scalar r = 1.09;  // Equilibrium C-H distance
    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r, 10.0, 10.0;

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    auto result = pot.compute(system, nl, false, false);

    // Energy should be negative (bound state)
    REQUIRE(result.energy < 0.0);
}

TEST_CASE("REBO2 H-H dimer", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 1;  // H
    system.atomic_numbers()(1) = 1;  // H

    Scalar r = 0.74;  // Near H-H equilibrium distance
    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 10.0 + r, 10.0, 10.0;

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    auto result = pot.compute(system, nl, false, false);

    // Energy should be negative (bound state)
    REQUIRE(result.energy < 0.0);
}

TEST_CASE("REBO2 numerical force test (dimer)", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    AtomicSystem system(2);

    Mat3 cell;
    cell << 20.0, 0.0, 0.0,
            0.0, 20.0, 0.0,
            0.0, 0.0, 20.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    system.atomic_numbers()(0) = 6;  // C
    system.atomic_numbers()(1) = 6;  // C

    system.positions().col(0) << 10.0, 10.0, 10.0;
    system.positions().col(1) << 11.5, 10.0, 10.0;

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    // Analytical forces
    system.zero_forces();
    auto result = pot.compute(system, nl, true, false);
    Array3X analytical_forces = system.forces();

    // Numerical forces
    const Scalar dx = 1e-6;
    Array3X numerical_forces = Array3X::Zero(3, 2);

    for (int atom = 0; atom < 2; ++atom) {
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
    for (int atom = 0; atom < 2; ++atom) {
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

TEST_CASE("REBO2 methane CH4", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    // Create methane molecule (tetrahedral geometry)
    AtomicSystem system(5);

    Mat3 cell;
    cell << 30.0, 0.0, 0.0,
            0.0, 30.0, 0.0,
            0.0, 0.0, 30.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    // Carbon at center
    system.atomic_numbers()(0) = 6;
    system.positions().col(0) << 15.0, 15.0, 15.0;

    // Four hydrogens in tetrahedral arrangement
    Scalar r_CH = 1.09;  // C-H bond length
    Scalar angle = std::acos(-1.0/3.0);  // Tetrahedral angle

    system.atomic_numbers()(1) = 1;
    system.positions().col(1) << 15.0 + r_CH, 15.0, 15.0;

    system.atomic_numbers()(2) = 1;
    system.positions().col(2) << 15.0 - r_CH/3, 15.0 + r_CH * std::sqrt(8.0/9.0), 15.0;

    system.atomic_numbers()(3) = 1;
    system.positions().col(3) << 15.0 - r_CH/3, 15.0 - r_CH * std::sqrt(2.0/9.0),
                                 15.0 + r_CH * std::sqrt(2.0/3.0);

    system.atomic_numbers()(4) = 1;
    system.positions().col(4) << 15.0 - r_CH/3, 15.0 - r_CH * std::sqrt(2.0/9.0),
                                 15.0 - r_CH * std::sqrt(2.0/3.0);

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    auto result = pot.compute(system, nl, true, false);

    // Methane should have negative total energy
    REQUIRE(result.energy < 0.0);

    // Forces should sum to zero (no net force on molecule)
    Vec3 total_force = Vec3::Zero();
    for (int i = 0; i < 5; ++i) {
        total_force += system.forces().col(i).matrix();
    }
    REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-8));
}

TEST_CASE("REBO2 Table2D interpolation", "[REBO2]") {
    // Test 2D bicubic interpolation
    Table2D table;

    // Create a simple quadratic function: f(x,y) = x^2 + y^2
    std::vector<std::vector<Scalar>> values(4, std::vector<Scalar>(4));
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            values[i][j] = i*i + j*j;
        }
    }

    table.init(3, 3, values);

    SECTION("Table is valid after init") {
        REQUIRE(table.is_valid());
    }

    SECTION("Values at grid points") {
        auto [v00, dx00, dy00] = table.eval(0.0, 0.0);
        REQUIRE_THAT(v00, WithinAbs(0.0, 1e-10));

        auto [v11, dx11, dy11] = table.eval(1.0, 1.0);
        REQUIRE_THAT(v11, WithinAbs(2.0, 1e-10));

        auto [v22, dx22, dy22] = table.eval(2.0, 2.0);
        REQUIRE_THAT(v22, WithinAbs(8.0, 1e-10));
    }

    SECTION("Interpolated values") {
        // At (1.5, 1.5), exact value would be 1.5^2 + 1.5^2 = 4.5
        // Interpolation may not be exact for quadratic function
        auto [v, dx, dy] = table.eval(1.5, 1.5);
        REQUIRE_THAT(v, WithinAbs(4.5, 1.0));  // Allow some interpolation error
    }
}

TEST_CASE("REBO2 Table3D interpolation", "[REBO2]") {
    // Test 3D tricubic interpolation
    Table3D table;

    // Create a simple function: f(x,y,z) = x + y + z
    std::vector<std::vector<std::vector<Scalar>>> values(3,
        std::vector<std::vector<Scalar>>(3, std::vector<Scalar>(3)));
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                values[i][j][k] = i + j + k;
            }
        }
    }

    table.init(2, 2, 2, values);

    SECTION("Table is valid after init") {
        REQUIRE(table.is_valid());
    }

    SECTION("Values at grid points") {
        auto [v000, dx, dy, dz] = table.eval(0.0, 0.0, 0.0);
        REQUIRE_THAT(v000, WithinAbs(0.0, 1e-10));

        auto [v111, dx1, dy1, dz1] = table.eval(1.0, 1.0, 1.0);
        REQUIRE_THAT(v111, WithinAbs(3.0, 1e-10));
    }
}

TEST_CASE("REBO2 bond order function", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    SECTION("b(0) = 1 for zero angular contribution") {
        auto [b, db] = pot.bond_order_func(REBO2_C, 0.0);
        REQUIRE_THAT(b, WithinAbs(1.0, 1e-10));
    }

    SECTION("b decreases with increasing z") {
        auto [b1, db1] = pot.bond_order_func(REBO2_C, 1.0);
        auto [b2, db2] = pot.bond_order_func(REBO2_C, 2.0);
        auto [b3, db3] = pot.bond_order_func(REBO2_C, 3.0);

        REQUIRE(b1 < 1.0);
        REQUIRE(b2 < b1);
        REQUIRE(b3 < b2);
    }

    SECTION("Derivative is negative") {
        auto [b, db] = pot.bond_order_func(REBO2_C, 1.0);
        REQUIRE(db < 0.0);
    }
}

TEST_CASE("REBO2 distance weight function", "[REBO2]") {
    REBO2 pot;
    pot.load_default_parameters();

    SECTION("C-C bonds have weight 1 for equal distances") {
        auto [h, dh] = pot.distance_weight(REBO2_C_C, REBO2_C_C, 0.0);
        REQUIRE_THAT(h, WithinAbs(1.0, 1e-10));
        REQUIRE_THAT(dh, WithinAbs(0.0, 1e-10));
    }

    SECTION("Weight applies when bond type sum > 4") {
        // C-C + C-H = 1 + 3 = 4, so no weight applied (returns 1.0)
        auto [h_cc_ch, dh] = pot.distance_weight(REBO2_C_C, REBO2_C_H, 0.0);
        REQUIRE_THAT(h_cc_ch, WithinAbs(1.0, 1e-10));

        // C-H + C-H = 3 + 3 = 6 > 4, so weight applied
        auto [h1, dh1] = pot.distance_weight(REBO2_C_H, REBO2_C_H, 0.0);
        REQUIRE_THAT(h1, WithinAbs(1.0, 1e-10));  // At dr=0, conear is 1.0 for same type

        auto [h2, dh2] = pot.distance_weight(REBO2_C_H, REBO2_C_H, 0.3);
        // Weight should increase with positive dr due to exp(lambda * dr)
        REQUIRE(h2 > h1);
    }
}

TEST_CASE("REBO2 ethane molecule", "[REBO2]") {
    // Test C2H6 (ethane) which has angular contributions
    REBO2 pot;
    pot.load_default_parameters();

    // 2 carbons + 6 hydrogens
    AtomicSystem system(8);

    Mat3 cell;
    cell << 30.0, 0.0, 0.0,
            0.0, 30.0, 0.0,
            0.0, 0.0, 30.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    Scalar r_CC = 1.54;  // C-C bond length
    Scalar r_CH = 1.09;  // C-H bond length

    // Two carbons along x-axis
    system.atomic_numbers()(0) = 6;
    system.positions().col(0) << 15.0, 15.0, 15.0;

    system.atomic_numbers()(1) = 6;
    system.positions().col(1) << 15.0 + r_CC, 15.0, 15.0;

    // Three hydrogens on first carbon (staggered)
    for (int i = 0; i < 3; ++i) {
        system.atomic_numbers()(2 + i) = 1;
        Scalar angle = 2.0 * M_PI * i / 3.0;
        system.positions().col(2 + i) << 15.0 - r_CH * std::cos(109.5*M_PI/180.0),
                                          15.0 + r_CH * std::sin(109.5*M_PI/180.0) * std::cos(angle),
                                          15.0 + r_CH * std::sin(109.5*M_PI/180.0) * std::sin(angle);
    }

    // Three hydrogens on second carbon
    for (int i = 0; i < 3; ++i) {
        system.atomic_numbers()(5 + i) = 1;
        Scalar angle = 2.0 * M_PI * i / 3.0 + M_PI/3.0;  // Staggered
        system.positions().col(5 + i) << 15.0 + r_CC + r_CH * std::cos(109.5*M_PI/180.0),
                                          15.0 + r_CH * std::sin(109.5*M_PI/180.0) * std::cos(angle),
                                          15.0 + r_CH * std::sin(109.5*M_PI/180.0) * std::sin(angle);
    }

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    auto result = pot.compute(system, nl, true, false);

    // Ethane should have negative energy
    REQUIRE(result.energy < 0.0);

    // Forces should sum to zero
    Vec3 total_force = Vec3::Zero();
    for (int i = 0; i < 8; ++i) {
        total_force += system.forces().col(i).matrix();
    }
    REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-7));
}

TEST_CASE("REBO2 graphene-like C trimer", "[REBO2]") {
    // Test 3-atom linear carbon chain to verify angular contributions
    REBO2 pot;
    pot.load_default_parameters();

    AtomicSystem system(3);

    Mat3 cell;
    cell << 30.0, 0.0, 0.0,
            0.0, 30.0, 0.0,
            0.0, 0.0, 30.0;
    system.set_cell(cell);
    system.pbc() = {false, false, false};

    Scalar r = 1.42;  // C-C bond length in graphene

    system.atomic_numbers()(0) = 6;
    system.positions().col(0) << 15.0, 15.0, 15.0;

    system.atomic_numbers()(1) = 6;
    system.positions().col(1) << 15.0 + r, 15.0, 15.0;

    system.atomic_numbers()(2) = 6;
    system.positions().col(2) << 15.0 + 2*r, 15.0, 15.0;

    NeighborList nl;
    nl.set_cutoff(pot.cutoff());
    nl.update(system);

    auto result = pot.compute(system, nl, true, false);

    // Chain should have negative energy
    REQUIRE(result.energy < 0.0);

    // Middle atom should feel different force than end atoms due to angular terms
    Scalar f_end = system.forces().col(0).matrix().norm();
    Scalar f_middle = system.forces().col(1).matrix().norm();

    // For symmetric chain, end atoms have outward forces, middle is balanced
    // The middle atom force should be smaller than end atoms for symmetric linear chain
    REQUIRE_THAT(f_middle, WithinAbs(0.0, 0.1));

    // Forces sum to zero
    Vec3 total_force = Vec3::Zero();
    for (int i = 0; i < 3; ++i) {
        total_force += system.forces().col(i).matrix();
    }
    REQUIRE_THAT(total_force.norm(), WithinAbs(0.0, 1e-8));
}
