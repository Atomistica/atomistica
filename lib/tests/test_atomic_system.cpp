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

#include <atomistica/core/atomic_system.hpp>

using namespace atomistica;
using Catch::Matchers::WithinRel;

TEST_CASE("AtomicSystem basic operations", "[AtomicSystem]") {
    AtomicSystem system(10);

    SECTION("Initial state") {
        REQUIRE(system.num_atoms() == 10);
        REQUIRE(system.cell().isApprox(Mat3::Identity()));
        REQUIRE(system.pbc() == std::array<bool, 3>{true, true, true});
    }

    SECTION("Resize") {
        system.resize(20);
        REQUIRE(system.num_atoms() == 20);
    }

    SECTION("Set cell") {
        Mat3 cell;
        cell << 5.0, 0.0, 0.0,
                0.0, 5.0, 0.0,
                0.0, 0.0, 5.0;
        system.set_cell(cell);
        REQUIRE(system.cell().isApprox(cell));
        REQUIRE_THAT(system.volume(), WithinRel(125.0, 1e-10));
    }

    SECTION("Non-orthorhombic cell") {
        Mat3 cell;
        cell << 5.0, 2.5, 0.0,
                0.0, 4.33, 0.0,
                0.0, 0.0, 5.0;
        system.set_cell(cell);
        REQUIRE_THAT(system.volume(), WithinRel(5.0 * 4.33 * 5.0, 1e-10));
    }

    SECTION("Set positions") {
        system.set_position(0, Vec3(1.0, 2.0, 3.0));
        REQUIRE_THAT(system.position(0)(0), WithinRel(1.0, 1e-10));
        REQUIRE_THAT(system.position(0)(1), WithinRel(2.0, 1e-10));
        REQUIRE_THAT(system.position(0)(2), WithinRel(3.0, 1e-10));
    }

    SECTION("Set atomic numbers") {
        system.atomic_numbers()(0) = 6;   // Carbon
        system.atomic_numbers()(1) = 14;  // Silicon
        REQUIRE(system.atomic_numbers()(0) == 6);
        REQUIRE(system.atomic_numbers()(1) == 14);
    }
}

TEST_CASE("AtomicSystem minimum image", "[AtomicSystem]") {
    AtomicSystem system(2);

    Mat3 cell;
    cell << 10.0, 0.0, 0.0,
            0.0, 10.0, 0.0,
            0.0, 0.0, 10.0;
    system.set_cell(cell);

    SECTION("No wrapping needed") {
        Vec3 dr(1.0, 2.0, 3.0);
        Vec3 result = system.minimum_image(dr);
        REQUIRE(result.isApprox(dr));
    }

    SECTION("Positive wrap") {
        Vec3 dr(8.0, 2.0, 3.0);
        Vec3 result = system.minimum_image(dr);
        REQUIRE_THAT(result(0), WithinRel(-2.0, 1e-10));
        REQUIRE_THAT(result(1), WithinRel(2.0, 1e-10));
        REQUIRE_THAT(result(2), WithinRel(3.0, 1e-10));
    }

    SECTION("Negative wrap") {
        Vec3 dr(-8.0, 2.0, 3.0);
        Vec3 result = system.minimum_image(dr);
        REQUIRE_THAT(result(0), WithinRel(2.0, 1e-10));
    }

    SECTION("Non-periodic direction") {
        system.pbc() = {true, false, true};
        Vec3 dr(8.0, 8.0, 8.0);
        Vec3 result = system.minimum_image(dr);
        REQUIRE_THAT(result(0), WithinRel(-2.0, 1e-10));
        REQUIRE_THAT(result(1), WithinRel(8.0, 1e-10));  // Not wrapped
        REQUIRE_THAT(result(2), WithinRel(-2.0, 1e-10));
    }
}

TEST_CASE("AtomicSystem wrap position", "[AtomicSystem]") {
    AtomicSystem system(1);

    Mat3 cell;
    cell << 10.0, 0.0, 0.0,
            0.0, 10.0, 0.0,
            0.0, 0.0, 10.0;
    system.set_cell(cell);

    SECTION("Position inside cell") {
        Vec3 r(5.0, 5.0, 5.0);
        Vec3 result = system.wrap_position(r);
        REQUIRE(result.isApprox(r));
    }

    SECTION("Position outside cell - positive") {
        Vec3 r(15.0, 5.0, 5.0);
        Vec3 result = system.wrap_position(r);
        REQUIRE_THAT(result(0), WithinRel(5.0, 1e-10));
    }

    SECTION("Position outside cell - negative") {
        Vec3 r(-3.0, 5.0, 5.0);
        Vec3 result = system.wrap_position(r);
        REQUIRE_THAT(result(0), WithinRel(7.0, 1e-10));
    }
}

TEST_CASE("PropertyMap", "[PropertyMap]") {
    PropertyMap props;

    SECTION("Add and get scalar property") {
        props.add<ArrayX>("charge", 10);
        REQUIRE(props.has("charge"));

        auto& charge = props.get<ArrayX>("charge");
        REQUIRE(charge.size() == 10);

        charge(0) = 1.5;
        REQUIRE_THAT(props.get<ArrayX>("charge")(0), WithinRel(1.5, 1e-10));
    }

    SECTION("Add and get vector property") {
        props.add<Array3X>("velocity", 10);
        REQUIRE(props.has("velocity"));

        auto& vel = props.get<Array3X>("velocity");
        REQUIRE(vel.cols() == 10);
        REQUIRE(vel.rows() == 3);
    }

    SECTION("Add and get integer property") {
        props.add<ArrayXi>("type", 10);
        REQUIRE(props.has("type"));

        auto& type = props.get<ArrayXi>("type");
        type(0) = 1;
        type(1) = 2;
        REQUIRE(props.get<ArrayXi>("type")(0) == 1);
        REQUIRE(props.get<ArrayXi>("type")(1) == 2);
    }

    SECTION("Remove property") {
        props.add<ArrayX>("temp", 10);
        REQUIRE(props.has("temp"));
        props.remove("temp");
        REQUIRE_FALSE(props.has("temp"));
    }
}
