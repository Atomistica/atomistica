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
#include <atomistica/core/neighbor_list.hpp>

#include <cmath>
#include <set>

using namespace atomistica;
using Catch::Matchers::WithinRel;

TEST_CASE("NeighborList basic operations", "[NeighborList]") {
    NeighborList nl;

    SECTION("Set cutoff") {
        nl.set_cutoff(5.0);
        REQUIRE_THAT(nl.cutoff(), WithinRel(5.0, 1e-10));
    }

    SECTION("Set Verlet shell") {
        nl.set_verlet_shell(0.5);
        REQUIRE_THAT(nl.verlet_shell(), WithinRel(0.5, 1e-10));
    }
}

TEST_CASE("NeighborList two atoms", "[NeighborList]") {
    AtomicSystem system(2);

    Mat3 cell;
    cell << 10.0, 0.0, 0.0,
            0.0, 10.0, 0.0,
            0.0, 0.0, 10.0;
    system.set_cell(cell);

    system.positions().col(0) << 0.0, 0.0, 0.0;
    system.positions().col(1) << 3.0, 0.0, 0.0;

    NeighborList nl;

    SECTION("Atoms within cutoff") {
        nl.set_cutoff(5.0);
        nl.update(system);

        REQUIRE(nl.num_atoms() == 2);
        REQUIRE(nl.num_neighbors(0) == 1);
        auto [begin, end] = nl.neighbors(0);
        REQUIRE(begin->index == 1);
    }

    SECTION("Atoms outside cutoff") {
        nl.set_cutoff(2.0);
        nl.update(system);

        REQUIRE(nl.num_neighbors(0) == 0);
        REQUIRE(nl.num_neighbors(1) == 0);
    }
}

TEST_CASE("NeighborList periodic boundary", "[NeighborList]") {
    AtomicSystem system(2);

    Mat3 cell;
    cell << 10.0, 0.0, 0.0,
            0.0, 10.0, 0.0,
            0.0, 0.0, 10.0;
    system.set_cell(cell);

    // Atoms near opposite boundaries
    system.positions().col(0) << 0.5, 5.0, 5.0;
    system.positions().col(1) << 9.5, 5.0, 5.0;  // Distance is 1.0 through PBC

    NeighborList nl;
    nl.set_cutoff(2.0);

    SECTION("With PBC - should find neighbor") {
        system.pbc() = {true, true, true};
        nl.update(system);

        REQUIRE(nl.num_neighbors(0) == 1);
        auto [begin, end] = nl.neighbors(0);
        const auto& neigh = *begin;
        REQUIRE(neigh.index == 1);
        // Cell shift should be non-zero
        REQUIRE((neigh.cell_shift[0] != 0 || neigh.cell_shift[1] != 0 || neigh.cell_shift[2] != 0));
    }

    SECTION("Without PBC in x - should not find neighbor") {
        system.pbc() = {false, true, true};
        nl.update(system);

        REQUIRE(nl.num_neighbors(0) == 0);
    }
}

TEST_CASE("NeighborList FCC lattice", "[NeighborList]") {
    // Create 2x2x2 FCC unit cells
    const Scalar a = 4.05;  // Aluminum lattice constant
    AtomicSystem system(32);  // 4 atoms per unit cell * 8 unit cells

    Mat3 cell;
    cell << 2*a, 0.0, 0.0,
            0.0, 2*a, 0.0,
            0.0, 0.0, 2*a;
    system.set_cell(cell);

    // FCC basis
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
                    system.atomic_numbers()(idx) = 13;  // Aluminum
                    ++idx;
                }
            }
        }
    }

    NeighborList nl;

    SECTION("First neighbor shell") {
        // First neighbor distance in FCC: a/sqrt(2)
        Scalar r1 = a / std::sqrt(2.0);
        nl.set_cutoff(r1 + 0.1);  // Slightly larger than first neighbor
        nl.update(system);

        // Each atom should have 12 first neighbors in FCC
        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            REQUIRE(nl.num_neighbors(i) == 12);
        }
    }

    SECTION("Second neighbor shell") {
        // Second neighbor distance in FCC: a
        nl.set_cutoff(a + 0.1);
        nl.update(system);

        // Each atom should have 12 + 6 = 18 neighbors (first + second shell)
        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            REQUIRE(nl.num_neighbors(i) == 18);
        }
    }
}

TEST_CASE("NeighborList symmetry", "[NeighborList]") {
    // Test that if j is neighbor of i, then i is neighbor of j
    AtomicSystem system(10);

    Mat3 cell;
    cell << 10.0, 0.0, 0.0,
            0.0, 10.0, 0.0,
            0.0, 0.0, 10.0;
    system.set_cell(cell);

    // Random positions
    srand(42);
    for (std::size_t i = 0; i < system.num_atoms(); ++i) {
        system.positions().col(i) << 10.0 * rand() / RAND_MAX,
                                     10.0 * rand() / RAND_MAX,
                                     10.0 * rand() / RAND_MAX;
    }

    NeighborList nl;
    nl.set_cutoff(5.0);
    nl.update(system);

    // Build set of unique pairs from neighbor list
    std::set<std::pair<std::size_t, std::size_t>> pairs;
    for (std::size_t i = 0; i < system.num_atoms(); ++i) {
        auto [begin, end] = nl.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;
            pairs.insert({std::min(i, j), std::max(i, j)});
        }
    }

    // Full neighbor list: num_pairs() counts all entries (both i->j and j->i)
    // unique pairs = num_pairs / 2
    REQUIRE(pairs.size() * 2 == nl.num_pairs());
}

// ============================================================================
// Self-image bond tests (atoms seeing their own periodic images)
// ============================================================================

TEST_CASE("NeighborList includes self-image bonds through PBC", "[neighbor_list][self_image]") {
    SECTION("Single atom sees its own periodic images") {
        // 1-atom simple cubic cell: atom should see itself at a, a*sqrt(2), etc.
        AtomicSystem sys(1);
        Mat3 cell = Mat3::Identity() * 5.0;  // 5 Å cube
        sys.set_cell(cell);
        sys.pbc() = {true, true, true};
        sys.positions() << 0.0, 0.0, 0.0;
        sys.atomic_numbers()(0) = 14;

        NeighborList nl;
        nl.set_cutoff(7.5);  // large enough to see multiple self-images
        nl.update(sys);

        // Should find 6 self-images at distance 5 Å (nearest along ±x,±y,±z)
        int count_5A = 0;
        for (const auto& nb : [&]{ auto [b,e] = nl.neighbors(0); return std::vector<Neighbor>(b,e); }()) {
            REQUIRE(nb.index == 0);  // all should be self-images
            Vec3 dr = (cell.col(0)*nb.cell_shift[0] +
                       cell.col(1)*nb.cell_shift[1] +
                       cell.col(2)*nb.cell_shift[2]);
            Scalar r = dr.norm();
            if (std::abs(r - 5.0) < 0.01) count_5A++;
        }
        REQUIRE(count_5A == 6);
    }

    SECTION("2-atom BCC cell: atom sees full 8+6 coordination shell") {
        // BCC: 2 atoms per conventional cell
        // Atom 0 at corner, atom 1 at body center
        // With small cell, 2nd NN of atom 0 are self-images via PBC
        const Scalar a = 3.165;  // BCC-W lattice constant
        AtomicSystem sys(2);
        Mat3 cell = Mat3::Identity() * a;
        sys.set_cell(cell);
        sys.pbc() = {true, true, true};
        sys.positions().col(0) << 0.0, 0.0, 0.0;
        sys.positions().col(1) << a/2, a/2, a/2;
        sys.atomic_numbers()(0) = 74;
        sys.atomic_numbers()(1) = 74;

        NeighborList nl;
        nl.set_cutoff(4.0);  // includes both 1st NN (2.74 Å) and 2nd NN (3.165 Å)
        nl.update(sys);

        // Atom 0 should have:
        // - 8 bonds to atom 1 at r = a*sqrt(3)/2 ≈ 2.741 Å (1st NN)
        // - 6 bonds to atom 0 (self) at r = a ≈ 3.165 Å (2nd NN, self-images)
        int count_1nn = 0, count_2nn_self = 0;
        for (const auto& nb : [&]{ auto [b,e] = nl.neighbors(0); return std::vector<Neighbor>(b,e); }()) {
            Vec3 dr = sys.positions().col(nb.index).matrix();
            dr -= sys.positions().col(0).matrix();
            dr += cell.col(0)*nb.cell_shift[0] + cell.col(1)*nb.cell_shift[1] + cell.col(2)*nb.cell_shift[2];
            Scalar r = dr.norm();

            if (nb.index == 1 && std::abs(r - a*std::sqrt(3.0)/2) < 0.01)
                count_1nn++;
            else if (nb.index == 0 && std::abs(r - a) < 0.01)
                count_2nn_self++;
        }

        REQUIRE(count_1nn == 8);    // 8 nearest neighbors
        REQUIRE(count_2nn_self == 6); // 6 self-image 2nd NN
    }
}
