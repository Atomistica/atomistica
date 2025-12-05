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
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <algorithm>

#include <atomistica/tightbinding/tightbinding.hpp>
#include <atomistica/core/atomic_system.hpp>
#include <atomistica/core/neighbor_list.hpp>

using namespace atomistica;
using namespace atomistica::tb;
using Catch::Approx;

// Test SK integrals (arbitrary values)
std::array<Scalar, NUM_SK_INTEGRALS> make_test_sk() {
    std::array<Scalar, NUM_SK_INTEGRALS> sk;
    sk[0] = 0.5;   // dds
    sk[1] = 0.3;   // ddp
    sk[2] = 0.1;   // ddd
    sk[3] = 0.4;   // pds
    sk[4] = 0.2;   // pdp
    sk[5] = 0.6;   // pps
    sk[6] = -0.2;  // ppp
    sk[7] = 0.3;   // sds
    sk[8] = 0.5;   // sps
    sk[9] = -0.8;  // sss
    return sk;
}

TEST_CASE("Slater-Koster: s-s interaction", "[tightbinding]") {
    Vec3 c_z(0.0, 0.0, 1.0);
    auto sk = make_test_sk();

    Scalar result = transform_orb(1, 1, c_z, sk);
    REQUIRE(result == Approx(sk[9]).epsilon(1e-10));
}

TEST_CASE("Slater-Koster: s-p interaction along z", "[tightbinding]") {
    Vec3 c_z(0.0, 0.0, 1.0);
    auto sk = make_test_sk();

    // s-pz along z-axis: should give sps
    Scalar result = transform_orb(1, 4, c_z, sk);
    REQUIRE(result == Approx(sk[8]).epsilon(1e-10));

    // s-px along z-axis: should give 0
    result = transform_orb(1, 2, c_z, sk);
    REQUIRE(result == Approx(0.0).margin(1e-10));

    // s-py along z-axis: should give 0
    result = transform_orb(1, 3, c_z, sk);
    REQUIRE(result == Approx(0.0).margin(1e-10));
}

TEST_CASE("Slater-Koster: s-p interaction along x", "[tightbinding]") {
    Vec3 c_x(1.0, 0.0, 0.0);
    auto sk = make_test_sk();

    // s-px along x-axis: should give sps
    Scalar result = transform_orb(1, 2, c_x, sk);
    REQUIRE(result == Approx(sk[8]).epsilon(1e-10));
}

TEST_CASE("Slater-Koster: p-p sigma interaction", "[tightbinding]") {
    Vec3 c_z(0.0, 0.0, 1.0);
    auto sk = make_test_sk();

    // pz-pz along z-axis: should give pps (sigma only)
    Scalar result = transform_orb(4, 4, c_z, sk);
    REQUIRE(result == Approx(sk[5]).epsilon(1e-10));
}

TEST_CASE("Slater-Koster: p-p pi interaction", "[tightbinding]") {
    Vec3 c_z(0.0, 0.0, 1.0);
    auto sk = make_test_sk();

    // px-px along z-axis: should give ppp (pi only)
    Scalar result = transform_orb(2, 2, c_z, sk);
    REQUIRE(result == Approx(sk[6]).epsilon(1e-10));
}

TEST_CASE("Slater-Koster: symmetry", "[tightbinding]") {
    Vec3 c_111 = Vec3(1.0, 1.0, 1.0).normalized();
    auto sk = make_test_sk();

    // s-p should be antisymmetric (odd parity)
    Scalar sp = transform_orb(1, 2, c_111, sk);
    Scalar ps = transform_orb(2, 1, c_111, sk);
    REQUIRE(sp == Approx(-ps).epsilon(1e-10));

    // p-p should be symmetric (even parity)
    Scalar pp12 = transform_orb(2, 3, c_111, sk);
    Scalar pp21 = transform_orb(3, 2, c_111, sk);
    REQUIRE(pp12 == Approx(pp21).epsilon(1e-10));
}

TEST_CASE("TBElementParams: Carbon MIO", "[tightbinding]") {
    auto c = parameters::carbon_mio();

    REQUIRE(c.symbol == "C");
    REQUIRE(c.atomic_number == 6);
    REQUIRE(c.num_orbitals == 4);  // sp basis
    REQUIRE(c.l_max == 1);
    REQUIRE(c.valence_electrons == Approx(4.0).epsilon(1e-10));
    REQUIRE(c.hubbard_U > 0.0);
}

TEST_CASE("TBElementParams: Hydrogen MIO", "[tightbinding]") {
    auto h = parameters::hydrogen_mio();

    REQUIRE(h.symbol == "H");
    REQUIRE(h.atomic_number == 1);
    REQUIRE(h.num_orbitals == 1);  // s only
    REQUIRE(h.l_max == 0);
    REQUIRE(h.valence_electrons == Approx(1.0).epsilon(1e-10));
}

TEST_CASE("DenseHamiltonian: resize", "[tightbinding]") {
    DenseHamiltonian ham;
    ham.resize(10, 40);

    REQUIRE(ham.num_atoms == 10);
    REQUIRE(ham.num_orbitals == 40);
    REQUIRE(ham.H.rows() == 40);
    REQUIRE(ham.H.cols() == 40);
    REQUIRE(ham.S.rows() == 40);
    REQUIRE(ham.S.cols() == 40);
    REQUIRE(ham.eigenvalues.size() == 40);
    REQUIRE(ham.charges.size() == 10);
}

TEST_CASE("DenseHamiltonian: clear matrices", "[tightbinding]") {
    DenseHamiltonian ham;
    ham.resize(5, 20);

    // Set some non-zero values
    ham.H(0, 0) = 1.0;
    ham.S(1, 1) = 2.0;

    ham.clear_matrices();

    REQUIRE(ham.H(0, 0) == 0.0);
    REQUIRE(ham.S(1, 1) == 0.0);
}

TEST_CASE("Fermi-Dirac: zero temperature", "[tightbinding]") {
    // Below Fermi level
    REQUIRE(fermi_dirac(-1.0, 0.0, 1e-10) == Approx(1.0).epsilon(1e-10));

    // Above Fermi level
    REQUIRE(fermi_dirac(1.0, 0.0, 1e-10) == Approx(0.0).margin(1e-10));
}

TEST_CASE("Fermi-Dirac: finite temperature", "[tightbinding]") {
    Scalar kT = 0.025;  // ~300 K
    Scalar mu = 0.0;

    // At Fermi level, occupation should be 0.5
    REQUIRE(fermi_dirac(mu, mu, kT) == Approx(0.5).epsilon(1e-10));

    // Symmetric around Fermi level
    REQUIRE(fermi_dirac(mu - 0.1, mu, kT) + fermi_dirac(mu + 0.1, mu, kT) == Approx(1.0).epsilon(1e-10));
}

TEST_CASE("TBSolver: diagonal matrix", "[tightbinding]") {
    DenseHamiltonian ham;
    int n = 4;
    ham.resize(1, n);

    // Simple diagonal H with identity S
    ham.H = MatX::Zero(n, n);
    ham.S = MatX::Identity(n, n);

    ham.H(0, 0) = -1.0;
    ham.H(1, 1) = -0.5;
    ham.H(2, 2) = 0.0;
    ham.H(3, 3) = 0.5;

    TBSolver solver;
    SolverParams params;
    params.electronic_temperature = 0.01;
    solver.set_params(params);

    solver.solve(ham);

    // Check eigenvalues (should be same as diagonal elements, sorted)
    REQUIRE(ham.eigenvalues[0] == Approx(-1.0).epsilon(1e-10));
    REQUIRE(ham.eigenvalues[1] == Approx(-0.5).epsilon(1e-10));
    REQUIRE(ham.eigenvalues[2] == Approx(0.0).margin(1e-10));
    REQUIRE(ham.eigenvalues[3] == Approx(0.5).epsilon(1e-10));
}

TEST_CASE("TBSolver: occupation", "[tightbinding]") {
    DenseHamiltonian ham;
    int n = 4;
    ham.resize(1, n);

    ham.H = MatX::Zero(n, n);
    ham.S = MatX::Identity(n, n);

    ham.H(0, 0) = -1.0;
    ham.H(1, 1) = -0.5;
    ham.H(2, 2) = 0.5;
    ham.H(3, 3) = 1.0;

    TBSolver solver;
    SolverParams params;
    params.electronic_temperature = 0.001;  // Very low temperature
    solver.set_params(params);

    solver.solve(ham);

    // With 4 electrons (spin degeneracy 2), lowest 2 states should be filled
    solver.compute_occupation(ham, 4.0, 2);

    // Check total occupation
    Scalar total = ham.occupation.sum();
    REQUIRE(total == Approx(4.0).epsilon(1e-6));

    // First two states should be nearly fully occupied
    REQUIRE(ham.occupation[0] > 1.9);
    REQUIRE(ham.occupation[1] > 1.9);

    // Last two states should be nearly empty
    REQUIRE(ham.occupation[2] < 0.1);
    REQUIRE(ham.occupation[3] < 0.1);
}

TEST_CASE("TBSolver: band energy", "[tightbinding]") {
    DenseHamiltonian ham;
    int n = 4;
    ham.resize(1, n);

    ham.eigenvalues = VecX::Zero(n);
    ham.eigenvalues << -1.0, -0.5, 0.5, 1.0;

    ham.occupation = VecX::Zero(n);
    ham.occupation << 2.0, 2.0, 0.0, 0.0;  // Fill lowest 2 states

    TBSolver solver;
    Scalar E_band = solver.compute_band_energy(ham);

    // E_band = 2 * (-1.0) + 2 * (-0.5) = -3.0
    REQUIRE(E_band == Approx(-3.0).epsilon(1e-10));
}

TEST_CASE("SCCParams: defaults", "[tightbinding]") {
    SCCParams params;

    REQUIRE(params.max_iterations == 200);
    REQUIRE(params.mixing_parameter > 0.0);
    REQUIRE(params.mixing_parameter < 1.0);
}

TEST_CASE("Orbital mapping: full basis", "[tightbinding]") {
    // Full spd basis (9 orbitals)
    REQUIRE(get_absolute_orbital(9, 1) == 1);
    REQUIRE(get_absolute_orbital(9, 4) == 4);
    REQUIRE(get_absolute_orbital(9, 9) == 9);
}

TEST_CASE("Orbital mapping: sp basis", "[tightbinding]") {
    // sp basis (4 orbitals)
    REQUIRE(get_absolute_orbital(4, 1) == 1);
    REQUIRE(get_absolute_orbital(4, 2) == 2);
    REQUIRE(get_absolute_orbital(4, 4) == 4);
}

TEST_CASE("Orbital mapping: s basis", "[tightbinding]") {
    // s only basis (1 orbital)
    REQUIRE(get_absolute_orbital(1, 1) == 1);
}

TEST_CASE("Required SK integrals: s-s", "[tightbinding]") {
    auto integrals = get_required_integrals(1, 1);  // s-s

    // Should only need sss
    REQUIRE(integrals.size() == 1);
    REQUIRE(std::find(integrals.begin(), integrals.end(), 9) != integrals.end());
}

TEST_CASE("Required SK integrals: sp-sp", "[tightbinding]") {
    auto integrals = get_required_integrals(4, 4);  // sp-sp

    // Should need sss, sps, pps, ppp
    REQUIRE(integrals.size() >= 4);
    REQUIRE(std::find(integrals.begin(), integrals.end(), 9) != integrals.end());  // sss
    REQUIRE(std::find(integrals.begin(), integrals.end(), 8) != integrals.end());  // sps
    REQUIRE(std::find(integrals.begin(), integrals.end(), 5) != integrals.end());  // pps
    REQUIRE(std::find(integrals.begin(), integrals.end(), 6) != integrals.end());  // ppp
}
