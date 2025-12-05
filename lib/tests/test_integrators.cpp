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

#include <atomistica/integrators/integrators.hpp>
#include <atomistica/core/atomic_system.hpp>

using namespace atomistica;
using Catch::Approx;

// Helper to create a simple 2-atom system
AtomicSystem create_two_atom_system() {
    AtomicSystem system;
    system.set_cell(Mat3::Identity() * 10.0);

    // Two atoms
    system.add_atom(6, Vec3(4.0, 5.0, 5.0), 12.0);  // Carbon
    system.add_atom(6, Vec3(6.0, 5.0, 5.0), 12.0);  // Carbon

    return system;
}

TEST_CASE("VelocityVerlet: free particle motion", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(0.1, 0.0, 0.0));
    system.set_velocity(1, Vec3(-0.1, 0.0, 0.0));

    VelocityVerlet verlet;
    verlet.set_timestep(1.0);

    MatX3 forces = MatX3::Zero(2, 3);

    Vec3 r0_initial = system.position(0);
    Vec3 v0_initial = system.velocity(0);

    // Step 1
    verlet.step1(system, forces);
    verlet.step2(system, forces);

    // Position should change by v * dt
    Vec3 r0_expected = r0_initial + v0_initial * 1.0;
    REQUIRE((system.position(0) - r0_expected).norm() == Approx(0.0).margin(1e-10));

    // Velocity should be unchanged
    REQUIRE((system.velocity(0) - v0_initial).norm() == Approx(0.0).margin(1e-10));
}

TEST_CASE("VelocityVerlet: constant force motion", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(0.1, 0.0, 0.0));
    system.set_velocity(1, Vec3(-0.1, 0.0, 0.0));

    VelocityVerlet verlet;
    Scalar dt = 0.1;
    verlet.set_timestep(dt);

    MatX3 forces = MatX3::Zero(2, 3);
    forces(0, 0) = 1.0;  // Force on atom 0 in x direction

    Vec3 r0_initial = system.position(0);
    Vec3 v0_initial = system.velocity(0);
    Scalar m = system.mass(0);

    // Do 10 steps
    for (int step = 0; step < 10; ++step) {
        verlet.step1(system, forces);
        verlet.step2(system, forces);
    }

    // Expected position: r = r0 + v0*t + 0.5*a*t^2
    Scalar t = 10.0 * dt;
    Vec3 r0_expected = r0_initial + v0_initial * t + 0.5 * (forces.row(0).transpose() / m) * t * t;

    REQUIRE((system.position(0) - r0_expected).norm() == Approx(0.0).margin(1e-8));
}

TEST_CASE("VelocityVerlet: kinetic energy", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(0.1, 0.0, 0.0));
    system.set_velocity(1, Vec3(-0.1, 0.0, 0.0));

    Scalar E_kin = VelocityVerlet::kinetic_energy(system);

    // E = 0.5 * m * v^2 for each atom
    Scalar m = system.mass(0);
    Scalar v = 0.1;
    Scalar expected = 2 * 0.5 * m * v * v;

    REQUIRE(E_kin == Approx(expected).epsilon(1e-10));
}

TEST_CASE("VelocityVerlet: temperature calculation", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(0.1, 0.0, 0.0));
    system.set_velocity(1, Vec3(-0.1, 0.0, 0.0));

    Scalar T = VelocityVerlet::temperature(system, 3);
    REQUIRE(T > 0.0);
}

TEST_CASE("BerendsenThermostat: construction", "[integrators]") {
    BerendsenThermostat thermo(300.0, 500.0);

    REQUIRE(thermo.target_temperature() == Approx(300.0).epsilon(1e-10));
    REQUIRE(thermo.tau() == Approx(500.0).epsilon(1e-10));
}

TEST_CASE("BerendsenThermostat: temperature control", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(1.0, 0.5, 0.0));
    system.set_velocity(1, Vec3(-1.0, -0.5, 0.0));

    BerendsenThermostat thermo(300.0, 100.0);

    Scalar T_initial = VelocityVerlet::temperature(system, 3);

    // Apply thermostat several times
    for (int i = 0; i < 100; ++i) {
        thermo.apply(system, 1.0, 3);
    }

    Scalar T_final = VelocityVerlet::temperature(system, 3);

    // Temperature should move toward target
    REQUIRE(std::abs(T_final - 300.0) < std::abs(T_initial - 300.0));
}

TEST_CASE("BerendsenThermostat: instant rescaling", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(1.0, 0.5, 0.0));
    system.set_velocity(1, Vec3(-1.0, -0.5, 0.0));

    // tau = 0 should give instant rescaling
    BerendsenThermostat thermo(300.0, 0.0);

    thermo.apply(system, 1.0, 3);

    Scalar T = VelocityVerlet::temperature(system, 3);
    REQUIRE(T == Approx(300.0).epsilon(0.01));  // Should be very close to target
}

TEST_CASE("LangevinThermostat: construction", "[integrators]") {
    LangevinThermostat langevin(300.0, 0.01, 12345);
    REQUIRE(true);  // Just verify construction doesn't crash
}

TEST_CASE("LangevinThermostat: multiple steps", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.set_velocity(0, Vec3(0.1, 0.0, 0.0));
    system.set_velocity(1, Vec3(-0.1, 0.0, 0.0));

    LangevinThermostat langevin(300.0, 0.01, 12345);
    MatX3 forces = MatX3::Zero(2, 3);

    // Run several steps
    for (int i = 0; i < 100; ++i) {
        langevin.step1(system, forces, 1.0);
        langevin.step2(system, forces, 1.0);
    }

    // Verify positions are still finite
    REQUIRE(system.position(0).allFinite());
    REQUIRE(system.position(1).allFinite());
}

TEST_CASE("NoseHooverThermostat: construction", "[integrators]") {
    NoseHooverThermostat nh(300.0, 100.0, 3);
    nh.init(6);
    REQUIRE(true);  // Just verify construction doesn't crash
}

TEST_CASE("Velocity initialization: zero mean velocity", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.add_atom(6, Vec3(5.0, 5.0, 6.0), 12.0);
    system.add_atom(6, Vec3(5.0, 5.0, 4.0), 12.0);

    initialize_velocities(system, 300.0, 12345, true);

    // With COM removal, mean velocity should be zero
    Vec3 v_mean = Vec3::Zero();
    Scalar total_mass = 0.0;
    for (int i = 0; i < system.num_atoms(); ++i) {
        Scalar m = system.mass(i);
        v_mean += m * system.velocity(i);
        total_mass += m;
    }
    v_mean /= total_mass;

    REQUIRE(v_mean.norm() == Approx(0.0).margin(1e-10));
}

TEST_CASE("Velocity initialization: target temperature", "[integrators]") {
    AtomicSystem system = create_two_atom_system();
    system.add_atom(6, Vec3(5.0, 5.0, 6.0), 12.0);
    system.add_atom(6, Vec3(5.0, 5.0, 4.0), 12.0);

    Scalar target_T = 300.0;
    initialize_velocities(system, target_T, 12345, true);

    // Temperature should be exactly target (due to rescaling)
    Scalar T = VelocityVerlet::temperature(system, 3);  // 3 constraints for COM
    REQUIRE(T == Approx(target_T).epsilon(0.01));  // Within 1%
}

TEST_CASE("BerendsenBarostat: construction", "[integrators]") {
    BerendsenBarostat baro(0.0, 1000.0);
    REQUIRE(true);
}

TEST_CASE("BerendsenBarostat: isotropic mode", "[integrators]") {
    BerendsenBarostat baro(0.0, 1000.0);
    baro.set_mode(PressureMode::Isotropic);

    AtomicSystem system = create_two_atom_system();

    // Create a stress tensor (isotropic, slightly above target)
    Mat3 stress = Mat3::Identity() * 0.001;  // Positive pressure

    Mat3 cell_initial = system.cell();

    // Apply barostat
    baro.apply(system, stress, 1.0);

    // Cell should change
    Mat3 cell_final = system.cell();
    REQUIRE(cell_final.determinant() != cell_initial.determinant());
}

TEST_CASE("AndersenBarostat: construction", "[integrators]") {
    AndersenBarostat baro(0.0, 1.0);
    REQUIRE(true);
}

TEST_CASE("AndersenBarostat: energy contribution", "[integrators]") {
    AndersenBarostat baro(0.001, 1.0);  // Small positive pressure

    AtomicSystem system = create_two_atom_system();

    Scalar E = baro.energy(system);

    // Energy should be positive (P * V with P > 0)
    REQUIRE(E > 0.0);
}

TEST_CASE("Stress tensor: kinetic contribution", "[integrators]") {
    AtomicSystem system = create_two_atom_system();

    // Set velocities
    system.set_velocity(0, Vec3(1.0, 0.0, 0.0));
    system.set_velocity(1, Vec3(-1.0, 0.0, 0.0));

    MatX3 forces = MatX3::Zero(2, 3);

    Mat3 stress = compute_stress_tensor(system, forces);

    // Stress should be non-zero due to kinetic contribution
    Scalar P = compute_pressure(stress);
    REQUIRE(P != 0.0);
}

TEST_CASE("ParrinelloRahmanBarostat: construction", "[integrators]") {
    ParrinelloRahmanBarostat pr(0.0, 1.0);
    REQUIRE(true);
}

TEST_CASE("NVE integration: energy conservation", "[integrators]") {
    AtomicSystem system = create_two_atom_system();

    // Set initial velocities
    system.set_velocity(0, Vec3(0.01, 0.0, 0.0));
    system.set_velocity(1, Vec3(-0.01, 0.0, 0.0));

    VelocityVerlet verlet;
    verlet.set_timestep(0.5);

    // Simple harmonic force between atoms
    auto compute_forces = [&](MatX3& forces) {
        Vec3 r01 = system.position(1) - system.position(0);
        Scalar r = r01.norm();
        Scalar r_eq = 2.0;
        Scalar k = 0.1;

        Vec3 f = -k * (r - r_eq) * r01.normalized();
        forces.row(0) = -f.transpose();
        forces.row(1) = f.transpose();

        // Return potential energy
        return 0.5 * k * (r - r_eq) * (r - r_eq);
    };

    MatX3 forces = MatX3::Zero(2, 3);
    Scalar U = compute_forces(forces);
    Scalar E_initial = VelocityVerlet::kinetic_energy(system) + U;

    // Run 1000 steps
    for (int step = 0; step < 1000; ++step) {
        verlet.step1(system, forces);
        U = compute_forces(forces);
        verlet.step2(system, forces);
    }

    Scalar E_final = VelocityVerlet::kinetic_energy(system) + U;

    // Energy should be conserved to within 1%
    REQUIRE(E_final == Approx(E_initial).epsilon(0.01));
}
