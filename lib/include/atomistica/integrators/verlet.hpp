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

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "../config.hpp"
#include "../core/atomic_system.hpp"

namespace atomistica {

/**
 * @brief Velocity Verlet integrator for molecular dynamics
 *
 * Implements the two-step velocity Verlet algorithm:
 *
 * Step 1 (before force calculation):
 *   v(t+dt/2) = v(t) + 0.5 * f(t)/m * dt
 *   r(t+dt) = r(t) + v(t+dt/2) * dt
 *
 * Step 2 (after force calculation):
 *   v(t+dt) = v(t+dt/2) + 0.5 * f(t+dt)/m * dt
 *
 * This algorithm is time-reversible and symplectic, providing excellent
 * energy conservation for constant energy (NVE) simulations.
 */
class VelocityVerlet {
public:
    VelocityVerlet() = default;

    /**
     * @brief Set time step
     *
     * @param dt Time step in femtoseconds
     */
    void set_timestep(Scalar dt) { dt_ = dt; }

    /**
     * @brief Get time step
     */
    Scalar timestep() const { return dt_; }

    /**
     * @brief Set maximum displacement for adaptive time stepping
     *
     * @param max_dr Maximum allowed displacement per step (0 = disabled)
     */
    void set_max_displacement(Scalar max_dr) { max_dr_ = max_dr; }

    /**
     * @brief Perform step 1: update positions and half-step velocities
     *
     * @param system Atomic system (positions and velocities modified)
     * @param forces Current forces
     * @return Maximum displacement of any atom
     */
    Scalar step1(AtomicSystem& system, const MatX3& forces) {
        int nat = system.num_atoms();
        Scalar max_disp_sq = 0.0;

        // Adaptive time stepping
        Scalar dt = dt_;
        if (max_dr_ > 0.0) {
            dt = compute_adaptive_dt(system, forces);
        }
        dt_actual_ = dt;

        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            Scalar inv_m = 1.0 / m;

            Vec3 v = system.velocity(i);
            Vec3 f = forces.row(i).transpose();

            // Half-step velocity update
            v += 0.5 * inv_m * f * dt;

            // Position update
            Vec3 dr = v * dt;
            Vec3 r = system.position(i) + dr;

            // Track maximum displacement
            max_disp_sq = std::max(max_disp_sq, dr.squaredNorm());

            system.set_velocity(i, v);
            system.set_position(i, r);
        }

        return std::sqrt(max_disp_sq);
    }

    /**
     * @brief Perform step 2: complete velocity update
     *
     * @param system Atomic system (velocities modified)
     * @param forces Updated forces at new positions
     */
    void step2(AtomicSystem& system, const MatX3& forces) {
        int nat = system.num_atoms();
        Scalar dt = dt_actual_;

        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            Scalar inv_m = 1.0 / m;

            Vec3 v = system.velocity(i);
            Vec3 f = forces.row(i).transpose();

            // Complete velocity update
            v += 0.5 * inv_m * f * dt;

            system.set_velocity(i, v);
        }
    }

    /**
     * @brief Get kinetic energy
     */
    static Scalar kinetic_energy(const AtomicSystem& system) {
        Scalar E_kin = 0.0;
        for (int i = 0; i < system.num_atoms(); ++i) {
            Scalar m = system.mass(i);
            Vec3 v = system.velocity(i);
            E_kin += 0.5 * m * v.squaredNorm();
        }
        return E_kin;
    }

    /**
     * @brief Get instantaneous temperature
     *
     * T = 2 * E_kin / (N_dof * k_B)
     */
    static Scalar temperature(const AtomicSystem& system, int n_constraints = 0) {
        Scalar E_kin = kinetic_energy(system);
        int n_dof = 3 * system.num_atoms() - n_constraints;
        if (n_dof <= 0) return 0.0;

        // k_B in eV/K
        constexpr Scalar kB = 8.617333262e-5;
        return 2.0 * E_kin / (n_dof * kB);
    }

    /**
     * @brief Get actual time step used (may differ due to adaptive stepping)
     */
    Scalar actual_timestep() const { return dt_actual_; }

private:
    Scalar dt_ = 1.0;         // Time step (fs)
    Scalar max_dr_ = 0.0;     // Max displacement for adaptive dt (0 = disabled)
    Scalar dt_actual_ = 1.0;  // Actual dt used (may be adapted)

    /**
     * @brief Compute adaptive time step based on maximum displacement
     */
    Scalar compute_adaptive_dt(const AtomicSystem& system, const MatX3& forces) {
        Scalar max_dr_sq = 0.0;

        for (int i = 0; i < system.num_atoms(); ++i) {
            Scalar m = system.mass(i);
            Scalar inv_m = 1.0 / m;
            Vec3 v = system.velocity(i);
            Vec3 f = forces.row(i).transpose();

            // Estimate displacement at maximum time step
            Vec3 dr = v * dt_ + 0.5 * inv_m * f * dt_ * dt_;
            max_dr_sq = std::max(max_dr_sq, dr.squaredNorm());
        }

        if (max_dr_sq < 1e-20) return dt_;

        // Scale time step to limit displacement
        Scalar scale = max_dr_ / std::sqrt(max_dr_sq);
        return std::min(scale * dt_, dt_);
    }
};

} // namespace atomistica
