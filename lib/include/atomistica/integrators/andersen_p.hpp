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

#include "../config.hpp"
#include "../core/atomic_system.hpp"

namespace atomistica {

/**
 * @brief Andersen pressure controller (NPH ensemble)
 *
 * Implements Andersen's extended system method for constant pressure MD.
 * See: H.C. Andersen, J. Chem. Phys. 72, 2384 (1980)
 *
 * The method introduces a fictitious "piston" mass W that controls how
 * quickly the simulation box can change volume in response to pressure
 * differences.
 */
class AndersenP {
public:
    AndersenP() = default;

    /**
     * @brief Set target pressure (3 components for anisotropic pressure)
     * @param Px Target pressure in x direction
     * @param Py Target pressure in y direction
     * @param Pz Target pressure in z direction
     */
    void set_target_pressure(Scalar Px, Scalar Py, Scalar Pz) {
        P_target_[0] = Px;
        P_target_[1] = Py;
        P_target_[2] = Pz;
    }

    /**
     * @brief Set target pressure (isotropic)
     * @param P Target pressure
     */
    void set_target_pressure(Scalar P) {
        P_target_[0] = P_target_[1] = P_target_[2] = P;
    }

    /**
     * @brief Set barostat mass
     * @param W Fictitious barostat mass (larger = slower box dynamics)
     */
    void set_barostat_mass(Scalar W) {
        W_ = W;
    }

    /**
     * @brief Set timestep
     * @param dt Timestep in femtoseconds
     */
    void set_timestep(Scalar dt) {
        dt_ = dt;
    }

    /**
     * @brief Initialize the barostat
     * @param system The atomic system
     */
    void initialize(AtomicSystem& system);

    /**
     * @brief First half of integration step
     *
     * Updates:
     * - eta (barostat momentum) using half-step
     * - velocities using full step
     * - positions using full step with cell scaling
     * - cell dimensions
     *
     * @param system The atomic system
     * @param forces Current forces on atoms (N x 3 matrix)
     * @param pressure Current pressure tensor (diagonal)
     */
    void step1(AtomicSystem& system, const MatX3& forces, const Vec3& pressure);

    /**
     * @brief Second half of integration step
     *
     * Updates:
     * - velocities using half-step
     * - eta (barostat momentum) using half-step
     *
     * @param system The atomic system
     * @param forces Current forces on atoms (N x 3 matrix)
     * @param pressure Updated pressure tensor (diagonal)
     */
    void step2(AtomicSystem& system, const MatX3& forces, const Vec3& pressure);

    /**
     * @brief Get barostat kinetic energy
     * @return Kinetic energy of the barostat degree of freedom
     */
    Scalar barostat_energy() const {
        return 0.5 * W_ * eta_.squaredNorm();
    }

    /**
     * @brief Get current eta (barostat velocity)
     */
    const Vec3& eta() const { return eta_; }

    /**
     * @brief Get current cell dimensions
     */
    Vec3 cell_lengths() const { return L_; }

private:
    Scalar dt_ = 0.1;           // Timestep
    Scalar W_ = 1e4;            // Barostat mass
    Vec3 P_target_ = Vec3::Zero();  // Target pressure (diagonal components)
    Vec3 eta_ = Vec3::Zero();   // Barostat velocity (dL/dt / some_factor)
    Vec3 L_ = Vec3::Zero();     // Current cell dimensions
    bool initialized_ = false;
};

// Implementation

inline void AndersenP::initialize(AtomicSystem& system) {
    // Get initial cell dimensions (assumes orthorhombic cell)
    L_(0) = system.cell()(0, 0);
    L_(1) = system.cell()(1, 1);
    L_(2) = system.cell()(2, 2);

    // Initialize eta to zero
    eta_.setZero();

    initialized_ = true;
}

inline void AndersenP::step1(AtomicSystem& system, const MatX3& forces, const Vec3& pressure) {
    if (!initialized_) {
        initialize(system);
    }

    Scalar dt2 = dt_ / 2.0;

    // Volume and cell dimensions
    Scalar V = L_(0) * L_(1) * L_(2);

    // Update barostat momentum (half step)
    // eta = eta + dt/2 * V/L * (P - P_target)
    for (int d = 0; d < 3; ++d) {
        eta_(d) += dt2 * (V / L_(d)) * (pressure(d) - P_target_(d));
    }

    // Update velocities (full step): v = v + dt * f / m
    auto& velocities = system.velocities();
    auto& masses = system.masses();
    int nat = system.num_atoms();

    for (int i = 0; i < nat; ++i) {
        Scalar m_inv = 1.0 / masses(i);
        for (int d = 0; d < 3; ++d) {
            velocities(d, i) += dt_ * forces(i, d) * m_inv;
        }
    }

    // Predict new cell dimensions
    // L_dt2 = L + dt/2 * eta / W
    // L_dt = L + dt * eta / W
    Vec3 L_dt2 = L_ + (dt_ / 2.0) * eta_ / W_;
    Vec3 L_dt = L_ + dt_ * eta_ / W_;

    // Position update with scaling factor
    // fac = (L / L_dt2)^2
    auto& positions = system.positions();
    for (int i = 0; i < nat; ++i) {
        for (int d = 0; d < 3; ++d) {
            Scalar fac = (L_(d) / L_dt2(d)) * (L_(d) / L_dt2(d));
            // r_new = r + dt * v * fac
            positions(d, i) += dt_ * velocities(d, i) * fac;
        }
    }

    // Scale velocities and positions to new cell
    // vfac = L / L_dt  (velocities scale inversely with cell)
    // rfac = L_dt / L  (positions scale with cell)
    for (int i = 0; i < nat; ++i) {
        for (int d = 0; d < 3; ++d) {
            Scalar vfac = L_(d) / L_dt(d);
            Scalar rfac = L_dt(d) / L_(d);
            velocities(d, i) *= vfac;
            positions(d, i) *= rfac;
        }
    }

    // Update cell
    Mat3 new_cell = Mat3::Zero();
    new_cell(0, 0) = L_dt(0);
    new_cell(1, 1) = L_dt(1);
    new_cell(2, 2) = L_dt(2);
    system.set_cell(new_cell);

    L_ = L_dt;
}

inline void AndersenP::step2(AtomicSystem& system, const MatX3& forces, const Vec3& pressure) {
    Scalar dt2 = dt_ / 2.0;

    // Update velocities (half step): v = v + dt/2 * f / m
    auto& velocities = system.velocities();
    auto& masses = system.masses();
    int nat = system.num_atoms();

    for (int i = 0; i < nat; ++i) {
        Scalar m_inv = 1.0 / masses(i);
        for (int d = 0; d < 3; ++d) {
            velocities(d, i) += dt2 * forces(i, d) * m_inv;
        }
    }

    // Update barostat momentum (half step)
    Scalar V = L_(0) * L_(1) * L_(2);
    for (int d = 0; d < 3; ++d) {
        eta_(d) += dt2 * (V / L_(d)) * (pressure(d) - P_target_(d));
    }
}

} // namespace atomistica
