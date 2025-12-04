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

#include "../config.hpp"
#include "../core/atomic_system.hpp"

namespace atomistica {

/**
 * @brief Pressure coupling modes for barostats
 */
enum class PressureMode {
    Isotropic,     // Same scaling in all directions (hydrostatic)
    Anisotropic,   // Independent scaling in x, y, z
    XY,            // Coupled x-y, free z
    XYZ_Fixed      // Fixed cell, compute pressure only
};

/**
 * @brief Compute stress tensor from forces and positions
 *
 * sigma = (1/V) * sum_i [m_i * v_i ⊗ v_i] + (1/V) * sum_{i<j} r_ij ⊗ f_ij
 *
 * The virial contribution from forces is:
 *   W = sum_i r_i ⊗ f_i
 *
 * @param system Atomic system
 * @param forces Forces on atoms
 * @param virial Optional precomputed virial tensor
 * @return 3x3 stress tensor
 */
inline Mat3 compute_stress_tensor(const AtomicSystem& system, const MatX3& forces,
                                  const Mat3* virial = nullptr) {
    int nat = system.num_atoms();
    Scalar volume = system.cell().determinant();

    // Kinetic contribution
    Mat3 stress = Mat3::Zero();
    for (int i = 0; i < nat; ++i) {
        Scalar m = system.mass(i);
        Vec3 v = system.velocity(i);
        stress += m * v * v.transpose();
    }

    // Virial contribution
    if (virial) {
        stress += *virial;
    } else {
        // Compute from positions and forces (site virial)
        for (int i = 0; i < nat; ++i) {
            Vec3 r = system.position(i);
            Vec3 f = forces.row(i).transpose();
            stress += r * f.transpose();
        }
    }

    // Convert to pressure units (eV/Å³ to GPa: 1 eV/Å³ = 160.2176634 GPa)
    return stress / volume;
}

/**
 * @brief Compute pressure from stress tensor
 *
 * P = (1/3) * Tr(sigma)
 */
inline Scalar compute_pressure(const Mat3& stress) {
    return stress.trace() / 3.0;
}

/**
 * @brief Berendsen barostat for pressure control
 *
 * Rescales the simulation cell according to:
 *   s = (1 + (dt/tau) * beta * (P_target - P))^(1/3)  (isotropic)
 *
 * or for each direction:
 *   s_ii = 1 + (dt/tau) * beta * (P_target_i - P_ii)  (anisotropic)
 *
 * The isothermal compressibility beta is typically 4.57e-5 bar^-1 for water.
 * In our units (eV, Å), we use beta in (Å³/eV).
 */
class BerendsenBarostat {
public:
    /**
     * @brief Constructor
     *
     * @param target_pressure Target pressure (eV/Å³)
     * @param tau Coupling time constant (fs)
     * @param beta Isothermal compressibility (Å³/eV)
     */
    BerendsenBarostat(Scalar target_pressure = 0.0, Scalar tau = 1000.0,
                      Scalar beta = 4.57e-5 * 160.2176634)  // Convert bar^-1 to Å³/eV
        : tau_(tau), beta_(beta) {
        P_target_[0] = P_target_[1] = P_target_[2] = target_pressure;
        tau_dir_[0] = tau_dir_[1] = tau_dir_[2] = tau;
    }

    /**
     * @brief Set target pressure (isotropic)
     */
    void set_pressure(Scalar P) {
        P_target_[0] = P_target_[1] = P_target_[2] = P;
    }

    /**
     * @brief Set target pressure (anisotropic)
     */
    void set_pressure(Scalar Px, Scalar Py, Scalar Pz) {
        P_target_[0] = Px;
        P_target_[1] = Py;
        P_target_[2] = Pz;
    }

    /**
     * @brief Set coupling time constant
     */
    void set_tau(Scalar tau) {
        tau_ = tau;
        tau_dir_[0] = tau_dir_[1] = tau_dir_[2] = tau;
    }

    /**
     * @brief Set coupling time constants per direction
     */
    void set_tau(Scalar tau_x, Scalar tau_y, Scalar tau_z) {
        tau_dir_[0] = tau_x;
        tau_dir_[1] = tau_y;
        tau_dir_[2] = tau_z;
    }

    /**
     * @brief Set pressure coupling mode
     */
    void set_mode(PressureMode mode) { mode_ = mode; }

    /**
     * @brief Apply barostat
     *
     * @param system Atomic system (positions and cell modified)
     * @param stress Current stress tensor
     * @param dt Time step
     */
    void apply(AtomicSystem& system, const Mat3& stress, Scalar dt) {
        Mat3 cell = system.cell();
        int nat = system.num_atoms();

        // Compute scaling matrix
        Mat3 scale = Mat3::Identity();

        switch (mode_) {
            case PressureMode::Isotropic: {
                // Hydrostatic pressure
                Scalar P = stress.trace() / 3.0;
                Scalar P_target = (P_target_[0] + P_target_[1] + P_target_[2]) / 3.0;
                Scalar s = std::cbrt(1.0 + dt * beta_ * (P_target - P) / tau_);
                scale = Mat3::Identity() * s;
                break;
            }

            case PressureMode::Anisotropic: {
                // Independent directions
                for (int d = 0; d < 3; ++d) {
                    if (tau_dir_[d] > 0.0) {
                        Scalar s = 1.0 + dt * beta_ * (P_target_[d] - stress(d, d)) / tau_dir_[d];
                        scale(d, d) = s;
                    }
                }
                break;
            }

            case PressureMode::XY: {
                // Coupled x-y
                Scalar P_xy = 0.5 * (stress(0, 0) + stress(1, 1));
                Scalar P_target_xy = 0.5 * (P_target_[0] + P_target_[1]);
                Scalar s = std::sqrt(1.0 + dt * beta_ * (P_target_xy - P_xy) / tau_);
                scale(0, 0) = s;
                scale(1, 1) = s;
                // z direction independent
                if (tau_dir_[2] > 0.0) {
                    scale(2, 2) = 1.0 + dt * beta_ * (P_target_[2] - stress(2, 2)) / tau_dir_[2];
                }
                break;
            }

            case PressureMode::XYZ_Fixed:
                // No scaling
                return;
        }

        // Scale cell vectors
        Mat3 new_cell = scale * cell;
        system.set_cell(new_cell);

        // Scale atomic positions
        for (int i = 0; i < nat; ++i) {
            Vec3 r = system.position(i);
            system.set_position(i, scale * r);
        }
    }

private:
    Scalar P_target_[3] = {0.0, 0.0, 0.0};  // Target pressure
    Scalar tau_ = 1000.0;                    // Coupling time (fs)
    Scalar tau_dir_[3] = {1000.0, 1000.0, 1000.0};  // Per-direction
    Scalar beta_ = 7.32e-3;                  // Compressibility (Å³/eV)
    PressureMode mode_ = PressureMode::Isotropic;
};

/**
 * @brief Andersen barostat (NPH ensemble)
 *
 * Extended system barostat with equations of motion:
 *   m_i * d²r_i/dt² = f_i - eta * m_i * v_i
 *   W * d²L/dt² = V * (P - P_target)
 *
 * where W is the fictitious barostat mass and eta = dL/dt / L.
 *
 * This generates a proper NPH ensemble with extended Hamiltonian:
 *   H_ext = K + U + (1/2)*W*eta² + P*V
 */
class AndersenBarostat {
public:
    /**
     * @brief Constructor
     *
     * @param target_pressure Target pressure (eV/Å³)
     * @param barostat_mass Fictitious mass for barostat
     */
    AndersenBarostat(Scalar target_pressure = 0.0, Scalar barostat_mass = 1.0)
        : P_target_(target_pressure), W_(barostat_mass) {
        eta_.setZero();
    }

    /**
     * @brief Set target pressure
     */
    void set_pressure(Scalar P) { P_target_ = P; }

    /**
     * @brief Set barostat mass
     */
    void set_mass(Scalar W) { W_ = W; }

    /**
     * @brief Apply barostat in step 1 (before force calculation)
     *
     * @param system Atomic system
     * @param stress Current stress tensor
     * @param dt Time step
     */
    void step1(AtomicSystem& system, const Mat3& stress, Scalar dt) {
        Mat3 cell = system.cell();
        int nat = system.num_atoms();

        // Cell lengths
        Vec3 L;
        for (int d = 0; d < 3; ++d) {
            L[d] = cell.col(d).norm();
        }
        Scalar V = cell.determinant();

        // Update barostat momentum
        Vec3 P_diag;
        for (int d = 0; d < 3; ++d) {
            P_diag[d] = stress(d, d);
        }
        eta_ += 0.5 * dt * (V / L.sum()) * (P_diag.array() - P_target_).matrix() / W_;

        // Update cell lengths
        Vec3 L_half = L + 0.5 * dt * eta_ / W_;
        Vec3 L_new = L_half + 0.5 * dt * eta_ / W_;

        // Scaling factors
        Vec3 scale = L_new.array() / L.array();

        // Scale cell
        Mat3 new_cell = cell;
        for (int d = 0; d < 3; ++d) {
            new_cell.col(d) *= scale[d];
        }
        system.set_cell(new_cell);

        // Scale positions and velocities
        Vec3 scale_v = L.array() / L_new.array();
        Vec3 scale_r = L_new.array() / L.array();

        for (int i = 0; i < nat; ++i) {
            Vec3 r = system.position(i);
            Vec3 v = system.velocity(i);

            // Scale each component
            for (int d = 0; d < 3; ++d) {
                r[d] *= scale_r[d];
                v[d] *= scale_v[d];
            }

            system.set_position(i, r);
            system.set_velocity(i, v);
        }
    }

    /**
     * @brief Apply barostat in step 2 (after force calculation)
     *
     * @param system Atomic system
     * @param stress Updated stress tensor
     * @param dt Time step
     */
    void step2(AtomicSystem& system, const Mat3& stress, Scalar dt) {
        Mat3 cell = system.cell();

        // Cell lengths
        Vec3 L;
        for (int d = 0; d < 3; ++d) {
            L[d] = cell.col(d).norm();
        }
        Scalar V = cell.determinant();

        // Update barostat momentum
        Vec3 P_diag;
        for (int d = 0; d < 3; ++d) {
            P_diag[d] = stress(d, d);
        }
        eta_ += 0.5 * dt * (V / L.sum()) * (P_diag.array() - P_target_).matrix() / W_;
    }

    /**
     * @brief Get barostat kinetic energy
     */
    Scalar kinetic_energy() const {
        return 0.5 * W_ * eta_.squaredNorm();
    }

    /**
     * @brief Get barostat energy contribution (kinetic + PV)
     */
    Scalar energy(const AtomicSystem& system) const {
        Scalar V = system.cell().determinant();
        return kinetic_energy() + P_target_ * V;
    }

private:
    Scalar P_target_ = 0.0;  // Target pressure
    Scalar W_ = 1.0;         // Barostat mass
    Vec3 eta_;               // Barostat momentum
};

/**
 * @brief Parrinello-Rahman barostat for fully flexible cell
 *
 * Allows the simulation cell shape to change, useful for
 * phase transitions and crystal simulations.
 *
 * Equations of motion include the full cell tensor h:
 *   W * d²h/dt² = V * (sigma - P_target * I) * (h^-1)^T
 */
class ParrinelloRahmanBarostat {
public:
    ParrinelloRahmanBarostat(Scalar target_pressure = 0.0, Scalar barostat_mass = 1.0)
        : P_target_(target_pressure), W_(barostat_mass) {
        h_dot_.setZero();
    }

    /**
     * @brief Set target pressure (isotropic)
     */
    void set_pressure(Scalar P) { P_target_ = P; }

    /**
     * @brief Set barostat mass
     */
    void set_mass(Scalar W) { W_ = W; }

    /**
     * @brief Apply barostat in step 1
     */
    void step1(AtomicSystem& system, const Mat3& stress, Scalar dt) {
        Mat3 h = system.cell();
        Scalar V = h.determinant();
        int nat = system.num_atoms();

        // Target stress tensor (isotropic)
        Mat3 sigma_target = Mat3::Identity() * P_target_;

        // Update cell velocity
        Mat3 h_inv_t = h.inverse().transpose();
        h_dot_ += 0.5 * dt * V * (stress - sigma_target) * h_inv_t / W_;

        // Update cell
        Mat3 h_new = h + dt * h_dot_;
        system.set_cell(h_new);

        // Update atomic positions (scaled coordinates preserved)
        Mat3 transform = h_new * h.inverse();
        for (int i = 0; i < nat; ++i) {
            Vec3 r = system.position(i);
            system.set_position(i, transform * r);
        }
    }

    /**
     * @brief Apply barostat in step 2
     */
    void step2(AtomicSystem& system, const Mat3& stress, Scalar dt) {
        Mat3 h = system.cell();
        Scalar V = h.determinant();

        Mat3 sigma_target = Mat3::Identity() * P_target_;
        Mat3 h_inv_t = h.inverse().transpose();

        h_dot_ += 0.5 * dt * V * (stress - sigma_target) * h_inv_t / W_;
    }

    /**
     * @brief Get barostat kinetic energy
     */
    Scalar kinetic_energy() const {
        return 0.5 * W_ * h_dot_.squaredNorm();
    }

private:
    Scalar P_target_ = 0.0;
    Scalar W_ = 1.0;
    Mat3 h_dot_;  // Cell velocity
};

} // namespace atomistica
