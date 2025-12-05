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
#include <random>
#include <vector>

#include "../config.hpp"
#include "../core/atomic_system.hpp"

namespace atomistica {

// Boltzmann constant in eV/K
constexpr Scalar kB_eV_K = 8.617333262e-5;

/**
 * @brief Berendsen thermostat for temperature control
 *
 * Rescales velocities according to:
 *   v_new = v * sqrt(1 + (dt/tau) * (T_target/T - 1))
 *
 * This is a weak coupling thermostat that drives the system towards
 * the target temperature with characteristic time tau.
 *
 * Note: Does not generate a proper canonical ensemble, but is simple
 * and effective for equilibration.
 */
class BerendsenThermostat {
public:
    /**
     * @brief Constructor
     *
     * @param target_temperature Target temperature in Kelvin
     * @param tau Coupling time constant in femtoseconds
     */
    BerendsenThermostat(Scalar target_temperature = 300.0, Scalar tau = 500.0)
        : T_target_(target_temperature), tau_(tau) {}

    /**
     * @brief Set target temperature
     */
    void set_temperature(Scalar T) { T_target_ = T; }

    /**
     * @brief Get target temperature
     */
    Scalar target_temperature() const { return T_target_; }

    /**
     * @brief Set coupling time constant
     */
    void set_tau(Scalar tau) { tau_ = tau; }

    /**
     * @brief Get coupling time constant
     */
    Scalar tau() const { return tau_; }

    /**
     * @brief Set temperature rate of change (for cooling/heating)
     *
     * @param dT Temperature change per femtosecond
     */
    void set_temperature_rate(Scalar dT) { dT_ = dT; }

    /**
     * @brief Enable exponential cooling to final temperature
     *
     * @param T_final Final temperature
     */
    void set_exponential_cooling(Scalar T_final) { T_final_ = T_final; }

    /**
     * @brief Apply thermostat to system
     *
     * @param system Atomic system (velocities modified)
     * @param dt Time step
     * @param n_constraints Number of constraints (for DOF calculation)
     */
    void apply(AtomicSystem& system, Scalar dt, int n_constraints = 0) {
        // Get current temperature
        Scalar E_kin = 0.0;
        int nat = system.num_atoms();
        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            Vec3 v = system.velocity(i);
            E_kin += 0.5 * m * v.squaredNorm();
        }

        int n_dof = 3 * nat - n_constraints;
        if (n_dof <= 0) return;

        Scalar T_current = 2.0 * E_kin / (n_dof * kB_eV_K);
        if (T_current < 1e-10) return;

        // Compute scaling factor
        Scalar scale;
        if (tau_ <= 0.0) {
            // Instant rescaling
            scale = std::sqrt(T_target_ / T_current);
        } else {
            // Berendsen coupling
            scale = std::sqrt(1.0 + (dt / tau_) * (T_target_ / T_current - 1.0));
        }

        // Rescale velocities
        for (int i = 0; i < nat; ++i) {
            Vec3 v = system.velocity(i);
            system.set_velocity(i, v * scale);
        }

        // Update target temperature if cooling/heating
        if (std::abs(dT_) > 1e-20) {
            T_target_ += dT_ * dt;
            T_target_ = std::max(T_target_, 0.0);
        }

        // Exponential cooling
        if (T_final_ > 0.0 && tau_ > 0.0) {
            T_target_ = T_final_ + (T_target_ - T_final_) * std::exp(-dt / tau_);
        }
    }

private:
    Scalar T_target_ = 300.0;  // Target temperature (K)
    Scalar tau_ = 500.0;       // Coupling time (fs)
    Scalar dT_ = 0.0;          // Temperature rate (K/fs)
    Scalar T_final_ = 0.0;     // Final temperature for exponential cooling
};

/**
 * @brief Langevin thermostat for temperature control
 *
 * Applies friction and random forces to maintain temperature:
 *   m*a = F - gamma*m*v + R
 *
 * where R is a random force satisfying the fluctuation-dissipation theorem:
 *   <R(t)*R(t')> = 2*gamma*m*kT*delta(t-t')
 *
 * This thermostat generates a proper canonical (NVT) ensemble.
 */
class LangevinThermostat {
public:
    /**
     * @brief Constructor
     *
     * @param target_temperature Target temperature in Kelvin
     * @param friction Friction coefficient (1/fs)
     * @param seed Random seed
     */
    LangevinThermostat(Scalar target_temperature = 300.0, Scalar friction = 0.01,
                       unsigned int seed = 12345)
        : T_target_(target_temperature), gamma_(friction), rng_(seed) {}

    /**
     * @brief Set target temperature
     */
    void set_temperature(Scalar T) { T_target_ = T; }

    /**
     * @brief Set friction coefficient
     */
    void set_friction(Scalar gamma) { gamma_ = gamma; }

    /**
     * @brief Set random seed
     */
    void set_seed(unsigned int seed) { rng_.seed(seed); }

    /**
     * @brief Compute coefficients for integration
     *
     * Precomputes coefficients c0, c1, c2 for the BBK integration scheme.
     */
    void compute_coefficients(Scalar dt) {
        Scalar gamma_dt = gamma_ * dt;

        if (gamma_dt < 1e-8) {
            // No friction - standard Verlet
            c0_ = 1.0;
            c1_ = 1.0;
            c2_ = 0.5;
        } else {
            c0_ = std::exp(-gamma_dt);
            c1_ = (1.0 - c0_) / gamma_dt;
            c2_ = (1.0 - c1_) / gamma_dt;
        }
    }

    /**
     * @brief Apply thermostat in step 1 (position update)
     *
     * @param system Atomic system
     * @param forces Current forces
     * @param dt Time step
     * @return Maximum displacement
     */
    Scalar step1(AtomicSystem& system, const MatX3& forces, Scalar dt) {
        compute_coefficients(dt);

        int nat = system.num_atoms();
        Scalar max_disp_sq = 0.0;

        // Standard normal distribution
        std::normal_distribution<Scalar> normal(0.0, 1.0);

        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            Scalar inv_m = 1.0 / m;
            Vec3 v = system.velocity(i);
            Vec3 f = forces.row(i).transpose();

            // Compute random force amplitudes
            Scalar gamma_dt = gamma_ * dt;
            Scalar sigma_r = 0.0, sigma_v = 0.0;

            if (gamma_dt > 1e-8 && T_target_ > 1e-10) {
                Scalar hlp = 2.0 - (3.0 + c0_*c0_ - 4.0*c0_) / gamma_dt;
                if (hlp > 0.0) {
                    sigma_r = std::sqrt(kB_eV_K * T_target_ / m * dt * dt / gamma_dt * hlp);
                }
                sigma_v = std::sqrt(kB_eV_K * T_target_ / m * (1.0 - c0_*c0_));
            }

            // Position update with random displacement
            Vec3 dr = c1_ * v * dt + c2_ * inv_m * f * dt * dt;

            // Add correlated random forces
            if (sigma_r > 0.0) {
                for (int d = 0; d < 3; ++d) {
                    dr[d] += sigma_r * normal(rng_);
                }
            }

            Vec3 r_new = system.position(i) + dr;
            max_disp_sq = std::max(max_disp_sq, dr.squaredNorm());

            // Velocity update
            Vec3 v_new = c0_ * v + (c1_ - c2_) * inv_m * f * dt;

            if (sigma_v > 0.0) {
                for (int d = 0; d < 3; ++d) {
                    v_new[d] += sigma_v * normal(rng_);
                }
            }

            system.set_position(i, r_new);
            system.set_velocity(i, v_new);
        }

        return std::sqrt(max_disp_sq);
    }

    /**
     * @brief Apply thermostat in step 2 (velocity completion)
     *
     * @param system Atomic system
     * @param forces Updated forces
     * @param dt Time step
     */
    void step2(AtomicSystem& system, const MatX3& forces, Scalar dt) {
        int nat = system.num_atoms();

        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            Scalar inv_m = 1.0 / m;
            Vec3 v = system.velocity(i);
            Vec3 f = forces.row(i).transpose();

            // Complete velocity update
            v += c2_ * inv_m * f * dt;
            system.set_velocity(i, v);
        }
    }

private:
    Scalar T_target_ = 300.0;  // Target temperature (K)
    Scalar gamma_ = 0.01;      // Friction coefficient (1/fs)

    // Integration coefficients
    Scalar c0_ = 1.0;
    Scalar c1_ = 1.0;
    Scalar c2_ = 0.5;

    std::mt19937 rng_;
};

/**
 * @brief Nose-Hoover thermostat (chain) for temperature control
 *
 * Extended system thermostat that generates a proper canonical ensemble.
 * Uses a chain of thermostats for better ergodicity.
 */
class NoseHooverThermostat {
public:
    /**
     * @brief Constructor
     *
     * @param target_temperature Target temperature in Kelvin
     * @param tau Coupling time constant
     * @param chain_length Number of thermostats in chain
     */
    NoseHooverThermostat(Scalar target_temperature = 300.0, Scalar tau = 100.0,
                         int chain_length = 3)
        : T_target_(target_temperature), tau_(tau) {
        eta_.resize(chain_length, 0.0);
        eta_dot_.resize(chain_length, 0.0);
        Q_.resize(chain_length, 0.0);
    }

    /**
     * @brief Initialize thermostat masses
     */
    void init(int n_dof) {
        n_dof_ = n_dof;
        Scalar kT = kB_eV_K * T_target_;

        // Thermostat mass Q = N_dof * k_B * T * tau^2
        Q_[0] = n_dof * kT * tau_ * tau_;
        for (size_t k = 1; k < Q_.size(); ++k) {
            Q_[k] = kT * tau_ * tau_;
        }
    }

    /**
     * @brief Apply thermostat (integrate thermostat chain)
     *
     * @param system Atomic system
     * @param dt Time step
     * @param n_constraints Number of constraints
     */
    void apply(AtomicSystem& system, Scalar dt, int n_constraints = 0) {
        int nat = system.num_atoms();
        n_dof_ = 3 * nat - n_constraints;
        if (n_dof_ <= 0) return;

        Scalar kT = kB_eV_K * T_target_;

        // Compute kinetic energy
        Scalar E_kin = 0.0;
        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            Vec3 v = system.velocity(i);
            E_kin += 0.5 * m * v.squaredNorm();
        }

        // Update thermostat chain (RESPA-style multiple time step)
        int n_chain = eta_.size();
        Scalar dt_chain = dt / 4.0;  // Sub-step

        for (int step = 0; step < 4; ++step) {
            // Update last thermostat
            if (n_chain > 1) {
                eta_dot_[n_chain-1] += dt_chain * (Q_[n_chain-2] * eta_dot_[n_chain-2] * eta_dot_[n_chain-2] - kT) / Q_[n_chain-1];
            }

            // Update middle thermostats
            for (int k = n_chain - 2; k > 0; --k) {
                Scalar factor = std::exp(-dt_chain * eta_dot_[k+1] / 2.0);
                eta_dot_[k] = eta_dot_[k] * factor;
                eta_dot_[k] += dt_chain * (Q_[k-1] * eta_dot_[k-1] * eta_dot_[k-1] - kT) / Q_[k];
                eta_dot_[k] = eta_dot_[k] * factor;
            }

            // Update first thermostat
            Scalar factor = std::exp(-dt_chain * eta_dot_[1] / 2.0);
            eta_dot_[0] = eta_dot_[0] * factor;
            eta_dot_[0] += dt_chain * (2.0 * E_kin - n_dof_ * kT) / Q_[0];
            eta_dot_[0] = eta_dot_[0] * factor;

            // Update thermostat positions
            for (int k = 0; k < n_chain; ++k) {
                eta_[k] += dt_chain * eta_dot_[k];
            }

            // Scale particle velocities
            Scalar scale = std::exp(-dt_chain * eta_dot_[0]);
            for (int i = 0; i < nat; ++i) {
                Vec3 v = system.velocity(i);
                system.set_velocity(i, v * scale);
            }

            // Update kinetic energy after scaling
            E_kin *= scale * scale;
        }
    }

    /**
     * @brief Get thermostat energy contribution
     */
    Scalar energy() const {
        Scalar E_therm = 0.0;
        Scalar kT = kB_eV_K * T_target_;

        // Kinetic energy of thermostat chain
        for (size_t k = 0; k < Q_.size(); ++k) {
            E_therm += 0.5 * Q_[k] * eta_dot_[k] * eta_dot_[k];
        }

        // Potential energy from constraint
        E_therm += n_dof_ * kT * eta_[0];
        for (size_t k = 1; k < eta_.size(); ++k) {
            E_therm += kT * eta_[k];
        }

        return E_therm;
    }

private:
    Scalar T_target_;
    Scalar tau_;
    int n_dof_ = 0;

    std::vector<Scalar> eta_;      // Thermostat positions
    std::vector<Scalar> eta_dot_;  // Thermostat velocities
    std::vector<Scalar> Q_;        // Thermostat masses
};

/**
 * @brief Initialize velocities from Maxwell-Boltzmann distribution
 *
 * @param system Atomic system (velocities modified)
 * @param temperature Target temperature in Kelvin
 * @param seed Random seed
 * @param remove_com_velocity Remove center of mass velocity
 */
inline void initialize_velocities(AtomicSystem& system, Scalar temperature,
                                  unsigned int seed = 12345,
                                  bool remove_com_velocity = true) {
    std::mt19937 rng(seed);
    std::normal_distribution<Scalar> normal(0.0, 1.0);

    int nat = system.num_atoms();

    // Generate random velocities with Maxwell-Boltzmann distribution
    for (int i = 0; i < nat; ++i) {
        Scalar m = system.mass(i);
        Scalar sigma = std::sqrt(kB_eV_K * temperature / m);

        Vec3 v;
        for (int d = 0; d < 3; ++d) {
            v[d] = sigma * normal(rng);
        }
        system.set_velocity(i, v);
    }

    // Remove center of mass velocity
    if (remove_com_velocity && nat > 0) {
        Vec3 v_com = Vec3::Zero();
        Scalar total_mass = 0.0;

        for (int i = 0; i < nat; ++i) {
            Scalar m = system.mass(i);
            v_com += m * system.velocity(i);
            total_mass += m;
        }
        v_com /= total_mass;

        for (int i = 0; i < nat; ++i) {
            system.set_velocity(i, system.velocity(i) - v_com);
        }
    }

    // Rescale to exact target temperature
    Scalar T_current = 0.0;
    for (int i = 0; i < nat; ++i) {
        Scalar m = system.mass(i);
        T_current += m * system.velocity(i).squaredNorm();
    }
    int n_dof = 3 * nat - (remove_com_velocity ? 3 : 0);
    T_current /= n_dof * kB_eV_K;

    if (T_current > 1e-10) {
        Scalar scale = std::sqrt(temperature / T_current);
        for (int i = 0; i < nat; ++i) {
            system.set_velocity(i, system.velocity(i) * scale);
        }
    }
}

} // namespace atomistica
