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

#include <cmath>
#include <random>
#include <string>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../simulation_context.hpp"
#include "../hook.hpp"
#include "../integrator.hpp"
#include "../md_utils.hpp"

namespace atomistica {

// NPT integrator combining velocity-Verlet with Berendsen pressure coupling
// and an optional Andersen thermostat.
//
// Pressure coupling: Berendsen rescaling with coupling constant tau_P.
//   mu = cbrt(1 - kappa * dt/tau_P * (P_target - P_current))
//   cell → mu * cell,  r_i → mu * r_i
//
// Temperature control (optional, Andersen stochastic thermostat):
//   Each atom is independently reassigned a Maxwell-Boltzmann velocity with
//   probability nu * dt per step.
//
// Config keys:
//   P      — target pressure in eV/Å³ (default 0.0)
//   tau    — pressure coupling time in internal time units (default 1000*dt)
//   W      — alias for tau (for compatibility with Fortran md.dat files)
//   kappa  — dimensionless compressibility prefactor (default 1.0)
//   d      — coupling direction: "all" / "x" / "y" / "z" (default "all")
//   T      — Andersen thermostat target temperature in K (default 0 = disabled)
//   nu     — Andersen collision frequency in 1/time_unit (default 0.01)
//   seed   — random seed (default 12345)
class AndersenP : public Integrator {
    double P_     = 0.0;
    double tau_   = 0.0;    // coupling time; 0 = set from dt on first step
    double kappa_ = 1.0;
    int    d_     = -1;     // -1=all, 0=x, 1=y, 2=z
    double T_     = 0.0;    // Andersen thermostat temperature (0 = disabled)
    double nu_    = 0.01;
    std::mt19937 rng_;
    bool   tau_from_config_ = false;

public:
    explicit AndersenP(const Config& cfg)
        : P_    (cfg.get_or<double>("P",     0.0))
        , kappa_(cfg.get_or<double>("kappa", 1.0))
        , T_    (cfg.get_or<double>("T",     0.0))
        , nu_   (cfg.get_or<double>("nu",    0.01))
        , rng_  (static_cast<unsigned>(cfg.get_or<int>("seed", 12345)))
    {
        // Accept either "tau" or "W" for the coupling parameter.
        if (cfg.has_key("tau")) {
            tau_            = cfg.get_or<double>("tau", 0.0);
            tau_from_config_ = true;
        } else if (cfg.has_key("W")) {
            tau_            = cfg.get_or<double>("W", 0.0);
            tau_from_config_ = true;
        }

        std::string ds = cfg.get_or<std::string>("d", "all");
        if      (ds == "x")   d_ = 0;
        else if (ds == "y")   d_ = 1;
        else if (ds == "z")   d_ = 2;
        else if (ds == "all") d_ = -1;
    }

    std::string name() const override { return "AndersenP"; }

    bool step(SimulationContext& ctx) override {
        // Default tau = 1000 * dt if not set by config.
        if (!tau_from_config_ && tau_ <= 0.0)
            tau_ = 1000.0 * ctx.dt;

        const double dt = ctx.dt;
        const size_t n  = ctx.system.num_atoms();

        // ---- Velocity-Verlet ------------------------------------------------

        // Half-kick (v += 0.5*dt*f/m)
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3   f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            ctx.system.set_velocity(i, ctx.system.velocity(i) + 0.5*dt*f/m);
        }

        // Drift (r += dt*v)
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            ctx.system.set_position(i,
                ctx.system.position(i) + dt * ctx.system.velocity(i));
        }
        wrap_positions(ctx.system);

        ctx.invoke_hooks(HookPoint::PRE_FORCE);
        ctx.nl.update(ctx.system);
        ctx.compute_forces();
        ctx.invoke_hooks(HookPoint::POST_FORCE);

        // Second half-kick
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3   f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            ctx.system.set_velocity(i, ctx.system.velocity(i) + 0.5*dt*f/m);
        }

        // ---- Berendsen pressure coupling ------------------------------------
        double P_curr = pressure(ctx.results.virial, ctx.system);
        double mu     = std::cbrt(1.0 - kappa_ * (dt / tau_) * (P_ - P_curr));

        if (d_ == -1) {
            // Isotropic scaling
            Mat3 new_cell = ctx.system.cell() * static_cast<Scalar>(mu);
            ctx.system.set_cell(new_cell);
            for (size_t i = 0; i < n; ++i)
                ctx.system.set_position(i, ctx.system.position(i) * static_cast<Scalar>(mu));
        } else {
            // Uniaxial scaling along direction d_
            Mat3 new_cell = ctx.system.cell();
            new_cell.col(d_) *= static_cast<Scalar>(mu);
            ctx.system.set_cell(new_cell);
            for (size_t i = 0; i < n; ++i) {
                Vec3 r = ctx.system.position(i);
                r[d_] *= static_cast<Scalar>(mu);
                ctx.system.set_position(i, r);
            }
        }
        wrap_positions(ctx.system);
        ctx.system.cell_changed();
        ctx.system.positions_changed();

        // ---- Andersen stochastic thermostat ---------------------------------
        if (T_ > 0.0 && nu_ > 0.0) {
            std::normal_distribution<double>  normal(0.0, 1.0);
            std::uniform_real_distribution<double> uni(0.0, 1.0);
            for (size_t i = 0; i < n; ++i) {
                if (is_frozen(ctx.system, i)) continue;
                if (uni(rng_) < nu_ * dt) {
                    double sigma = std::sqrt(kB_eV * T_ / ctx.system.mass(i));
                    Vec3 v(sigma*normal(rng_), sigma*normal(rng_), sigma*normal(rng_));
                    ctx.system.set_velocity(i, v);
                }
            }
        }

        ctx.invoke_hooks(HookPoint::POST_STEP);
        ctx.time += ctx.dt;
        ++ctx.step;
        return true;
    }
};

} // namespace atomistica
