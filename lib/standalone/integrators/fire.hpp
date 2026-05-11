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
#include <string>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../simulation_context.hpp"
#include "../hook.hpp"
#include "../integrator.hpp"
#include "../md_utils.hpp"

namespace atomistica {

// Fast Inertial Relaxation Engine (FIRE) structure relaxation.
// Bitzek et al., PRL 97, 170201 (2006).
class FIRE : public Integrator {
public:
    explicit FIRE(const Config& cfg)
        : dt_(0.0),
          dt_max_(0.0),
          dt_init_(0.0),
          alpha_(cfg.get_or<double>("alpha_start", 0.1)),
          alpha_start_(cfg.get_or<double>("alpha_start", 0.1)),
          f_inc_(cfg.get_or<double>("f_inc", 1.1)),
          f_dec_(cfg.get_or<double>("f_dec", 0.5)),
          f_alpha_(cfg.get_or<double>("f_alpha", 0.99)),
          n_min_(cfg.get_or<int>("n_min", 5)),
          n_positive_(0),
          fmax_tol_(cfg.get_or<double>("fmax_tol", 0.01)),
          initialized_(false) {}

    std::string name() const override { return "FIRE"; }

    bool step(SimulationContext& ctx) override {
        if (!initialized_) {
            dt_      = ctx.dt;
            dt_max_  = 10.0 * dt_;
            dt_init_ = dt_;
            initialized_ = true;
        }

        const size_t n = ctx.system.num_atoms();

        // Half-kick
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3 f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            Vec3 v = ctx.system.velocity(i);
            ctx.system.set_velocity(i, v + 0.5 * dt_ * f / m);
        }

        // Drift
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            Vec3 r = ctx.system.position(i);
            Vec3 v = ctx.system.velocity(i);
            ctx.system.set_position(i, r + dt_ * v);
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
            Vec3 f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            Vec3 v = ctx.system.velocity(i);
            ctx.system.set_velocity(i, v + 0.5 * dt_ * f / m);
        }

        // Compute P = F·v, |v|, |F|
        double P = 0.0;
        double v_norm2 = 0.0;
        double f_norm2 = 0.0;
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            Vec3 v = ctx.system.velocity(i);
            Vec3 f = ctx.system.forces().col(static_cast<int>(i)).matrix();
            P       += f.dot(v);
            v_norm2 += v.squaredNorm();
            f_norm2 += f.squaredNorm();
        }
        double v_norm = std::sqrt(v_norm2);
        double f_norm = std::sqrt(f_norm2);

        // FIRE velocity mixing: v_i = (1-alpha)*v_i + alpha*(|v|/|F|)*F_i
        if (f_norm > 0.0) {
            double scale = alpha_ * v_norm / f_norm;
            for (size_t i = 0; i < n; ++i) {
                if (is_frozen(ctx.system, i)) continue;
                Vec3 v = ctx.system.velocity(i);
                Vec3 f = ctx.system.forces().col(static_cast<int>(i)).matrix();
                ctx.system.set_velocity(i, (1.0 - alpha_) * v + scale * f);
            }
        }

        // FIRE parameter update
        if (P > 0.0) {
            ++n_positive_;
            if (n_positive_ >= n_min_) {
                dt_    = std::min(dt_ * f_inc_, dt_max_);
                alpha_ *= f_alpha_;
            }
        } else {
            dt_         *= f_dec_;
            alpha_       = alpha_start_;
            n_positive_  = 0;
            for (size_t i = 0; i < n; ++i) {
                if (is_frozen(ctx.system, i)) continue;
                ctx.system.set_velocity(i, Vec3::Zero());
            }
        }

        ctx.invoke_hooks(HookPoint::POST_STEP);

        ctx.dt    = dt_;
        ctx.time += dt_;
        ++ctx.step;

        if (fmax(ctx.system) < fmax_tol_)
            return false;

        return true;
    }

private:
    double dt_;
    double dt_max_;
    double dt_init_;
    double alpha_;
    double alpha_start_;
    double f_inc_;
    double f_dec_;
    double f_alpha_;
    int    n_min_;
    int    n_positive_;
    double fmax_tol_;
    bool   initialized_;
};

} // namespace atomistica
