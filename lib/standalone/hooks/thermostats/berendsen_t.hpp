// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
// GPL-2.0-or-later — see https://www.gnu.org/licenses/
// ======================================================================

#pragma once

#include <cmath>
#include <stdexcept>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"

namespace atomistica {

class BerendsenT : public Hook {
    double T_;
    double tau_;
    int    d_;   // -1=all, 0=x, 1=y, 2=z

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 20; }
    std::string name() const override { return "BerendsenT"; }

    explicit BerendsenT(const Config& cfg)
        : T_(cfg.get_or<double>("T", 300.0))
        , tau_(cfg.get_or<double>("tau", 100.0))
        , d_(-1)
    {
        std::string ds = cfg.get_or<std::string>("d", "all");
        if      (ds == "x") d_ = 0;
        else if (ds == "y") d_ = 1;
        else if (ds == "z") d_ = 2;
        else if (ds == "all") d_ = -1;
        else throw std::runtime_error("BerendsenT: unknown dimension '" + ds + "'");
    }

    void invoke(SimulationContext& ctx) override {
        double T_curr = temperature_K(ctx.system, /*skip_frozen=*/true);
        if (T_curr < 1e-10) return;

        double scale;
        if (tau_ <= 0.0 || T_curr < 1e-10) {
            scale = std::sqrt(T_ / T_curr);
        } else {
            double ratio = 1.0 + (ctx.dt / tau_) * (T_ / T_curr - 1.0);
            scale = std::sqrt(ratio > 0.0 ? ratio : 0.0);
        }

        size_t n = ctx.system.num_atoms();
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            Vec3 v = ctx.system.velocity(i);
            if (d_ == -1) {
                v *= scale;
            } else {
                v[d_] *= scale;
            }
            ctx.system.set_velocity(i, v);
        }
    }
};

} // namespace atomistica
