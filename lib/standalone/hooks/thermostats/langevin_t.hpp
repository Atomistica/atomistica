// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
// GPL-2.0-or-later — see https://www.gnu.org/licenses/
// ======================================================================

#pragma once

#include <cmath>
#include <random>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"

namespace atomistica {

class LangevinT : public Hook {
    double       T_;
    double       gamma_;
    std::mt19937 rng_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 21; }
    std::string name() const override { return "LangevinT"; }

    explicit LangevinT(const Config& cfg)
        : T_(cfg.get_or<double>("T", 300.0))
        , gamma_(cfg.get_or<double>("gamma", 0.01))
        , rng_(static_cast<unsigned>(cfg.get_or<int>("seed", 12345)))
    {}

    void invoke(SimulationContext& ctx) override {
        std::normal_distribution<double> normal(0.0, 1.0);
        double dt = ctx.dt;
        size_t n  = ctx.system.num_atoms();

        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m_i   = ctx.system.mass(i);
            double sigma = std::sqrt(2.0 * gamma_ * kB_eV * T_ / m_i * dt);
            Vec3 v = ctx.system.velocity(i);
            v[0] += -gamma_ * dt * v[0] + sigma * normal(rng_);
            v[1] += -gamma_ * dt * v[1] + sigma * normal(rng_);
            v[2] += -gamma_ * dt * v[2] + sigma * normal(rng_);
            ctx.system.set_velocity(i, v);
        }
    }
};

} // namespace atomistica
