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

class BerendsenP : public Hook {
    double P_;
    double tau_;
    double kappa_;
    int    d_;   // -1=isotropic, 0=x, 1=y, 2=z

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 30; }
    std::string name() const override { return "BerendsenP"; }

    explicit BerendsenP(const Config& cfg)
        : P_(cfg.get_or<double>("P", 0.0))
        , tau_(cfg.get_or<double>("tau", 1000.0))
        , kappa_(cfg.get_or<double>("kappa", 1.0))
        , d_(-1)
    {
        std::string ds = cfg.get_or<std::string>("d", "all");
        if      (ds == "x")   d_ = 0;
        else if (ds == "y")   d_ = 1;
        else if (ds == "z")   d_ = 2;
        else if (ds == "all") d_ = -1;
        else throw std::runtime_error("BerendsenP: unknown dimension '" + ds + "'");
    }

    void invoke(SimulationContext& ctx) override {
        double P_curr = pressure(ctx.results.virial, ctx.system);
        double mu = std::cbrt(1.0 - kappa_ * (ctx.dt / tau_) * (P_ - P_curr));

        size_t n = ctx.system.num_atoms();

        if (d_ == -1) {
            Mat3 new_cell = ctx.system.cell() * mu;
            ctx.system.set_cell(new_cell);
            for (size_t i = 0; i < n; ++i) {
                Vec3 r = ctx.system.position(i) * mu;
                ctx.system.set_position(i, r);
            }
        } else {
            Mat3 new_cell = ctx.system.cell();
            new_cell.col(d_) *= mu;
            ctx.system.set_cell(new_cell);
            for (size_t i = 0; i < n; ++i) {
                Vec3 r = ctx.system.position(i);
                r[d_] *= mu;
                ctx.system.set_position(i, r);
            }
        }

        // Wrap positions into the new cell
        wrap_positions(ctx.system);

        ctx.system.cell_changed();
        ctx.system.positions_changed();
    }
};

} // namespace atomistica
