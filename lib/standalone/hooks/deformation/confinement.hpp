// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <stdexcept>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"

#include "../../registry.hpp"

namespace atomistica {

class Confinement : public Hook {
    double k_;
    double z_lo_;
    double z_hi_;
    int    d_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_FORCE; }
    int priority() const override { return 12; }
    std::string name() const override { return "Confinement"; }

    explicit Confinement(const Config& cfg)
        : k_(cfg.get_or<double>("k", 10.0))
        , z_lo_(cfg.get_or<double>("z_lo", 0.0))
        , z_hi_(cfg.get_or<double>("z_hi", 0.0))
        , d_(2)
    {
        std::string ds = cfg.get_or<std::string>("d", "z");
        if      (ds == "x") d_ = 0;
        else if (ds == "y") d_ = 1;
        else if (ds == "z") d_ = 2;
        else throw std::runtime_error("Confinement: unknown direction '" + ds + "'");
    }

    void bind_to(AtomicSystem& sys, NeighborList&) override {
        if (z_hi_ == 0.0)
            z_hi_ = sys.cell().col(d_).norm();
    }

    void invoke(SimulationContext& ctx) override {
        size_t n = ctx.system.num_atoms();
        for (size_t i = 0; i < n; ++i) {
            double r_d = ctx.system.position(i)[d_];
            if (r_d < z_lo_)
                ctx.system.forces().col(i)[d_] += k_ * (z_lo_ - r_d);
            else if (r_d > z_hi_)
                ctx.system.forces().col(i)[d_] -= k_ * (r_d - z_hi_);
        }
    }
};

REGISTER_HOOK("Confinement", Confinement)

} // namespace atomistica
