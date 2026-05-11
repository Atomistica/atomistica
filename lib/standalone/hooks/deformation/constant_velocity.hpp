// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"

#include "../../registry.hpp"

namespace atomistica {

class ConstantVelocity : public Hook {
    int  group_id_;
    Vec3 vel_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 25; }
    std::string name() const override { return "ConstantVelocity"; }

    explicit ConstantVelocity(const Config& cfg)
        : group_id_(cfg.get_or<int>("group", 1))
    {
        vel_[0] = cfg.get_or<double>("vx", 0.0);
        vel_[1] = cfg.get_or<double>("vy", 0.0);
        vel_[2] = cfg.get_or<double>("vz", 0.0);
    }

    void invoke(SimulationContext& ctx) override {
        if (!ctx.system.properties().has("group")) return;
        const ArrayXi& groups = ctx.system.properties().get<ArrayXi>("group");
        size_t n = ctx.system.num_atoms();
        for (size_t i = 0; i < n; ++i) {
            if (groups[static_cast<int>(i)] == group_id_)
                ctx.system.set_velocity(i, vel_);
        }
    }
};

REGISTER_HOOK("ConstantVelocity", ConstantVelocity)

} // namespace atomistica
