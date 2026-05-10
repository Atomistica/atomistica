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

class ConstantForce : public Hook {
    int  group_id_;
    Vec3 force_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_FORCE; }
    int priority() const override { return 10; }
    std::string name() const override { return "ConstantForce"; }

    explicit ConstantForce(const Config& cfg)
        : group_id_(cfg.get_or<int>("group", 0))
    {
        force_[0] = cfg.get_or<double>("fx", 0.0);
        force_[1] = cfg.get_or<double>("fy", 0.0);
        force_[2] = cfg.get_or<double>("fz", 0.0);
    }

    void invoke(SimulationContext& ctx) override {
        size_t n = ctx.system.num_atoms();
        if (group_id_ == 0) {
            for (size_t i = 0; i < n; ++i)
                ctx.system.forces().col(i) += force_.array();
        } else {
            if (!ctx.system.properties().has("group")) return;
            const ArrayXi& groups = ctx.system.properties().get<ArrayXi>("group");
            for (size_t i = 0; i < n; ++i) {
                if (groups[static_cast<int>(i)] == group_id_)
                    ctx.system.forces().col(i) += force_.array();
            }
        }
    }
};

REGISTER_HOOK("ConstantForce", ConstantForce)

} // namespace atomistica
