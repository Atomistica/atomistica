// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <string>
#include <vector>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"

#include "../../registry.hpp"

namespace atomistica {

class HarmonicHook : public Hook {
    double k_;
    int    group_id_;
    std::vector<Vec3> r_ref_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_FORCE; }
    int priority() const override { return 11; }
    std::string name() const override { return "HarmonicHook"; }

    explicit HarmonicHook(const Config& cfg)
        : k_(cfg.get_or<double>("k", 1.0))
        , group_id_(cfg.get_or<int>("group", 0))
    {}

    void bind_to(AtomicSystem& sys, NeighborList&) override {
        size_t n = sys.num_atoms();
        r_ref_.resize(n);
        for (size_t i = 0; i < n; ++i)
            r_ref_[i] = sys.position(i);
    }

    void invoke(SimulationContext& ctx) override {
        size_t n = ctx.system.num_atoms();
        if (group_id_ == 0) {
            for (size_t i = 0; i < n; ++i) {
                Vec3 dr = ctx.system.position(i) - r_ref_[i];
                ctx.system.forces().col(i) += (-k_ * dr).array();
            }
        } else {
            if (!ctx.system.properties().has("group")) return;
            const ArrayXi& groups = ctx.system.properties().get<ArrayXi>("group");
            for (size_t i = 0; i < n; ++i) {
                if (groups[static_cast<int>(i)] == group_id_) {
                    Vec3 dr = ctx.system.position(i) - r_ref_[i];
                    ctx.system.forces().col(i) += (-k_ * dr).array();
                }
            }
        }
    }
};

REGISTER_HOOK("HarmonicHook", HarmonicHook)

} // namespace atomistica
