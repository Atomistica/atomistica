// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
// GPL-2.0-or-later — see https://www.gnu.org/licenses/
// ======================================================================

#pragma once

#include <cmath>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"

namespace atomistica {

class RemoveRotation : public Hook {
public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 23; }
    std::string name() const override { return "RemoveRotation"; }

    explicit RemoveRotation(const Config&) {}

    void invoke(SimulationContext& ctx) override {
        size_t n = ctx.system.num_atoms();

        double M = 0.0;
        Vec3   r_cm = Vec3::Zero();
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            M    += m;
            r_cm += m * ctx.system.position(i);
        }
        if (M <= 0.0) return;
        r_cm /= M;

        Vec3 L = Vec3::Zero();
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3   r = ctx.system.position(i) - r_cm;
            L += m * r.cross(ctx.system.velocity(i));
        }
        if (L.norm() < 1e-20) return;

        Mat3 I = Mat3::Zero();
        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            double m = ctx.system.mass(i);
            Vec3   r = ctx.system.position(i) - r_cm;
            I += m * (r.squaredNorm() * Mat3::Identity() - r * r.transpose());
        }

        Vec3 omega = I.inverse() * L;

        for (size_t i = 0; i < n; ++i) {
            if (is_frozen(ctx.system, i)) continue;
            Vec3 r = ctx.system.position(i) - r_cm;
            ctx.system.set_velocity(i,
                ctx.system.velocity(i) - omega.cross(r));
        }
    }
};

} // namespace atomistica
