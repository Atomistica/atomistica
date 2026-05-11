// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cmath>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"

#include "../../registry.hpp"

namespace atomistica {

class ConstantStrainRate : public Hook {
    double gamma_dot_;
    double shear_disp_ = 0.0;
    int    d_shear_;
    int    d_normal_;

public:
    HookPoint hook_point() const override { return HookPoint::PRE_FORCE; }
    int priority() const override { return 0; }
    std::string name() const override { return "ConstantStrainRate"; }

    explicit ConstantStrainRate(const Config& cfg)
        : gamma_dot_(cfg.get_or<double>("gamma_dot", 0.0))
        , d_shear_(cfg.get_or<int>("d_shear", 0))
        , d_normal_(cfg.get_or<int>("d_normal", 1))
    {}

    void invoke(SimulationContext& ctx) override {
        shear_disp_ += gamma_dot_ * ctx.dt;

        Mat3& H = ctx.system.cell();
        double Ln = H.col(d_normal_).norm();
        H(d_shear_, d_normal_) = shear_disp_ * Ln;
        ctx.system.cell_changed();

        // Re-wrap atoms with Lees-Edwards velocity correction using fractional coords.
        Mat3 H_inv = ctx.system.inverse_cell();
        size_t n = ctx.system.num_atoms();

        for (size_t i = 0; i < n; ++i) {
            Vec3 r = ctx.system.position(i);
            // fractional coordinates
            Vec3 s = H_inv * r;

            double sn = s[d_normal_];
            double sn_wrapped = sn - std::floor(sn);
            double dsn = sn_wrapped - sn;

            if (std::abs(dsn) > 0.5) {
                // Wrapping occurred — apply velocity correction
                Vec3 v = ctx.system.velocity(i);
                double correction = gamma_dot_ * Ln;
                if (dsn < 0.0)
                    v[d_shear_] -= correction;
                else
                    v[d_shear_] += correction;
                ctx.system.set_velocity(i, v);
            }

            if (std::abs(dsn) > 1e-12) {
                // Reconstruct wrapped position
                Vec3 s_new = s;
                s_new[d_normal_] = sn_wrapped;
                ctx.system.set_position(i, H * s_new);
            }
        }
    }
};

REGISTER_HOOK("ConstantStrainRate", ConstantStrainRate)

} // namespace atomistica
