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
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"

namespace atomistica {

// Peters DPD thermostat as a Hook.
// Reference: E.A.J.F. Peters, Europhys. Lett. 66, 311 (2004).
class PetersT : public Hook {
    double       T_;
    double       gamma_;
    double       cutoff_;
    std::mt19937 rng_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 22; }
    std::string name() const override { return "PetersT"; }

    explicit PetersT(const Config& cfg)
        : T_(cfg.get_or<double>("T", 300.0))
        , gamma_(cfg.get_or<double>("gamma", 1.0))
        , cutoff_(cfg.get_or<double>("cutoff", 0.0))
        , rng_(static_cast<unsigned>(cfg.get_or<int>("seed", 12345)))
    {}

    void bind_to(AtomicSystem& /*system*/, NeighborList& nl) override {
        if (cutoff_ <= 0.0)
            cutoff_ = nl.cutoff();
    }

    void invoke(SimulationContext& ctx) override {
        std::normal_distribution<double> normal(0.0, 1.0);
        double dt        = ctx.dt;
        double cutoff2   = cutoff_ * cutoff_;
        size_t nat       = ctx.system.num_atoms();

        for (size_t i = 0; i < nat; ++i) {
            auto [begin, end] = ctx.nl.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                size_t j = it->index;
                if (j <= i) continue;

                Vec3 shift;
                shift << static_cast<double>(it->cell_shift[0]),
                         static_cast<double>(it->cell_shift[1]),
                         static_cast<double>(it->cell_shift[2]);
                Vec3 dr = ctx.system.position(j)
                          + ctx.system.cell() * shift
                          - ctx.system.position(i);

                double r2 = dr.squaredNorm();
                if (r2 >= cutoff2) continue;

                double r    = std::sqrt(r2);
                Vec3   rhat = dr / r;

                double m_i    = ctx.system.mass(i);
                double m_j    = ctx.system.mass(j);
                double r_muij = 1.0 / m_i + 1.0 / m_j;
                double mu     = 1.0 / r_muij;

                double kernel = 1.0 - r / cutoff_;
                double w      = dt * r_muij * gamma_ * kernel;

                double a     = mu * (1.0 - std::exp(-w));
                double b_var = kB_eV * T_ * mu * (1.0 - std::exp(-2.0 * w));
                double b     = std::sqrt(b_var > 0.0 ? b_var : 0.0) * normal(rng_);

                Vec3 vij  = ctx.system.velocity(i) - ctx.system.velocity(j);
                Vec3 dmom = (-a * vij.dot(rhat) + b) * rhat;

                ctx.system.set_velocity(i, ctx.system.velocity(i) + dmom / m_i);
                ctx.system.set_velocity(j, ctx.system.velocity(j) - dmom / m_j);
            }
        }
    }
};

} // namespace atomistica
