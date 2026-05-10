// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"

#include "../../registry.hpp"

namespace atomistica {

class HeatFlux : public Hook {
    int         freq_;
    int         n_slabs_;
    std::string file_;
    double      exchanged_energy_ = 0.0;
    int         start_step_ = -1;
    std::ofstream out_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 42; }
    std::string name() const override { return "HeatFlux"; }

    explicit HeatFlux(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 50))
        , n_slabs_(cfg.get_or<int>("n_slabs", 20))
        , file_(cfg.get_or<std::string>("file", "heatflux.dat"))
    {}

    void bind_to(AtomicSystem&, NeighborList&) override {
        out_.open(file_);
        if (!out_)
            throw std::runtime_error("HeatFlux: cannot open '" + file_ + "'");
        out_ << "# step time cumulative_exchanged_energy\n";
        out_.flush();
    }

    void invoke(SimulationContext& ctx) override {
        if (ctx.step % freq_ != 0) return;

        if (start_step_ < 0)
            start_step_ = ctx.step;

        double Lz = ctx.system.cell().col(2).norm();
        size_t n  = ctx.system.num_atoms();
        int hot_slab  = 0;
        int cold_slab = n_slabs_ / 2;

        size_t hot_idx  = static_cast<size_t>(-1);
        size_t cold_idx = static_cast<size_t>(-1);
        double min_ekin = std::numeric_limits<double>::max();
        double max_ekin = -std::numeric_limits<double>::max();

        for (size_t i = 0; i < n; ++i) {
            double rz = ctx.system.position(i)[2];
            int slab = static_cast<int>(rz / Lz * n_slabs_) % n_slabs_;
            if (slab < 0) slab += n_slabs_;

            double m  = ctx.system.mass(i);
            double vz = ctx.system.velocity(i)[2];
            double ek = 0.5 * m * vz * vz;

            if (slab == hot_slab) {
                if (ek > max_ekin) {
                    max_ekin = ek;
                    hot_idx  = i;
                }
            } else if (slab == cold_slab) {
                if (ek < min_ekin) {
                    min_ekin = ek;
                    cold_idx = i;
                }
            }
        }

        if (hot_idx == static_cast<size_t>(-1) || cold_idx == static_cast<size_t>(-1))
            return;

        Vec3 v_hot  = ctx.system.velocity(hot_idx);
        Vec3 v_cold = ctx.system.velocity(cold_idx);
        double m_hot  = ctx.system.mass(hot_idx);
        double m_cold = ctx.system.mass(cold_idx);

        double ek_hot_old  = 0.5 * m_hot  * v_hot[2]  * v_hot[2];
        double ek_cold_old = 0.5 * m_cold * v_cold[2] * v_cold[2];

        // Swap vz components
        double tmp  = v_hot[2];
        v_hot[2]    = v_cold[2];
        v_cold[2]   = tmp;

        double ek_hot_new = 0.5 * m_hot * v_hot[2] * v_hot[2];

        // Only swap if it transfers energy from hot to cold
        // (hot loses energy, cold gains energy)
        double delta_hot = ek_hot_new - ek_hot_old;
        if (delta_hot < 0.0) {
            ctx.system.set_velocity(hot_idx,  v_hot);
            ctx.system.set_velocity(cold_idx, v_cold);
            exchanged_energy_ += -delta_hot;
        }
        // else: leave velocities unchanged

        double t = static_cast<double>(ctx.step - start_step_) * ctx.dt;
        out_ << ctx.step << " " << t << " " << exchanged_energy_ << "\n";
        out_.flush();
    }
};

REGISTER_HOOK("HeatFlux", HeatFlux)

} // namespace atomistica
