// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../../include/atomistica/core/neighbor_list.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../md_utils.hpp"
#include "../../atoms_io.hpp"

#include "../../registry.hpp"

namespace atomistica {

class DiffusionCoefficient : public Hook {
    int         freq_;
    std::string file_;
    std::ofstream out_;
    std::vector<Vec3> r0_;
    int         start_step_ = -1;
    bool        per_element_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 40; }
    std::string name() const override { return "DiffusionCoefficient"; }

    explicit DiffusionCoefficient(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 10))
        , file_(cfg.get_or<std::string>("file", "diffusion.dat"))
        , per_element_(cfg.get_or<bool>("per_element", false))
    {}

    void bind_to(AtomicSystem& sys, NeighborList&) override {
        out_.open(file_);
        if (!out_)
            throw std::runtime_error("DiffusionCoefficient: cannot open '" + file_ + "'");
        out_ << "# step time MSD D";
        if (per_element_) {
            std::map<int,bool> seen;
            for (size_t i = 0; i < sys.num_atoms(); ++i) {
                int Z = sys.atomic_number(i);
                if (!seen[Z]) {
                    seen[Z] = true;
                    std::string sym = Z_to_symbol(Z);
                    out_ << " MSD_" << sym << " D_" << sym;
                }
            }
        }
        out_ << "\n";
        out_.flush();
    }

    void invoke(SimulationContext& ctx) override {
        if (ctx.step % freq_ != 0) return;

        size_t n = ctx.system.num_atoms();

        if (start_step_ < 0) {
            r0_.resize(n);
            for (size_t i = 0; i < n; ++i)
                r0_[i] = ctx.system.position(i);
            start_step_ = ctx.step;
        }

        double t = static_cast<double>(ctx.step - start_step_) * ctx.dt;

        double msd = 0.0;
        for (size_t i = 0; i < n; ++i) {
            Vec3 dr = ctx.system.position(i) - r0_[i];
            msd += dr.squaredNorm();
        }
        msd /= static_cast<double>(n);

        out_ << ctx.step << " " << t;
        out_ << " " << msd;

        if (t < 1e-20) {
            out_ << " 0.0";
        } else {
            out_ << " " << msd / (6.0 * t);
        }

        if (per_element_) {
            std::map<int, std::pair<double,int>> elem_msd;
            for (size_t i = 0; i < n; ++i) {
                int Z = ctx.system.atomic_number(i);
                Vec3 dr = ctx.system.position(i) - r0_[i];
                elem_msd[Z].first  += dr.squaredNorm();
                elem_msd[Z].second += 1;
            }
            for (auto& kv : elem_msd) {
                double emsd = kv.second.second > 0
                    ? kv.second.first / static_cast<double>(kv.second.second)
                    : 0.0;
                double eD = (t < 1e-20) ? 0.0 : emsd / (6.0 * t);
                out_ << " " << emsd << " " << eD;
            }
        }

        out_ << "\n";
        out_.flush();
    }
};

REGISTER_HOOK("DiffusionCoefficient", DiffusionCoefficient)

} // namespace atomistica
