// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cmath>
#include <fstream>
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

class Slicing : public Hook {
    int         freq_;
    std::string file_;
    int         n_bins_;
    int         d_;
    std::vector<double> density_;
    std::vector<double> vx_, vy_, vz_;
    std::vector<double> T_bin_;
    std::vector<int>    count_;
    int         acc_count_ = 0;
    std::ofstream out_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 41; }
    std::string name() const override { return "Slicing"; }

    explicit Slicing(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 100))
        , file_(cfg.get_or<std::string>("file", "slicing.dat"))
        , n_bins_(cfg.get_or<int>("n_bins", 100))
        , d_(2)
    {
        std::string ds = cfg.get_or<std::string>("d", "z");
        if      (ds == "x") d_ = 0;
        else if (ds == "y") d_ = 1;
        else if (ds == "z") d_ = 2;
        else throw std::runtime_error("Slicing: unknown direction '" + ds + "'");
    }

    void bind_to(AtomicSystem&, NeighborList&) override {
        density_.assign(n_bins_, 0.0);
        vx_.assign(n_bins_, 0.0);
        vy_.assign(n_bins_, 0.0);
        vz_.assign(n_bins_, 0.0);
        T_bin_.assign(n_bins_, 0.0);
        count_.assign(n_bins_, 0);
        acc_count_ = 0;

        out_.open(file_);
        if (!out_)
            throw std::runtime_error("Slicing: cannot open '" + file_ + "'");
        out_ << "# bin center density vx vy vz T\n";
        out_.flush();
    }

    void invoke(SimulationContext& ctx) override {
        double cell_length = ctx.system.cell().col(d_).norm();
        size_t n = ctx.system.num_atoms();

        for (size_t i = 0; i < n; ++i) {
            double r_d = ctx.system.position(i)[d_];
            int b = static_cast<int>((r_d / cell_length) * n_bins_) % n_bins_;
            if (b < 0) b += n_bins_;
            if (b >= n_bins_) b = n_bins_ - 1;

            Vec3 v = ctx.system.velocity(i);
            double m = ctx.system.mass(i);

            density_[b] += 1.0;
            vx_[b]      += v[0];
            vy_[b]      += v[1];
            vz_[b]      += v[2];
            T_bin_[b]   += m * v.squaredNorm();
            count_[b]   += 1;
        }
        ++acc_count_;

        if (acc_count_ < freq_) return;

        double bin_vol = cell_length / n_bins_;
        double cell_area = ctx.system.volume() / cell_length;
        double vol_bin = bin_vol * cell_area;

        for (int b = 0; b < n_bins_; ++b) {
            int c = count_[b];
            double center = (b + 0.5) / n_bins_;
            double dens = density_[b] / (acc_count_ * vol_bin);
            double avg_vx = c > 0 ? vx_[b] / c : 0.0;
            double avg_vy = c > 0 ? vy_[b] / c : 0.0;
            double avg_vz = c > 0 ? vz_[b] / c : 0.0;
            // T = (1/3N) sum m*v^2 / kB  — but here we use (2/3) ekin / (N kB)
            double T = 0.0;
            if (c > 0) {
                double ekin = 0.5 * T_bin_[b] / c;
                T = 2.0 * ekin / (3.0 * kB_eV);
            }
            out_ << b << " " << center << " " << dens << " "
                 << avg_vx << " " << avg_vy << " " << avg_vz << " " << T << "\n";
        }
        out_ << "\n";
        out_.flush();

        density_.assign(n_bins_, 0.0);
        vx_.assign(n_bins_, 0.0);
        vy_.assign(n_bins_, 0.0);
        vz_.assign(n_bins_, 0.0);
        T_bin_.assign(n_bins_, 0.0);
        count_.assign(n_bins_, 0);
        acc_count_ = 0;
    }
};

REGISTER_HOOK("Slicing", Slicing)

} // namespace atomistica
