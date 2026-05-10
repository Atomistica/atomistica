// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

#include "../../config.hpp"
#include "../../../include/atomistica/core/atomic_system.hpp"
#include "../../simulation_context.hpp"
#include "../../hook.hpp"
#include "../../atoms_io.hpp"
#include "../../registry.hpp"

namespace atomistica {

class OutputCFG : public Hook {
    int           freq_;
    std::string   file_;
    bool          write_velocities_;
    std::ofstream out_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 52; }
    std::string name() const override { return "OutputCFG"; }

    explicit OutputCFG(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 10))
        , file_(cfg.get_or<std::string>("file", "trajectory.cfg"))
        , write_velocities_(cfg.get_or<bool>("velocities", false))
    {}

    void bind_to(AtomicSystem&, NeighborList&) override {
        out_.open(file_);
        if (!out_)
            throw std::runtime_error("OutputCFG: cannot open '" + file_ + "'");
    }

    void invoke(SimulationContext& ctx) override {
        if (ctx.step % freq_ != 0) return;

        size_t n = ctx.system.num_atoms();
        const Mat3& H = ctx.system.cell();
        Mat3 H_inv = H.inverse();

        out_ << std::scientific << std::setprecision(10);

        out_ << "Number of particles = " << n << "\n";
        out_ << "A = 1.0 Angstrom (basic length-scale)\n";

        // H0(i,j): column j is lattice vector a_j, so H0(i,j) = cell(i-1,j-1)
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                out_ << "H0(" << (i+1) << "," << (j+1) << ") = "
                     << H(i, j) << "\n";

        if (!write_velocities_)
            out_ << ".NO_VELOCITY.\n";

        int entry_count = write_velocities_ ? 6 : 3;
        out_ << "entry_count = " << entry_count << "\n";

        for (size_t i = 0; i < n; ++i) {
            int Z = ctx.system.atomic_number(i);
            double mass_amu = standard_atomic_mass(Z);
            std::string sym = Z_to_symbol(Z);
            Vec3 s = H_inv * ctx.system.position(i);

            out_ << mass_amu << "\n";
            out_ << sym << "\n";
            out_ << s[0] << " " << s[1] << " " << s[2];
            if (write_velocities_) {
                Vec3 v = ctx.system.velocity(i);
                out_ << " " << v[0] << " " << v[1] << " " << v[2];
            }
            out_ << "\n";
        }
        out_.flush();
    }
};

REGISTER_HOOK("OutputCFG", OutputCFG)

} // namespace atomistica
