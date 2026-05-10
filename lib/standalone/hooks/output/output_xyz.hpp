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

class OutputXYZ : public Hook {
    int           freq_;
    std::string   file_;
    bool          write_velocities_;
    bool          write_forces_;
    std::ofstream out_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 51; }
    std::string name() const override { return "OutputXYZ"; }

    explicit OutputXYZ(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 10))
        , file_(cfg.get_or<std::string>("file", "trajectory.xyz"))
        , write_velocities_(cfg.get_or<bool>("velocities", false))
        , write_forces_(cfg.get_or<bool>("forces", false))
    {}

    void bind_to(AtomicSystem&, NeighborList&) override {
        out_.open(file_);
        if (!out_)
            throw std::runtime_error("OutputXYZ: cannot open '" + file_ + "'");
    }

    void invoke(SimulationContext& ctx) override {
        if (ctx.step % freq_ != 0) return;

        size_t n = ctx.system.num_atoms();

        out_ << n << "\n";

        // Properties string
        out_ << "Properties=species:S:1:pos:R:3";
        if (write_velocities_) out_ << ":vel:R:3";
        if (write_forces_)     out_ << ":forces:R:3";
        out_ << " Time=" << std::scientific << std::setprecision(6)
             << ctx.time * ctx.time_display_scale
             << " Step=" << ctx.step << "\n";

        out_ << std::scientific << std::setprecision(8);
        for (size_t i = 0; i < n; ++i) {
            std::string sym = Z_to_symbol(ctx.system.atomic_number(i));
            Vec3 r = ctx.system.position(i);
            out_ << sym
                 << "  " << r[0] << "  " << r[1] << "  " << r[2];
            if (write_velocities_) {
                Vec3 v = ctx.system.velocity(i);
                out_ << "  " << v[0] << "  " << v[1] << "  " << v[2];
            }
            if (write_forces_) {
                Vec3 f = ctx.system.forces().col(i).matrix();
                out_ << "  " << f[0] << "  " << f[1] << "  " << f[2];
            }
            out_ << "\n";
        }
        out_.flush();
    }
};

REGISTER_HOOK("OutputXYZ", OutputXYZ)

} // namespace atomistica
