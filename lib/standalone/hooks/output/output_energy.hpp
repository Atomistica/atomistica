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
#include "../../md_utils.hpp"
#include "../../registry.hpp"

namespace atomistica {

class OutputEnergy : public Hook {
    int           freq_;
    std::string   file_;
    std::ofstream out_;

public:
    HookPoint hook_point() const override { return HookPoint::POST_STEP; }
    int priority() const override { return 50; }
    std::string name() const override { return "OutputEnergy"; }

    explicit OutputEnergy(const Config& cfg)
        : freq_(cfg.get_or<int>("freq", 10))
        , file_(cfg.get_or<std::string>("file", "thermo.dat"))
    {}

    void bind_to(AtomicSystem&, NeighborList&) override {
        out_.open(file_);
        if (!out_)
            throw std::runtime_error("OutputEnergy: cannot open '" + file_ + "'");
        out_ << "# step  time[" << "unit]  Ekin[eV]  Epot[eV]  Etot[eV]"
                "  fmax[eV/A]  T[K]  P[eV/A3]\n";
        out_.flush();
    }

    void invoke(SimulationContext& ctx) override {
        if (ctx.step % freq_ != 0) return;

        double ekin  = kinetic_energy(ctx.system);
        double epot  = ctx.results.energy;
        double etot  = ekin + epot;
        double T_K   = temperature_K(ctx.system);
        double P     = pressure(ctx.results.virial, ctx.system)
                       * ctx.pressure_display_scale;
        double fm    = fmax(ctx.system);
        double t_disp = ctx.time * ctx.time_display_scale;

        out_ << ctx.step
             << "  " << std::scientific << std::setprecision(6) << t_disp
             << "  " << ekin
             << "  " << epot
             << "  " << etot
             << "  " << fm
             << "  " << std::fixed << std::setprecision(3) << T_K
             << "  " << std::scientific << std::setprecision(6) << P
             << "\n";
        out_.flush();
    }
};

REGISTER_HOOK("OutputEnergy", OutputEnergy)

} // namespace atomistica
