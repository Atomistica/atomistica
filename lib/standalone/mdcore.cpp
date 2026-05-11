// ======================================================================
// Atomistica - Interatomic potential library and molecular dynamics code
// https://github.com/Atomistica/atomistica
//
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// and others. See the AUTHORS file in the top-level Atomistica directory.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
// ======================================================================

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include "../include/atomistica/config.hpp"
#include "../include/atomistica/core/atomic_system.hpp"
#include "../include/atomistica/core/neighbor_list.hpp"
#include "../include/atomistica/potentials/potential_base.hpp"

#include "config.hpp"
#include "atoms_io.hpp"
#include "md_utils.hpp"
#include "simulation_context.hpp"
#include "hook.hpp"
#include "integrator.hpp"
#include "registry.hpp"

using namespace atomistica;

int main() {
    // -----------------------------------------------------------------------
    // 1. Parse configuration
    // -----------------------------------------------------------------------
    Config cfg;
    cfg.parse_file("md.dat");

    std::string sou = cfg.get_string("system_of_units", "eV/A");
    bool fs_mode = (sou.find("fs") != std::string::npos);

    double dt_raw     = cfg.get_double("dt", 0.1);
    double max_time_raw = cfg.get_double("max_time", -1.0);
    long   max_iter   = static_cast<long>(cfg.get_int("n_iterations", -1));
    int    scr_freq   = cfg.get_int("scr_freq", 10);
    int    file_freq  = cfg.get_int("file_freq", 10);
    double cutoff_add = cfg.get_double("cutoff_add", 0.5);

    // If neither max_time nor n_iterations specified, default to 100 time units
    bool has_max_time = (max_time_raw > 0.0);
    bool has_max_iter = (max_iter > 0);
    if (!has_max_time && !has_max_iter) {
        max_time_raw = 100.0;
        has_max_time = true;
    }

    // -----------------------------------------------------------------------
    // 2. Read atoms
    // -----------------------------------------------------------------------
    AtomsData atoms_data = read_atoms_dat("atoms.dat");
    AtomicSystem& sys = atoms_data.system;

    // -----------------------------------------------------------------------
    // 3. Setup unit system
    // -----------------------------------------------------------------------
    const double sqrt_c = SQRT_AMU_AFSQ_PER_EV;  // ≈ 10.18

    double dt = fs_mode ? dt_raw / sqrt_c : dt_raw;
    double max_time_internal = fs_mode ? max_time_raw / sqrt_c : max_time_raw;

    const char* time_label     = fs_mode ? "fs" : "10fs";
    const char* pressure_label = "eV/A^3";

    // -----------------------------------------------------------------------
    // 4. Create SimulationContext
    // -----------------------------------------------------------------------
    NeighborList nl;
    SimulationContext ctx(sys, nl);

    // -----------------------------------------------------------------------
    // 5. Configure context timing and display units
    // -----------------------------------------------------------------------
    ctx.dt = dt;
    ctx.time_display_scale = fs_mode ? sqrt_c : 1.0;
    ctx.time_unit_label    = time_label;
    ctx.pressure_display_scale = 1.0;
    ctx.pressure_unit_label    = pressure_label;

    // -----------------------------------------------------------------------
    // 6. Dispatch all sections to Registry
    // -----------------------------------------------------------------------
    Registry& reg = Registry::instance();
    std::unique_ptr<Integrator> integrator;

    for (const auto& [name, section_cfg] : cfg.all_sections()) {
        if (reg.is_integrator(name)) {
            if (integrator)
                throw std::runtime_error(
                    "Multiple integrator sections in md.dat — only one allowed");
            integrator = reg.make_integrator(name, section_cfg);
        } else if (reg.is_potential(name)) {
            ctx.potentials.push_back(reg.make_potential(name, section_cfg));
        } else if (reg.is_coulomb(name)) {
            if (ctx.coulomb)
                throw std::runtime_error(
                    "Multiple Coulomb solver sections in md.dat — only one allowed");
            ctx.coulomb = reg.make_coulomb(name, section_cfg);
        } else if (reg.is_hook(name)) {
            ctx.hooks.push_back(reg.make_hook(name, section_cfg));
        } else {
            std::fprintf(stderr,
                "Warning: unknown section '%s' in md.dat (ignored).\n"
                "  Known names: %s\n",
                name.c_str(), reg.known_names().c_str());
        }
    }

    // -----------------------------------------------------------------------
    // 7. Default integrator (VelocityVerlet) if none was specified
    // -----------------------------------------------------------------------
    if (!integrator)
        integrator = reg.make_integrator("VelocityVerlet", Config{});

    // -----------------------------------------------------------------------
    // 8. Sort hooks by (hook_point, priority)
    // -----------------------------------------------------------------------
    ctx.sort_hooks();

    // -----------------------------------------------------------------------
    // 9. Two-pass neighbor-list initialisation
    // -----------------------------------------------------------------------
    // Pass A: use a generous initial cutoff so bind_all() can load parameters
    //         and report the real cutoffs (e.g. BOP, PetersT).
    nl.set_cutoff(10.0 + cutoff_add);
    nl.set_verlet_shell(cutoff_add * 0.5);
    nl.update(sys);

    integrator->bind_to(sys, nl);

    // Pass A: potentials only — loads BOP/DFTB parameters so cutoffs are known.
    ctx.bind_potentials_only();

    // Determine the correct NL cutoff from potentials + hook requirements.
    // PetersT cutoff comes from set_potential_cutoff (not yet called), but it
    // defaults to the potential cutoff, so ctx.max_cutoff() covers it.
    double actual_cutoff = ctx.max_cutoff();

    // Pass B: rebuild NL at the correct cutoff, then do the full bind (once).
    if (actual_cutoff + cutoff_add < nl.cutoff() - 0.05) {
        nl.set_cutoff(actual_cutoff + cutoff_add);
        nl.set_verlet_shell(cutoff_add * 0.5);
        nl.update(sys);
    }
    // Full bind: hooks open files, store reference positions, etc.
    ctx.bind_all();

    // -----------------------------------------------------------------------
    // 10. Print setup summary
    // -----------------------------------------------------------------------
    std::fprintf(stdout, "Read %zu atoms\n", sys.num_atoms());
    std::fprintf(stdout, "Unit mode: %s\n", fs_mode ? "eV/A/fs" : "eV/A");
    {
        std::set<int> elems;
        for (size_t i = 0; i < sys.num_atoms(); ++i)
            elems.insert(sys.atomic_number(i));
        std::fprintf(stdout, "Elements:");
        for (int Z : elems)
            std::fprintf(stdout, " %s", Z_to_symbol(Z).c_str());
        std::fprintf(stdout, "\n");
    }
    std::fprintf(stdout, "NL cutoff: %.3f A\n", nl.cutoff());

    // -----------------------------------------------------------------------
    // 11. Compute initial forces
    // -----------------------------------------------------------------------
    ctx.compute_forces();

    // -----------------------------------------------------------------------
    // 12. Status printing helpers
    // -----------------------------------------------------------------------
    int dof = n_dof(sys, /*skip_frozen=*/true);

    int nout = 0;

    auto print_header = [&]() {
        std::printf("%10s  %10s  %10s  %12s  %12s  %12s  %12s  %10s  %12s\n",
            "it",
            (std::string("t[") + time_label + "]").c_str(),
            (std::string("dt[") + time_label + "]").c_str(),
            "ekin[eV]", "epot[eV]", "etot[eV]", "fmax[eV/A]",
            "T[K]",
            (std::string("P[") + pressure_label + "]").c_str());
    };

    auto print_status = [&](int& nout_ref) {
        double ekin = kinetic_energy(sys);
        double T_K  = (dof > 0) ? 2.0 * ekin / (dof * kB_eV) : 0.0;
        double P    = pressure(ctx.results.virial, sys);
        double fm   = fmax(sys);
        double ti_d = ctx.time * ctx.time_display_scale;
        double dt_d = ctx.dt   * ctx.time_display_scale;

        if (nout_ref % 20 == 0)
            print_header();

        std::printf(
            "%10ld  %10.1f  %10.6f  %12.5E  %12.5E  %12.5E  %12.5E  %10.3f"
            "  %12.3E\n",
            static_cast<long>(ctx.step),
            ti_d, dt_d,
            ekin, ctx.results.energy, ekin + ctx.results.energy,
            fm, T_K, P);
        std::fflush(stdout);
        ++nout_ref;
    };

    // -----------------------------------------------------------------------
    // 13. Stopping predicate
    // -----------------------------------------------------------------------
    auto should_stop = [&]() -> bool {
        if (has_max_iter && static_cast<long>(ctx.step) >= max_iter)
            return true;
        if (has_max_time && ctx.time >= max_time_internal)
            return true;
        return false;
    };

    // -----------------------------------------------------------------------
    // 14. Initial status line
    // -----------------------------------------------------------------------
    print_status(nout);

    // -----------------------------------------------------------------------
    // 15. Main MD loop
    // -----------------------------------------------------------------------
    while (!should_stop()) {
        // Checkpoint before each step (at the requested frequency)
        if (ctx.step > 0 && static_cast<int>(ctx.step) % file_freq == 0) {
            std::string ckpt = ((ctx.step / file_freq) % 2 == 0)
                ? "atomsA.out" : "atomsB.out";
            write_atoms_dat(ckpt, sys, atoms_data.unit_mode);
        }

        bool cont = integrator->step(ctx);
        if (!cont) break;

        if (static_cast<int>(ctx.step) % scr_freq == 0)
            print_status(nout);
    }

    // -----------------------------------------------------------------------
    // 16. Final status, output, and DONE marker
    // -----------------------------------------------------------------------
    print_status(nout);

    write_atoms_dat("atoms.out", sys, atoms_data.unit_mode);
    std::fprintf(stdout, "Wrote atoms.out\n");

    { std::ofstream done("DONE"); }

    return 0;
}
