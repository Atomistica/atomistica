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
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>

#include "../include/atomistica/config.hpp"
#include "../include/atomistica/core/atomic_system.hpp"
#include "../include/atomistica/core/neighbor_list.hpp"
#include "../include/atomistica/integrators/verlet.hpp"
#include "../include/atomistica/integrators/thermostats.hpp"
#include "../include/atomistica/potentials/potential_base.hpp"
#include "../include/atomistica/potentials/bop/tersoff.hpp"
#include "../include/atomistica/potentials/bop/brenner.hpp"
#include "../include/atomistica/potentials/bop/kumagai.hpp"
#include "../include/atomistica/potentials/bop/juslin.hpp"
#include "../include/atomistica/potentials/bop/rebo2.hpp"

#include "config.hpp"
#include "atoms_io.hpp"
#include "peters_t.hpp"

using namespace atomistica;

// Collect unique element atomic numbers in a system
static std::set<int> unique_elements(const AtomicSystem& sys) {
    std::set<int> elems;
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        elems.insert(sys.atomic_number(i));
    return elems;
}

// Create a potential based on md.dat configuration and element types
static std::unique_ptr<Potential> create_potential(const Config& cfg,
                                                    const AtomicSystem& sys) {
    auto elems = unique_elements(sys);

    // Tersoff
    if (cfg.has_section("Tersoff")) {
        auto* w = new PotentialWrapper<Tersoff<false>>();
        std::string param_set;
        if (elems.count(14)) {
            param_set = "Tersoff_PRB_39_5566_Si_C";
        } else if (elems.count(13) && elems.count(7)) {
            param_set = "Goumri_Said_ChemPhys_302_135_Al_N";
        } else if ((elems.count(5) || elems.count(6) || elems.count(7)) &&
                   elems.size() <= 3) {
            param_set = "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N";
        } else {
            delete w;
            throw std::runtime_error(
                "Tersoff: cannot auto-select parameters for the given elements");
        }
        w->get().load_parameters(param_set);
        std::fprintf(stdout, "Tersoff: loaded parameter set '%s'\n",
                     param_set.c_str());
        return std::unique_ptr<Potential>(w);
    }

    // TersoffScr
    if (cfg.has_section("TersoffScr")) {
        auto* w = new PotentialWrapper<Tersoff<true>>();
        w->get().load_parameters("Tersoff_PRB_39_5566_Si_C");
        return std::unique_ptr<Potential>(w);
    }

    // Brenner
    if (cfg.has_section("Brenner")) {
        auto* w = new PotentialWrapper<Brenner<false>>();
        std::string param_set;
        if (elems.count(14) && elems.count(6)) {
            param_set = "Erhart_PRB_71_035211_SiC";
        } else if (elems.count(78) && elems.count(6)) {
            param_set = "Albe_PRB_65_195124_PtC";
        } else if (elems.count(26) && elems.count(6)) {
            param_set = "Henriksson_PRB_79_144107_FeC";
        } else if (elems.count(6) && elems.size() <= 2) {
            param_set = "Brenner_PRB_42_9458_C_II";
        } else {
            delete w;
            throw std::runtime_error(
                "Brenner: cannot auto-select parameters for the given elements");
        }
        w->get().load_parameters(param_set);
        std::fprintf(stdout, "Brenner: loaded parameter set '%s'\n",
                     param_set.c_str());
        return std::unique_ptr<Potential>(w);
    }

    // Kumagai
    if (cfg.has_section("Kumagai")) {
        auto* w = new PotentialWrapper<Kumagai<false>>();
        w->get().load_parameters("Kumagai_CompMaterSci_39_457_Si");
        return std::unique_ptr<Potential>(w);
    }

    // Juslin
    if (cfg.has_section("Juslin")) {
        auto* w = new PotentialWrapper<Juslin<false>>();
        w->get().load_parameters("Juslin_JAP_98_123520_WCH");
        return std::unique_ptr<Potential>(w);
    }

    // REBO2 / Rebo2
    if (cfg.has_section("Rebo2") || cfg.has_section("REBO2") ||
        cfg.has_section("rebo2")) {
        auto* w = new PotentialWrapper<REBO2>();
        w->get().load_default_parameters();
        std::fprintf(stdout, "REBO2: loaded default parameters\n");
        return std::unique_ptr<Potential>(w);
    }

    throw std::runtime_error(
        "No supported potential found in md.dat. "
        "Supported: Tersoff, TersoffScr, Brenner, Kumagai, Juslin, Rebo2");
}

// Wrap all atom positions back into the simulation cell
static void wrap_positions(AtomicSystem& sys) {
    for (size_t i = 0; i < sys.num_atoms(); ++i) {
        Vec3 r = sys.wrap_position(sys.position(i));
        sys.set_position(i, r);
    }
}

// Compute fmax (maximum force magnitude)
static double compute_fmax(const AtomicSystem& sys) {
    double fmax2 = 0.0;
    for (size_t i = 0; i < sys.num_atoms(); ++i) {
        double f2 = sys.forces().col(i).matrix().squaredNorm();
        fmax2 = std::max(fmax2, f2);
    }
    return std::sqrt(fmax2);
}

// Compute kinetic energy (eV)
static double kinetic_energy(const AtomicSystem& sys) {
    double ekin = 0.0;
    for (size_t i = 0; i < sys.num_atoms(); ++i) {
        double m = sys.mass(i);
        Vec3 v = sys.velocity(i);
        ekin += 0.5 * m * v.squaredNorm();
    }
    return ekin;
}

// Compute kinetic virial: W_kin = sum_i m_i * v_i ⊗ v_i
static Mat3 kinetic_virial(const AtomicSystem& sys) {
    Mat3 W = Mat3::Zero();
    for (size_t i = 0; i < sys.num_atoms(); ++i) {
        double m = sys.mass(i);
        Vec3 v = sys.velocity(i);
        W += m * v * v.transpose();
    }
    return W;
}

// Print column headers (every 10 output lines)
static void print_header(const char* time_label, const char* pressure_label) {
    std::printf("%10s  %10s  %10s  %12s  %12s  %12s  %12s  %10s  %12s\n",
        "it",
        (std::string("t[") + time_label + "]").c_str(),
        (std::string("dt[") + time_label + "]").c_str(),
        "ekin[eV]", "epot[eV]", "etot[eV]", "fmax[eV/A]",
        "T[K]",
        (std::string("P[") + pressure_label + "]").c_str());
}

// Print one status line
static void print_status(long it, double ti, double dt,
                          double ekin, double epot, double fmax,
                          double T, double P,
                          const char* time_label, const char* pressure_label,
                          int& nout) {
    if (nout % 10 == 0) {
        print_header(time_label, pressure_label);
    }
    std::printf("%10ld  %10.1f  %10.6f  %12.5E  %12.5E  %12.5E  %12.5E  %10.3f  %12.3E\n",
        it, ti, dt, ekin, epot, ekin + epot, fmax, T, P);
    ++nout;
}

int main() {
    // -----------------------------------------------------------------------
    // 1. Parse configuration
    // -----------------------------------------------------------------------
    Config cfg;
    cfg.parse_file("md.dat");

    std::string sou = cfg.get_string("system_of_units", "eV/A");
    bool fs_mode = (sou.find("fs") != std::string::npos);  // eV/A/fs

    double dt_raw    = cfg.get_double("dt", 0.1);
    double max_time  = cfg.get_double("max_time", 100.0);
    int    scr_freq  = cfg.get_int("scr_freq", 10);
    int    file_freq = cfg.get_int("file_freq", 10);
    double cutoff_add = cfg.get_double("cutoff_add", 0.5);

    // Convert dt to internal units
    const double sqrt_c = std::sqrt(AMU_AFSQ_PER_EV);  // ≈ 10.18
    double dt = fs_mode ? dt_raw / sqrt_c : dt_raw;
    double max_time_internal = fs_mode ? max_time / sqrt_c : max_time;

    const char* time_label     = fs_mode ? "fs" : "10fs";
    const char* pressure_label = "eV/A^3";

    // -----------------------------------------------------------------------
    // 2. Read atoms
    // -----------------------------------------------------------------------
    AtomsData atoms_data = read_atoms_dat("atoms.dat");
    AtomicSystem& sys = atoms_data.system;

    std::fprintf(stdout, "Read %zu atoms\n", sys.num_atoms());
    std::fprintf(stdout, "Unit mode: %s\n", fs_mode ? "eV/A/fs" : "eV/A");

    // Print elements
    {
        auto elems = unique_elements(sys);
        std::fprintf(stdout, "Elements:");
        for (int Z : elems) std::fprintf(stdout, " %s", Z_to_symbol(Z).c_str());
        std::fprintf(stdout, "\n");
    }

    int n_dof = 3 * static_cast<int>(sys.num_atoms()) - 3;  // remove COM

    // -----------------------------------------------------------------------
    // 3. Create potential
    // -----------------------------------------------------------------------
    std::unique_ptr<Potential> pot = create_potential(cfg, sys);

    // -----------------------------------------------------------------------
    // 4. Set up neighbor list
    // -----------------------------------------------------------------------
    double nl_cutoff = pot->cutoff() + cutoff_add;
    NeighborList nl;
    nl.set_cutoff(nl_cutoff);
    nl.set_verlet_shell(cutoff_add * 0.5);

    nl.update(sys);
    pot->bind_to(sys, nl);

    // -----------------------------------------------------------------------
    // 5. Thermostat
    // -----------------------------------------------------------------------
    std::unique_ptr<PetersT> thermostat;
    double peters_cutoff_nl = nl_cutoff;  // default: use potential cutoff

    if (cfg.has_section("PetersT")) {
        const Config& pc = cfg.section("PetersT");
        double T      = pc.get_double("T", 300.0);
        double gamma  = pc.get_double("gamma", 1.0);
        double cutoff = pc.get_double("cutoff", pot->cutoff());

        // Convert gamma for eV_A_fs mode
        if (fs_mode) gamma /= sqrt_c;

        thermostat = std::make_unique<PetersT>(T, gamma, cutoff);
        peters_cutoff_nl = cutoff;

        // Ensure NL is built with enough cutoff for thermostat too
        if (peters_cutoff_nl + cutoff_add > nl.cutoff()) {
            nl.set_cutoff(peters_cutoff_nl + cutoff_add);
            nl.invalidate();
            nl.update(sys);
        }

        std::fprintf(stdout, "PetersT thermostat: T=%.1f K, gamma=%.4f, cutoff=%.2f\n",
                     T, gamma, cutoff);
    }

    // -----------------------------------------------------------------------
    // 6. Compute initial forces
    // -----------------------------------------------------------------------
    sys.zero_forces();
    nl.update(sys);
    PotentialResults results = pot->compute(sys, nl);

    // -----------------------------------------------------------------------
    // 7. Main loop
    // -----------------------------------------------------------------------
    VelocityVerlet verlet;
    verlet.set_timestep(dt);

    double ti = 0.0;
    long it = 0;
    int nout = 0;

    // Print initial status
    {
        double ekin = kinetic_energy(sys);
        double T_K  = (n_dof > 0) ? 2.0 * ekin / (n_dof * kB_eV_K) : 0.0;
        double vol  = sys.volume();
        Mat3 Wkin   = kinetic_virial(sys);
        double P    = (vol > 0.0)
            ? (results.virial + Wkin).trace() / (3.0 * vol) : 0.0;
        double fmax = compute_fmax(sys);

        // Display time depends on mode
        double ti_disp = fs_mode ? ti * sqrt_c : ti;
        double dt_disp = fs_mode ? dt * sqrt_c : dt;

        print_status(it, ti_disp, dt_disp, ekin, results.energy, fmax,
                     T_K, P, time_label, pressure_label, nout);
    }

    while (ti < max_time_internal) {
        ++it;

        // Checkpoint to alternating files
        if (it == 1 || it % file_freq == 0) {
            std::string fname = ((it / file_freq) % 2 == 0)
                ? "atomsA.out" : "atomsB.out";
            write_atoms_dat(fname, sys, atoms_data.unit_mode);
        }

        // Verlet step 1: update v(+dt/2) and r(+dt)
        verlet.step1(sys, sys.forces().matrix().transpose());
        wrap_positions(sys);

        // Rebuild NL if needed and compute new forces
        nl.update(sys);
        sys.zero_forces();
        results = pot->compute(sys, nl);

        // Verlet step 2: complete velocity update
        verlet.step2(sys, sys.forces().matrix().transpose());

        // Apply thermostat
        if (thermostat) {
            thermostat->apply(sys, nl, dt);
        }

        ti += dt;

        // Screen output
        if (it % scr_freq == 0) {
            double ekin = kinetic_energy(sys);
            double T_K  = (n_dof > 0) ? 2.0 * ekin / (n_dof * kB_eV_K) : 0.0;
            double vol  = sys.volume();
            Mat3 Wkin   = kinetic_virial(sys);
            double P    = (vol > 0.0)
                ? (results.virial + Wkin).trace() / (3.0 * vol) : 0.0;
            double fmax = compute_fmax(sys);

            double ti_disp = fs_mode ? ti * sqrt_c : ti;
            double dt_disp = fs_mode ? dt * sqrt_c : dt;

            print_status(it, ti_disp, dt_disp, ekin, results.energy, fmax,
                         T_K, P, time_label, pressure_label, nout);
        }
    }

    // -----------------------------------------------------------------------
    // 8. Final output
    // -----------------------------------------------------------------------
    // Final status
    {
        double ekin = kinetic_energy(sys);
        double T_K  = (n_dof > 0) ? 2.0 * ekin / (n_dof * kB_eV_K) : 0.0;
        double vol  = sys.volume();
        Mat3 Wkin   = kinetic_virial(sys);
        double P    = (vol > 0.0)
            ? (results.virial + Wkin).trace() / (3.0 * vol) : 0.0;
        double fmax = compute_fmax(sys);

        double ti_disp = fs_mode ? ti * sqrt_c : ti;
        double dt_disp = fs_mode ? dt * sqrt_c : dt;

        print_status(it, ti_disp, dt_disp, ekin, results.energy, fmax,
                     T_K, P, time_label, pressure_label, nout);
    }

    write_atoms_dat("atoms.out", sys, atoms_data.unit_mode);
    std::fprintf(stdout, "Wrote atoms.out\n");

    // Create DONE file
    {
        std::ofstream done("DONE");
    }

    return 0;
}
