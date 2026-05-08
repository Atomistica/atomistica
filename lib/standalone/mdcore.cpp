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
#include "../include/atomistica/integrators/verlet.hpp"
#include "../include/atomistica/integrators/thermostats.hpp"
#include "../include/atomistica/potentials/potential_base.hpp"
#include "../include/atomistica/potentials/bop/tersoff.hpp"
#include "../include/atomistica/potentials/bop/brenner.hpp"
#include "../include/atomistica/potentials/bop/kumagai.hpp"
#include "../include/atomistica/potentials/bop/juslin.hpp"
#include "../include/atomistica/potentials/bop/rebo2.hpp"
#include "../include/atomistica/tightbinding/dftb.hpp"

#include "config.hpp"
#include "atoms_io.hpp"
#include "peters_t.hpp"

using namespace atomistica;

// ---------------------------------------------------------------------------
// DFTBPotential: wraps tb::DFTB to implement the Potential virtual interface
// ---------------------------------------------------------------------------
class DFTBPotential : public Potential {
public:
    explicit DFTBPotential(const std::string& skf_path, bool enable_scc,
                           const tb::SCCParams& scc_params,
                           const tb::SolverParams& solver_params)
        : dftb_(skf_path, enable_scc)
    {
        dftb_.set_scc_params(scc_params);
        dftb_.set_solver_params(solver_params);
    }

    Scalar cutoff() const override { return dftb_.cutoff(); }

    // Initialise DFTB for the given system (call before first compute).
    // Sets element list; must be called once after system is known.
    void pre_init(const AtomicSystem& system) {
        dftb_.init(system);
    }

    void bind_to(AtomicSystem& system, NeighborList& /*nl*/) override {
        // Init is done via pre_init; re-init here is a no-op guard
        if (dftb_.cutoff() < 1e-6)
            dftb_.init(system);
    }

    PotentialResults compute(AtomicSystem& system,
                             NeighborList& neighbors,
                             bool /*compute_forces*/ = true,
                             bool compute_virial = true) override {
        int nat = static_cast<int>(system.num_atoms());
        MatX3 forces(nat, 3);
        forces.setZero();

        PotentialResults results;

        if (compute_virial) {
            Mat3 stress = Mat3::Zero();
            results.energy = dftb_.compute_with_stress(system, neighbors,
                                                        forces, stress);
            // compute_with_stress returns stress = virial/volume;
            // we store virial = stress * volume to match other potentials
            results.virial = stress * system.volume();
        } else {
            results.energy = dftb_.compute(system, neighbors, forces);
        }

        // Copy forces (N×3) into system.forces() (3×N)
        for (int i = 0; i < nat; ++i)
            system.forces().col(i) = forces.row(i).transpose();

        return results;
    }

private:
    tb::DFTB dftb_;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Collect unique element atomic numbers in a system
static std::set<int> unique_elements(const AtomicSystem& sys) {
    std::set<int> elems;
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        elems.insert(sys.atomic_number(i));
    return elems;
}

// Auto-select an SKF directory for the given element set.
// Checks (in order): TBPARAM env var, then a set of known ~/Databases paths.
static std::string auto_select_skf_path(const std::set<int>& elems) {
    // 1. Explicit env var (matches Fortran TBPARAM convention)
    const char* tbparam = std::getenv("TBPARAM");
    if (tbparam && *tbparam) return tbparam;

    // 2. Try standard database directories in $HOME/Databases
    const char* home = std::getenv("HOME");
    if (!home) home = "";
    std::string home_str(home);

    // Candidate databases in preference order
    struct Candidate {
        std::string subdir;
        std::set<int> supported_Z;  // empty means "try anyway"
    };
    std::vector<Candidate> candidates = {
        // mio-1-1: C(6) H(1) N(7) O(8) S(16) P(15)
        {"mio-1-1",   {1, 6, 7, 8, 15, 16}},
        // 3ob-3-1: broader organic + Mg, Zn, Ca, Na, K, Cl, Br, etc.
        {"3ob-3-1",   {1, 6, 7, 8, 11, 12, 15, 16, 17, 19, 20, 30, 35}},
        // pbc-0-3: solid-state: C H N O F Si Fe Ni Cu
        {"pbc-0-3",   {1, 6, 7, 8, 9, 14, 26, 28, 29}},
        // matsci-0-3: broader solid state
        {"matsci-0-3", {}},
    };

    for (auto& cand : candidates) {
        std::string path = home_str + "/Databases/" + cand.subdir;
        // Check if directory exists (try to open a dummy path)
        std::ifstream probe(path + "/.check_exists_dummy");
        // probe.good() will be false, but the directory check via stat is
        // simpler via a known file. Instead: try to open a plausible skf file.
        // If elems set matches, try the path regardless.
        bool matches = cand.supported_Z.empty();
        if (!matches) {
            matches = true;
            for (int Z : elems) {
                if (!cand.supported_Z.count(Z)) { matches = false; break; }
            }
        }
        if (!matches) continue;

        // Verify the directory actually exists by trying one possible file
        // Use a generic check: directory must contain at least one .skf file
        // We do a simple existence check via trying to open a known-format path
        // for the first element.
        if (!elems.empty()) {
            int Z1 = *elems.begin();
            // get_element_symbol is internal; use atoms_io's Z_to_symbol
            std::string sym = Z_to_symbol(Z1);
            std::string test_file = path + "/" + sym + "-" + sym + ".skf";
            std::ifstream f(test_file);
            if (f.good()) return path;
        } else {
            return path;  // no elements yet — will fail later
        }
    }

    return "";  // not found — caller will produce a useful error
}

// Create a potential based on md.dat configuration and element types.
// For DFTB, pre_init is called here so cutoff() is valid after return.
static std::unique_ptr<Potential> create_potential(const Config& cfg,
                                                    const AtomicSystem& sys) {
    auto elems = unique_elements(sys);

    // ----- Tersoff -----
    if (cfg.has_section("Tersoff")) {
        auto* w = new PotentialWrapper<Tersoff<false>>();
        std::string param_set;
        if (elems.count(14))
            param_set = "Tersoff_PRB_39_5566_Si_C";
        else if (elems.count(13) && elems.count(7))
            param_set = "Goumri_Said_ChemPhys_302_135_Al_N";
        else if ((elems.count(5) || elems.count(6) || elems.count(7))
                 && elems.size() <= 3)
            param_set = "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N";
        else {
            delete w;
            throw std::runtime_error(
                "Tersoff: cannot auto-select parameters for the given elements");
        }
        w->get().load_parameters(param_set);
        std::fprintf(stdout, "Tersoff: loaded parameter set '%s'\n",
                     param_set.c_str());
        return std::unique_ptr<Potential>(w);
    }

    // ----- TersoffScr -----
    if (cfg.has_section("TersoffScr")) {
        auto* w = new PotentialWrapper<Tersoff<true>>();
        w->get().load_parameters("Tersoff_PRB_39_5566_Si_C");
        return std::unique_ptr<Potential>(w);
    }

    // ----- Brenner -----
    if (cfg.has_section("Brenner")) {
        auto* w = new PotentialWrapper<Brenner<false>>();
        std::string param_set;
        if (elems.count(14) && elems.count(6))
            param_set = "Erhart_PRB_71_035211_SiC";
        else if (elems.count(78) && elems.count(6))
            param_set = "Albe_PRB_65_195124_PtC";
        else if (elems.count(26) && elems.count(6))
            param_set = "Henriksson_PRB_79_144107_FeC";
        else if (elems.count(6) && elems.size() <= 2)
            param_set = "Brenner_PRB_42_9458_C_II";
        else {
            delete w;
            throw std::runtime_error(
                "Brenner: cannot auto-select parameters for the given elements");
        }
        w->get().load_parameters(param_set);
        std::fprintf(stdout, "Brenner: loaded parameter set '%s'\n",
                     param_set.c_str());
        return std::unique_ptr<Potential>(w);
    }

    // ----- Kumagai -----
    if (cfg.has_section("Kumagai")) {
        auto* w = new PotentialWrapper<Kumagai<false>>();
        w->get().load_parameters("Kumagai_CompMaterSci_39_457_Si");
        return std::unique_ptr<Potential>(w);
    }

    // ----- Juslin -----
    if (cfg.has_section("Juslin")) {
        auto* w = new PotentialWrapper<Juslin<false>>();
        w->get().load_parameters("Juslin_JAP_98_123520_WCH");
        return std::unique_ptr<Potential>(w);
    }

    // ----- REBO2 / Rebo2 -----
    if (cfg.has_section("Rebo2") || cfg.has_section("REBO2") ||
        cfg.has_section("rebo2")) {
        auto* w = new PotentialWrapper<REBO2>();
        w->get().load_default_parameters();
        std::fprintf(stdout, "REBO2: loaded default parameters\n");
        return std::unique_ptr<Potential>(w);
    }

    // ----- TightBinding (DFTB) -----
    if (cfg.has_section("TightBinding") || cfg.has_section("tightbinding")) {
        const Config& tb_cfg = cfg.has_section("TightBinding")
                               ? cfg.section("TightBinding")
                               : cfg.section("tightbinding");

        // SKF database path: md.dat key > TBPARAM env > auto-detect
        std::string skf_path = tb_cfg.get_string("database", "");
        if (skf_path.empty())
            skf_path = auto_select_skf_path(elems);
        if (skf_path.empty())
            throw std::runtime_error(
                "DFTB: cannot find SKF database. Set 'database' in md.dat "
                "TightBinding section, or set TBPARAM environment variable, "
                "or place SKF files in ~/Databases/mio-1-1/");

        // SCC parameters
        bool enable_scc = tb_cfg.has_section("SCC");
        tb::SCCParams scc_params;
        if (enable_scc) {
            const Config& scc = tb_cfg.section("SCC");
            if (!scc.get_string("dq_crit", "").empty())
                scc_params.convergence_threshold =
                    scc.get_double("dq_crit", 1e-4);
            if (!scc.get_string("maximum_iterations", "").empty())
                scc_params.max_iterations =
                    scc.get_int("maximum_iterations", 200);
            if (!scc.get_string("mixing", "").empty())
                scc_params.mixing_parameter = scc.get_double("mixing", 0.2);
            if (!scc.get_string("andersen_memory", "").empty())
                scc_params.anderson_memory =
                    scc.get_int("andersen_memory", 3);
        }

        // Solver parameters (from SolverLAPACK or SolverCP sections)
        tb::SolverParams solver_params;
        for (const char* sec : {"SolverLAPACK", "SolverCP", "Solver"}) {
            if (tb_cfg.has_section(sec)) {
                const Config& slv = tb_cfg.section(sec);
                if (!slv.get_string("electronic_T", "").empty())
                    solver_params.electronic_temperature =
                        slv.get_double("electronic_T", 0.01);
                break;
            }
        }

        auto* dftb_pot = new DFTBPotential(skf_path, enable_scc,
                                            scc_params, solver_params);
        // Pre-init to discover element basis → cutoff becomes valid
        dftb_pot->pre_init(sys);

        std::fprintf(stdout,
                     "DFTB: SKF path '%s', SCC=%s, electronic_T=%.4f eV\n",
                     skf_path.c_str(),
                     enable_scc ? "yes" : "no",
                     solver_params.electronic_temperature);
        return std::unique_ptr<Potential>(dftb_pot);
    }

    throw std::runtime_error(
        "No supported potential found in md.dat. "
        "Supported: Tersoff, TersoffScr, Brenner, Kumagai, Juslin, "
        "Rebo2, TightBinding");
}

// Wrap all atom positions back into the simulation cell
static void wrap_positions(AtomicSystem& sys) {
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        sys.set_position(i, sys.wrap_position(sys.position(i)));
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
    if (nout % 10 == 0)
        print_header(time_label, pressure_label);
    std::printf(
        "%10ld  %10.1f  %10.6f  %12.5E  %12.5E  %12.5E  %12.5E  %10.3f"
        "  %12.3E\n",
        it, ti, dt, ekin, epot, ekin + epot, fmax, T, P);
    std::fflush(stdout);
    ++nout;
}

int main() {
    // -----------------------------------------------------------------------
    // 1. Parse configuration
    // -----------------------------------------------------------------------
    Config cfg;
    cfg.parse_file("md.dat");

    std::string sou = cfg.get_string("system_of_units", "eV/A");
    bool fs_mode = (sou.find("fs") != std::string::npos);

    double dt_raw      = cfg.get_double("dt", 0.1);
    double max_time    = cfg.get_double("max_time", -1.0);  // -1 = not set
    long   max_iter    = static_cast<long>(cfg.get_int("n_iterations", -1));
    int    scr_freq    = cfg.get_int("scr_freq", 10);
    int    file_freq   = cfg.get_int("file_freq", 10);
    double cutoff_add  = cfg.get_double("cutoff_add", 0.5);

    // If neither max_time nor n_iterations specified, default to 100 time units
    bool has_max_time = (max_time > 0.0);
    bool has_max_iter = (max_iter > 0);
    if (!has_max_time && !has_max_iter) {
        max_time = 100.0;
        has_max_time = true;
    }

    // Convert dt and max_time to internal units
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
    {
        auto elems = unique_elements(sys);
        std::fprintf(stdout, "Elements:");
        for (int Z : elems) std::fprintf(stdout, " %s", Z_to_symbol(Z).c_str());
        std::fprintf(stdout, "\n");
    }

    int n_dof = 3 * static_cast<int>(sys.num_atoms()) - 3;

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
    double peters_cutoff_nl = nl_cutoff;

    if (cfg.has_section("PetersT")) {
        const Config& pc = cfg.section("PetersT");
        double T_therm = pc.get_double("T", 300.0);
        double gamma   = pc.get_double("gamma", 1.0);
        double cutoff  = pc.get_double("cutoff", pot->cutoff());

        if (fs_mode) gamma /= sqrt_c;

        thermostat = std::make_unique<PetersT>(T_therm, gamma, cutoff);
        peters_cutoff_nl = cutoff;

        if (peters_cutoff_nl + cutoff_add > nl.cutoff()) {
            nl.set_cutoff(peters_cutoff_nl + cutoff_add);
            nl.invalidate();
            nl.update(sys);
        }

        std::fprintf(stdout,
                     "PetersT thermostat: T=%.1f K, gamma=%.4f, cutoff=%.2f\n",
                     T_therm, gamma, cutoff);
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

    // Initial status
    {
        double ekin = kinetic_energy(sys);
        double T_K  = (n_dof > 0) ? 2.0 * ekin / (n_dof * kB_eV_K) : 0.0;
        double vol  = sys.volume();
        Mat3   Wkin = kinetic_virial(sys);
        double P    = (vol > 0.0)
            ? (results.virial + Wkin).trace() / (3.0 * vol) : 0.0;
        double fmax = compute_fmax(sys);
        double ti_d = fs_mode ? ti * sqrt_c : ti;
        double dt_d = fs_mode ? dt * sqrt_c : dt;
        print_status(it, ti_d, dt_d, ekin, results.energy, fmax, T_K, P,
                     time_label, pressure_label, nout);
    }

    auto should_stop = [&]() -> bool {
        if (has_max_time && ti >= max_time_internal) return true;
        if (has_max_iter && it >= max_iter) return true;
        return false;
    };

    while (!should_stop()) {
        ++it;

        // Checkpoint to alternating files
        if (it == 1 || it % file_freq == 0) {
            std::string fname = ((it / file_freq) % 2 == 0)
                ? "atomsA.out" : "atomsB.out";
            write_atoms_dat(fname, sys, atoms_data.unit_mode);
        }

        // Verlet step 1: v(+dt/2), r(+dt)
        verlet.step1(sys, sys.forces().matrix().transpose());
        wrap_positions(sys);

        // Forces at new positions
        nl.update(sys);
        sys.zero_forces();
        results = pot->compute(sys, nl);

        // Verlet step 2: v(+dt)
        verlet.step2(sys, sys.forces().matrix().transpose());

        // Apply thermostat
        if (thermostat)
            thermostat->apply(sys, nl, dt);

        ti += dt;

        if (it % scr_freq == 0) {
            double ekin = kinetic_energy(sys);
            double T_K  = (n_dof > 0) ? 2.0 * ekin / (n_dof * kB_eV_K) : 0.0;
            double vol  = sys.volume();
            Mat3   Wkin = kinetic_virial(sys);
            double P    = (vol > 0.0)
                ? (results.virial + Wkin).trace() / (3.0 * vol) : 0.0;
            double fmax = compute_fmax(sys);
            double ti_d = fs_mode ? ti * sqrt_c : ti;
            double dt_d = fs_mode ? dt * sqrt_c : dt;
            print_status(it, ti_d, dt_d, ekin, results.energy, fmax, T_K, P,
                         time_label, pressure_label, nout);
        }
    }

    // -----------------------------------------------------------------------
    // 8. Final output
    // -----------------------------------------------------------------------
    {
        double ekin = kinetic_energy(sys);
        double T_K  = (n_dof > 0) ? 2.0 * ekin / (n_dof * kB_eV_K) : 0.0;
        double vol  = sys.volume();
        Mat3   Wkin = kinetic_virial(sys);
        double P    = (vol > 0.0)
            ? (results.virial + Wkin).trace() / (3.0 * vol) : 0.0;
        double fmax = compute_fmax(sys);
        double ti_d = fs_mode ? ti * sqrt_c : ti;
        double dt_d = fs_mode ? dt * sqrt_c : dt;
        print_status(it, ti_d, dt_d, ekin, results.energy, fmax, T_K, P,
                     time_label, pressure_label, nout);
    }

    write_atoms_dat("atoms.out", sys, atoms_data.unit_mode);
    std::fprintf(stdout, "Wrote atoms.out\n");

    { std::ofstream done("DONE"); }

    return 0;
}
