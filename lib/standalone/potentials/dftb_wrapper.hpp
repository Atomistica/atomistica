// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <cstdlib>
#include <fstream>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../../include/atomistica/potentials/potential_base.hpp"
#include "../../include/atomistica/tightbinding/dftb.hpp"

#include "../config.hpp"
#include "../registry.hpp"
#include "../atoms_io.hpp"

namespace atomistica {

static std::set<int> dftb_unique_elements(const AtomicSystem& sys) {
    std::set<int> elems;
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        elems.insert(sys.atomic_number(i));
    return elems;
}

static std::string auto_select_skf_path(const std::set<int>& elems) {
    const char* tbparam = std::getenv("TBPARAM");
    if (tbparam && *tbparam) return tbparam;

    const char* home = std::getenv("HOME");
    if (!home) home = "";
    std::string home_str(home);

    struct Candidate {
        std::string subdir;
        std::set<int> supported_Z;
    };
    std::vector<Candidate> candidates = {
        {"mio-1-1",    {1, 6, 7, 8, 15, 16}},
        {"3ob-3-1",    {1, 6, 7, 8, 11, 12, 15, 16, 17, 19, 20, 30, 35}},
        {"pbc-0-3",    {1, 6, 7, 8, 9, 14, 26, 28, 29}},
        {"matsci-0-3", {}},
    };

    for (auto& cand : candidates) {
        std::string path = home_str + "/Databases/" + cand.subdir;
        bool matches = cand.supported_Z.empty();
        if (!matches) {
            matches = true;
            for (int Z : elems) {
                if (!cand.supported_Z.count(Z)) { matches = false; break; }
            }
        }
        if (!matches) continue;

        if (!elems.empty()) {
            int Z1 = *elems.begin();
            std::string sym = Z_to_symbol(Z1);
            std::string test_file = path + "/" + sym + "-" + sym + ".skf";
            std::ifstream f(test_file);
            if (f.good()) return path;
        } else {
            return path;
        }
    }

    return "";
}

class DFTBPotential : public Potential {
    // Stored config for deferred DFTB construction (skf_path may need elements)
    std::string skf_path_;
    bool enable_scc_;
    tb::SCCParams scc_params_;
    tb::SolverParams solver_params_;

    std::unique_ptr<tb::DFTB> dftb_;
    bool initialized_ = false;

public:
    explicit DFTBPotential(const Config& cfg)
        : skf_path_(cfg.get_or<std::string>("database", ""))
        , enable_scc_(cfg.has_section("SCC"))
    {
        if (enable_scc_) {
            const Config& scc = cfg.section("SCC");
            if (scc.has_key("dq_crit"))
                scc_params_.convergence_threshold =
                    scc.get_or<double>("dq_crit", 1e-4);
            if (scc.has_key("maximum_iterations"))
                scc_params_.max_iterations =
                    scc.get_or<int>("maximum_iterations", 200);
            if (scc.has_key("mixing"))
                scc_params_.mixing_parameter =
                    scc.get_or<double>("mixing", 0.2);
            if (scc.has_key("andersen_memory"))
                scc_params_.anderson_memory =
                    scc.get_or<int>("andersen_memory", 3);
        }

        for (const char* sec : {"SolverLAPACK", "SolverCP", "Solver"}) {
            if (cfg.has_section(sec)) {
                const Config& slv = cfg.section(sec);
                if (slv.has_key("electronic_T"))
                    solver_params_.electronic_temperature =
                        slv.get_or<double>("electronic_T", 0.01);
                break;
            }
        }

        // If an explicit path was given, construct DFTB now
        if (!skf_path_.empty()) {
            dftb_ = std::make_unique<tb::DFTB>(skf_path_, enable_scc_);
            dftb_->set_scc_params(scc_params_);
            dftb_->set_solver_params(solver_params_);
        }
    }

    Scalar cutoff() const override {
        return dftb_ ? dftb_->cutoff() : 0.0;
    }

    void bind_to(AtomicSystem& sys, NeighborList& /*nl*/) override {
        if (!initialized_) {
            if (!dftb_) {
                auto elems = dftb_unique_elements(sys);
                skf_path_ = auto_select_skf_path(elems);
                if (skf_path_.empty())
                    throw std::runtime_error(
                        "DFTB: cannot find SKF database. Set 'database' in md.dat "
                        "TightBinding section, or set TBPARAM environment variable, "
                        "or place SKF files in ~/Databases/mio-1-1/");
                dftb_ = std::make_unique<tb::DFTB>(skf_path_, enable_scc_);
                dftb_->set_scc_params(scc_params_);
                dftb_->set_solver_params(solver_params_);
            }
            dftb_->init(sys);
            initialized_ = true;
        }
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool /*cf*/ = true, bool cv = true) override {
        int nat = static_cast<int>(sys.num_atoms());
        MatX3 forces(nat, 3);
        forces.setZero();

        PotentialResults results;

        if (cv) {
            Mat3 stress = Mat3::Zero();
            results.energy = dftb_->compute_with_stress(sys, nl, forces, stress);
            results.virial = stress * sys.volume();
        } else {
            results.energy = dftb_->compute(sys, nl, forces);
        }

        for (int i = 0; i < nat; ++i)
            sys.forces().col(i) = forces.row(i).transpose();

        return results;
    }
};

REGISTER_POTENTIAL("TightBinding", DFTBPotential)

} // namespace atomistica
