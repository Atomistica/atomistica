// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <set>
#include <stdexcept>
#include <string>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../../include/atomistica/potentials/potential_base.hpp"
#include "../../include/atomistica/potentials/bop/tersoff.hpp"
#include "../../include/atomistica/potentials/bop/brenner.hpp"
#include "../../include/atomistica/potentials/bop/kumagai.hpp"
#include "../../include/atomistica/potentials/bop/juslin.hpp"
#include "../../include/atomistica/potentials/bop/rebo2.hpp"

#include "../config.hpp"
#include "../registry.hpp"

namespace atomistica {

static std::set<int> unique_Z(const AtomicSystem& sys) {
    std::set<int> elems;
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        elems.insert(sys.atomic_number(i));
    return elems;
}

// ---------------------------------------------------------------------------
// TersoffWrapper
// ---------------------------------------------------------------------------

class TersoffWrapper : public Potential {
    PotentialWrapper<Tersoff<false>> pot_;
    bool loaded_ = false;
    std::string param_set_;

public:
    explicit TersoffWrapper(const Config& cfg)
        : param_set_(cfg.get_or<std::string>("param_set", ""))
    {
        if (!param_set_.empty()) {
            pot_.get().load_parameters(param_set_);
            loaded_ = true;
        }
    }

    Scalar cutoff() const override { return loaded_ ? pot_.cutoff() : 0.0; }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        if (!loaded_) {
            auto elems = unique_Z(sys);
            if (elems.count(14))
                param_set_ = "Tersoff_PRB_39_5566_Si_C";
            else if (elems.count(13) && elems.count(7))
                param_set_ = "Goumri_Said_ChemPhys_302_135_Al_N";
            else if ((elems.count(5) || elems.count(6) || elems.count(7))
                     && elems.size() <= 3)
                param_set_ = "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N";
            else
                throw std::runtime_error(
                    "Tersoff: cannot auto-select param_set for given elements");
            pot_.get().load_parameters(param_set_);
            loaded_ = true;
        }
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("Tersoff", TersoffWrapper)

// ---------------------------------------------------------------------------
// TersoffScrWrapper
// ---------------------------------------------------------------------------

class TersoffScrWrapper : public Potential {
    PotentialWrapper<Tersoff<true>> pot_;
    bool loaded_ = false;

public:
    explicit TersoffScrWrapper(const Config& cfg) {
        std::string param_set =
            cfg.get_or<std::string>("param_set", "Tersoff_PRB_39_5566_Si_C");
        pot_.get().load_parameters(param_set);
        loaded_ = true;
    }

    Scalar cutoff() const override { return loaded_ ? pot_.cutoff() : 0.0; }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("TersoffScr", TersoffScrWrapper)

// ---------------------------------------------------------------------------
// BrennerWrapper
// ---------------------------------------------------------------------------

class BrennerWrapper : public Potential {
    PotentialWrapper<Brenner<false>> pot_;
    bool loaded_ = false;
    std::string param_set_;

public:
    explicit BrennerWrapper(const Config& cfg)
        : param_set_(cfg.get_or<std::string>("param_set", ""))
    {
        if (!param_set_.empty()) {
            pot_.get().load_parameters(param_set_);
            loaded_ = true;
        }
    }

    Scalar cutoff() const override { return loaded_ ? pot_.cutoff() : 0.0; }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        if (!loaded_) {
            auto elems = unique_Z(sys);
            if (elems.count(14) && elems.count(6))
                param_set_ = "Erhart_PRB_71_035211_SiC";
            else if (elems.count(78) && elems.count(6))
                param_set_ = "Albe_PRB_65_195124_PtC";
            else if (elems.count(26) && elems.count(6))
                param_set_ = "Henriksson_PRB_79_144107_FeC";
            else if (elems.count(6) && elems.size() <= 2)
                param_set_ = "Brenner_PRB_42_9458_C_II";
            else
                throw std::runtime_error(
                    "Brenner: cannot auto-select param_set for given elements");
            pot_.get().load_parameters(param_set_);
            loaded_ = true;
        }
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("Brenner", BrennerWrapper)

// ---------------------------------------------------------------------------
// KumagaiWrapper
// ---------------------------------------------------------------------------

class KumagaiWrapper : public Potential {
    PotentialWrapper<Kumagai<false>> pot_;
    bool loaded_ = false;

public:
    explicit KumagaiWrapper(const Config& /*cfg*/) {
        pot_.get().load_parameters("Kumagai_CompMaterSci_39_457_Si");
        loaded_ = true;
    }

    Scalar cutoff() const override { return loaded_ ? pot_.cutoff() : 0.0; }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("Kumagai", KumagaiWrapper)

// ---------------------------------------------------------------------------
// JuslinWrapper
// ---------------------------------------------------------------------------

class JuslinWrapper : public Potential {
    PotentialWrapper<Juslin<false>> pot_;
    bool loaded_ = false;

public:
    explicit JuslinWrapper(const Config& /*cfg*/) {
        pot_.get().load_parameters("Juslin_JAP_98_123520_WCH");
        loaded_ = true;
    }

    Scalar cutoff() const override { return loaded_ ? pot_.cutoff() : 0.0; }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("Juslin", JuslinWrapper)

// ---------------------------------------------------------------------------
// REBO2Wrapper
// ---------------------------------------------------------------------------

class REBO2Wrapper : public Potential {
    PotentialWrapper<REBO2> pot_;
    bool loaded_ = false;

public:
    explicit REBO2Wrapper(const Config& /*cfg*/) {
        pot_.get().load_default_parameters();
        loaded_ = true;
    }

    Scalar cutoff() const override { return loaded_ ? pot_.cutoff() : 0.0; }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("REBO2", REBO2Wrapper)
// "Rebo2" case-insensitive alias handled by Registry::to_lower in make_potential

} // namespace atomistica
