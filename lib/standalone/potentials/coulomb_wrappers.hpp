// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <stdexcept>
#include <vector>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../../include/atomistica/potentials/potential_base.hpp"
#include "../../include/atomistica/potentials/coulomb/coulomb.hpp"

#include "../config.hpp"
#include "../registry.hpp"
#include "../coulomb_solver.hpp"

namespace atomistica {

// ---------------------------------------------------------------------------
// DirectCoulombWrapper
// ---------------------------------------------------------------------------

class DirectCoulombWrapper : public CoulombSolver {
    DirectCoulomb solver_;

public:
    explicit DirectCoulombWrapper(const Config& cfg) {
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_.set_epsilon_r(eps_r);
    }

    Scalar cutoff() const override { return solver_.cutoff(); }

    void bind_to(AtomicSystem& /*sys*/, NeighborList& /*nl*/) override {}

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            std::vector<Scalar> charges(q.data(), q.data() + q.size());
            solver_.set_charges(charges);
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    void compute_potential(AtomicSystem& /*sys*/, NeighborList& /*nl*/,
                           const ArrayX& /*q*/, ArrayX& /*phi*/) override {
        throw std::runtime_error("DirectCoulomb::compute_potential: not implemented");
    }
};

REGISTER_COULOMB("DirectCoulomb", DirectCoulombWrapper)

// ---------------------------------------------------------------------------
// CutoffCoulombWrapper
// ---------------------------------------------------------------------------

class CutoffCoulombWrapper : public CoulombSolver {
    CutoffCoulomb solver_;

public:
    explicit CutoffCoulombWrapper(const Config& cfg) {
        double rc    = cfg.get_or<double>("cutoff", 10.0);
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_.set_cutoff(rc);
        solver_.set_epsilon_r(eps_r);
    }

    Scalar cutoff() const override { return solver_.cutoff(); }

    void bind_to(AtomicSystem& /*sys*/, NeighborList& /*nl*/) override {}

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            std::vector<Scalar> charges(q.data(), q.data() + q.size());
            solver_.set_charges(charges);
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    void compute_potential(AtomicSystem& /*sys*/, NeighborList& /*nl*/,
                           const ArrayX& /*q*/, ArrayX& /*phi*/) override {
        throw std::runtime_error("CutoffCoulomb::compute_potential: not implemented");
    }
};

REGISTER_COULOMB("CutoffCoulomb", CutoffCoulombWrapper)

// ---------------------------------------------------------------------------
// WolfCoulombWrapper
// ---------------------------------------------------------------------------

class WolfCoulombWrapper : public CoulombSolver {
    WolfCoulomb solver_;

public:
    explicit WolfCoulombWrapper(const Config& cfg) {
        double rc    = cfg.get_or<double>("cutoff", 10.0);
        double alpha = cfg.get_or<double>("alpha", 0.2);
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_.set_cutoff(rc);
        solver_.set_alpha(alpha);
        solver_.set_epsilon_r(eps_r);
    }

    Scalar cutoff() const override { return solver_.cutoff(); }

    void bind_to(AtomicSystem& /*sys*/, NeighborList& /*nl*/) override {}

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            std::vector<Scalar> charges(q.data(), q.data() + q.size());
            solver_.set_charges(charges);
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    void compute_potential(AtomicSystem& /*sys*/, NeighborList& /*nl*/,
                           const ArrayX& /*q*/, ArrayX& /*phi*/) override {
        throw std::runtime_error("WolfCoulomb::compute_potential: not implemented");
    }
};

REGISTER_COULOMB("WolfCoulomb", WolfCoulombWrapper)

} // namespace atomistica
