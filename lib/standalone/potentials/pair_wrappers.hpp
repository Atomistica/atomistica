// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../../include/atomistica/potentials/potential_base.hpp"
#include "../../include/atomistica/potentials/pair/lj.hpp"
#include "../../include/atomistica/potentials/pair/simple_pairs.hpp"

#include "../config.hpp"
#include "../registry.hpp"

namespace atomistica {

// ---------------------------------------------------------------------------
// LJCutWrapper
// ---------------------------------------------------------------------------

class LJCutWrapper : public Potential {
    PotentialWrapper<LJCut> pot_;

public:
    explicit LJCutWrapper(const Config& cfg) {
        int Z1       = cfg.get_or<int>("Z1", 0);
        int Z2       = cfg.get_or<int>("Z2", 0);
        double eps   = cfg.get_or<double>("epsilon", 1.0);
        double sigma = cfg.get_or<double>("sigma", 1.0);
        double rc    = cfg.get_or<double>("cutoff", 5.0);
        pot_.get().set_params(Z1, Z2, eps, sigma, rc);
    }

    Scalar cutoff() const override { return pot_.cutoff(); }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("LJCut", LJCutWrapper)

// ---------------------------------------------------------------------------
// LJCutShiftWrapper
// ---------------------------------------------------------------------------

class LJCutShiftWrapper : public Potential {
    PotentialWrapper<LJCutShift> pot_;

public:
    explicit LJCutShiftWrapper(const Config& cfg) {
        int Z1       = cfg.get_or<int>("Z1", 0);
        int Z2       = cfg.get_or<int>("Z2", 0);
        double eps   = cfg.get_or<double>("epsilon", 1.0);
        double sigma = cfg.get_or<double>("sigma", 1.0);
        double rc    = cfg.get_or<double>("cutoff", 5.0);
        pot_.get().set_params(Z1, Z2, eps, sigma, rc);
    }

    Scalar cutoff() const override { return pot_.cutoff(); }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("LJCutShift", LJCutShiftWrapper)

// ---------------------------------------------------------------------------
// BornMayerWrapper
// ---------------------------------------------------------------------------

class BornMayerWrapper : public Potential {
    PotentialWrapper<BornMayer> pot_;

public:
    explicit BornMayerWrapper(const Config& cfg) {
        int Z1     = cfg.get_or<int>("Z1", 0);
        int Z2     = cfg.get_or<int>("Z2", 0);
        double A   = cfg.get_or<double>("A", 1000.0);
        double rho = cfg.get_or<double>("rho", 0.3);
        double rc  = cfg.get_or<double>("cutoff", 5.0);
        pot_.get().Z1             = Z1;
        pot_.get().Z2             = Z2;
        pot_.get().A              = A;
        pot_.get().rho            = rho;
        pot_.get().cutoff_radius  = rc;
    }

    Scalar cutoff() const override { return pot_.cutoff(); }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("BornMayer", BornMayerWrapper)

// ---------------------------------------------------------------------------
// HarmonicPairWrapper
// ---------------------------------------------------------------------------

class HarmonicPairWrapper : public Potential {
    PotentialWrapper<Harmonic> pot_;

public:
    explicit HarmonicPairWrapper(const Config& cfg) {
        int Z1    = cfg.get_or<int>("Z1", 0);
        int Z2    = cfg.get_or<int>("Z2", 0);
        double k  = cfg.get_or<double>("k", 1.0);
        double r0 = cfg.get_or<double>("r0", 1.0);
        double rc = cfg.get_or<double>("cutoff", 3.0);
        bool shift = cfg.get_or<bool>("shift", false);
        pot_.get().Z1            = Z1;
        pot_.get().Z2            = Z2;
        pot_.get().k             = k;
        pot_.get().r0            = r0;
        pot_.get().cutoff_radius = rc;
        pot_.get().shift         = shift;
    }

    Scalar cutoff() const override { return pot_.cutoff(); }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("HarmonicPair", HarmonicPairWrapper)

// ---------------------------------------------------------------------------
// R6Wrapper
// ---------------------------------------------------------------------------

class R6Wrapper : public Potential {
    PotentialWrapper<R6> pot_;

public:
    explicit R6Wrapper(const Config& cfg) {
        int Z1    = cfg.get_or<int>("Z1", 0);
        int Z2    = cfg.get_or<int>("Z2", 0);
        double A  = cfg.get_or<double>("A", 1.0);
        double r0 = cfg.get_or<double>("r0", 1.0);
        double rc = cfg.get_or<double>("cutoff", 5.0);
        pot_.get().Z1            = Z1;
        pot_.get().Z2            = Z2;
        pot_.get().A             = A;
        pot_.get().r0            = r0;
        pot_.get().cutoff_radius = rc;
    }

    Scalar cutoff() const override { return pot_.cutoff(); }

    void bind_to(AtomicSystem& sys, NeighborList& nl) override {
        pot_.bind_to(sys, nl);
    }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                             bool cf = true, bool cv = true) override {
        return pot_.compute(sys, nl, cf, cv);
    }
};

REGISTER_POTENTIAL("R6", R6Wrapper)

} // namespace atomistica
