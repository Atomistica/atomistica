// ======================================================================
// Atomistica — GPL-2.0-or-later
// Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
// ======================================================================

#pragma once

#include <stdexcept>
#include <string>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../../include/atomistica/potentials/potential_base.hpp"
#include "../../include/atomistica/potentials/eam/eam.hpp"

#include "../config.hpp"
#include "../registry.hpp"

namespace atomistica {

class TabulatedEAMWrapper : public Potential {
    PotentialWrapper<TabulatedEAM> pot_;

public:
    explicit TabulatedEAMWrapper(const Config& cfg) {
        std::string file = cfg.get_or<std::string>("file", "");
        if (file.empty())
            throw std::runtime_error("TabulatedEAM: 'file' required");
        pot_.get().load(file);
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

REGISTER_POTENTIAL("TabulatedEAM", TabulatedEAMWrapper)

class TabulatedAlloyEAMWrapper : public Potential {
    PotentialWrapper<TabulatedAlloyEAM> pot_;

public:
    explicit TabulatedAlloyEAMWrapper(const Config& cfg) {
        std::string file = cfg.get_or<std::string>("file", "");
        if (file.empty())
            throw std::runtime_error("TabulatedAlloyEAM: 'file' required");
        pot_.get().load(file);
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

REGISTER_POTENTIAL("TabulatedAlloyEAM", TabulatedAlloyEAMWrapper)

} // namespace atomistica
