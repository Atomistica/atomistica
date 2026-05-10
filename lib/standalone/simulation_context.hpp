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

#pragma once

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <vector>

#include "../include/atomistica/core/atomic_system.hpp"
#include "../include/atomistica/core/neighbor_list.hpp"
#include "../include/atomistica/potentials/potential_base.hpp"

#include "coulomb_solver.hpp"
#include "hook.hpp"

namespace atomistica {

// Central hub for one MD simulation.
//
// SimulationContext is the single object passed through every step,
// hook invocation, and force calculation. It owns all simulation
// components (potentials, Coulomb solver, hooks) but holds non-owning
// references to AtomicSystem and NeighborList, which are owned by the
// caller (mdcore.cpp main()).
//
// Typical setup sequence in mdcore.cpp:
//   ctx.potentials.push_back(registry.make_potential("REBO2", cfg));
//   ctx.coulomb = registry.make_coulomb("DirectCoulomb", cfg);
//   ctx.hooks.push_back(registry.make_hook("BerendsenT", cfg));
//   ctx.sort_hooks();
//   ctx.nl.set_cutoff(ctx.max_cutoff() + cutoff_add);
//   ctx.nl.update(ctx.system);
//   ctx.bind_all();
//   ctx.compute_forces();          // initial forces
//
//   while (integrator->step(ctx)) { }
struct SimulationContext {
    // Non-owning references — owned by the caller
    AtomicSystem& system;
    NeighborList& nl;

    // Owned simulation components (populated before the loop starts)
    std::vector<std::unique_ptr<Potential>> potentials;
    std::unique_ptr<CoulombSolver>          coulomb;   // nullptr if not used

    // Hooks sorted by (hook_point, priority) — call sort_hooks() after all
    // hooks are added.
    std::vector<std::unique_ptr<Hook>>      hooks;

    // Accumulated energy and virial from the last compute_forces() call.
    PotentialResults results;

    // Time-stepping state
    double dt   = 0.0;
    double time = 0.0;
    int    step = 0;

    explicit SimulationContext(AtomicSystem& sys, NeighborList& nl_in)
        : system(sys), nl(nl_in) {}

    // Non-copyable, non-movable (holds references)
    SimulationContext(const SimulationContext&)            = delete;
    SimulationContext& operator=(const SimulationContext&) = delete;

    // -----------------------------------------------------------------------
    // Setup helpers
    // -----------------------------------------------------------------------

    // Sort hooks by (hook_point ascending, priority ascending).
    // Must be called once after all hooks have been pushed.
    void sort_hooks() {
        std::stable_sort(hooks.begin(), hooks.end(),
            [](const std::unique_ptr<Hook>& a, const std::unique_ptr<Hook>& b) {
                if (a->hook_point() != b->hook_point())
                    return static_cast<int>(a->hook_point())
                         < static_cast<int>(b->hook_point());
                return a->priority() < b->priority();
            });
    }

    // Return the maximum cutoff across all potentials and the Coulomb solver.
    // Use this to set nl.set_cutoff() before bind_all().
    double max_cutoff() const {
        double c = 0.0;
        for (const auto& p : potentials)
            c = std::max(c, static_cast<double>(p->cutoff()));
        if (coulomb)
            c = std::max(c, static_cast<double>(coulomb->cutoff()));
        return c;
    }

    // Call bind_to() on every potential, Coulomb solver, and hook.
    // The neighbor list must already have its cutoff set and be up to date
    // before this is called.
    void bind_all() {
        for (auto& p : potentials)
            p->bind_to(system, nl);
        if (coulomb)
            coulomb->bind_to(system, nl);
        for (auto& h : hooks)
            h->bind_to(system, nl);
    }

    // -----------------------------------------------------------------------
    // Per-step operations (called from Integrator::step)
    // -----------------------------------------------------------------------

    // Zero forces, then compute and accumulate E/F/W from every potential
    // and (if present) the Coulomb solver.
    //
    // Forces are accumulated into system.forces() by each potential's
    // compute() call. Energy and virial are accumulated in results.
    //
    // If a Coulomb solver is present, it reads charges from
    //   system.properties().get<ArrayX>("charges")
    // which must have been set by the ChargeEquilibration hook (PRE_FORCE)
    // before this call.
    void compute_forces() {
        results.energy = 0.0;
        results.virial = Mat3::Zero();
        system.zero_forces();

        for (auto& p : potentials) {
            PotentialResults r = p->compute(system, nl);
            results.energy += r.energy;
            results.virial += r.virial;
        }

        if (coulomb) {
            PotentialResults r = coulomb->compute(system, nl);
            results.energy += r.energy;
            results.virial += r.virial;
        }
    }

    // Invoke all hooks registered at the given point, in priority order.
    // Hooks are assumed to already be sorted (call sort_hooks() at setup).
    void invoke_hooks(HookPoint point) {
        for (auto& h : hooks) {
            if (h->hook_point() == point)
                h->invoke(*this);
        }
    }
};

} // namespace atomistica
