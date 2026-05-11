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

#include <string>

#include "../include/atomistica/core/atomic_system.hpp"
#include "../include/atomistica/core/neighbor_list.hpp"

namespace atomistica {

// Forward declaration — full definition in simulation_context.hpp
struct SimulationContext;

// Abstract base for time integrators.
//
// An integrator owns the inner structure of a single MD step. It is
// responsible for firing hook points at the right moments via
// ctx.invoke_hooks() and calling ctx.compute_forces() at the appropriate
// time. This lets each integrator (VelocityVerlet, FIRE, …) define its own
// algorithm without exposing those details to the rest of the code.
//
// Typical velocity-Verlet implementation:
//   1. half_kick(ctx)          — v += 0.5*dt*f/m
//   2. drift(ctx)              — r += dt*v; wrap positions
//   3. ctx.invoke_hooks(PRE_FORCE)
//   4. ctx.nl.update(ctx.system)
//   5. ctx.compute_forces()
//   6. ctx.invoke_hooks(POST_FORCE)
//   7. half_kick(ctx)          — v += 0.5*dt*f/m
//   8. ctx.invoke_hooks(POST_STEP)
class Integrator {
public:
    virtual ~Integrator() = default;

    // Perform one integration step.
    // Returns true to continue, false to stop (e.g., FIRE has converged).
    virtual bool step(SimulationContext& ctx) = 0;

    // Called once after system and neighbor list are fully initialised.
    virtual void bind_to(AtomicSystem& /*system*/, NeighborList& /*nl*/) {}

    // Human-readable name for logging/error messages.
    virtual std::string name() const { return "Integrator"; }
};

} // namespace atomistica
