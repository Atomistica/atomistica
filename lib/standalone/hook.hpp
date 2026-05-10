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

// The three points in the velocity-Verlet loop where hooks are called:
//
//   PRE_FORCE  — after position update, before force evaluation
//                Use for: position constraints (SETTLE), cell modifications
//                         (ConstantStrainRate), coordinate resets
//   POST_FORCE — after force evaluation, before velocity half-step
//                Use for: force constraints, force logging
//   POST_STEP  — after velocity half-step (end of full MD step)
//                Use for: thermostats, barostats, output, analysis
//
// Within each point, hooks are invoked in ascending priority order.
// Suggested priority bands:
//   PRE_FORCE  0-9   : cell modifiers (ConstantStrainRate, Lees-Edwards)
//              10-19 : position/force constraints (SETTLE, ConstantForce)
//              5     : charge equilibration
//   POST_STEP  20-29 : thermostats (BerendsenT, LangevinT, PetersT)
//              30-39 : barostats (BerendsenP)
//              40-49 : analysis (DiffusionCoefficient, Slicing, HeatFlux)
//              50-59 : output (OutputEnergy, OutputXYZ, OutputCFG, OutputNC)
enum class HookPoint { PRE_FORCE = 0, POST_FORCE = 1, POST_STEP = 2 };

class Hook {
public:
    virtual ~Hook() = default;

    // Called at the hook's registered point each MD step.
    virtual void invoke(SimulationContext& ctx) = 0;

    // Which point in the step this hook fires at. Default: POST_STEP.
    virtual HookPoint hook_point() const { return HookPoint::POST_STEP; }

    // Within a given hook point, lower priority fires first. Default: 50.
    virtual int priority() const { return 50; }

    // Called once after system and neighbor list are fully initialised.
    // Override to cache atom indices, allocate per-atom storage, etc.
    virtual void bind_to(AtomicSystem& /*system*/, NeighborList& /*nl*/) {}

    // Human-readable name for logging/error messages.
    virtual std::string name() const { return "Hook"; }
};

} // namespace atomistica
