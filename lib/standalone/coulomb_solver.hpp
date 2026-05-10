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

#include "../include/atomistica/potentials/potential_base.hpp"

namespace atomistica {

// Abstract base for electrostatic (Coulomb) solvers.
//
// Extends Potential so it participates in the normal force-evaluation loop
// inside SimulationContext::compute_forces(). Adds compute_potential() for
// use in charge self-consistency (SCF) loops inside the ChargeEquilibration
// hook (PRE_FORCE, priority 5).
//
// Lifecycle:
//   1. Registered in the Registry as a Coulomb type ("DirectCoulomb", "PME", …)
//   2. Stored in SimulationContext::coulomb (unique ownership, separate from
//      the potentials vector)
//   3. bind_to() called once during setup
//   4. Each step: ChargeEquilibration invokes compute_potential() iteratively
//      until charges converge, writing them to
//      system.properties().get<ArrayX>("charges")
//   5. SimulationContext::compute_forces() calls compute() which reads charges
//      from properties and adds electrostatic E/F/W to the accumulated results
//
// Implementors must override Potential::compute() to read charges from
//   system.properties().get<ArrayX>("charges")
// and accumulate forces into system.forces().
class CoulombSolver : public Potential {
public:
    // Compute the electrostatic potential phi_i at each atom site given an
    // explicit charge array q.  Used by ChargeEquilibration during the SCF
    // loop; charges are not yet committed to system.properties() at this point.
    //
    // q   — per-atom charges, size == system.num_atoms()
    // phi — output, per-atom electrostatic potential, same size; overwritten
    virtual void compute_potential(AtomicSystem& system,
                                   NeighborList& nl,
                                   const ArrayX& q,
                                   ArrayX& phi) = 0;
};

} // namespace atomistica
