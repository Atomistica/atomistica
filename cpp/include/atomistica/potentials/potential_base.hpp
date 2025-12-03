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

#include <optional>

#include "../config.hpp"
#include "../core/atomic_system.hpp"
#include "../core/neighbor_list.hpp"

namespace atomistica {

/**
 * @brief Results from potential energy/force computation
 */
struct PotentialResults {
    Scalar energy = 0.0;           // Total potential energy
    Mat3 virial = Mat3::Zero();    // Virial stress tensor

    // Optional per-atom decomposition
    std::optional<ArrayX> energy_per_atom;
};

/**
 * @brief CRTP base class for interatomic potentials
 *
 * Provides common interface and shared functionality.
 * Derived classes implement the actual physics.
 *
 * @tparam Derived The derived potential class (CRTP pattern)
 */
template<typename Derived>
class PotentialBase {
public:
    /**
     * @brief Get the interaction cutoff
     */
    Scalar cutoff() const {
        return static_cast<const Derived*>(this)->cutoff_impl();
    }

    /**
     * @brief Bind potential to a particle system
     *
     * Called before the first energy/force computation.
     * Can be used to allocate per-atom storage, set up element-specific parameters, etc.
     */
    void bind_to(AtomicSystem& system, NeighborList& neighbors) {
        static_cast<Derived*>(this)->bind_to_impl(system, neighbors);
    }

    /**
     * @brief Compute energy and forces
     *
     * Forces are accumulated into system.forces().
     *
     * @param system Atomic system (positions, cell, etc.)
     * @param neighbors Neighbor list
     * @param compute_forces Whether to compute forces
     * @param compute_virial Whether to compute virial stress
     * @return PotentialResults containing energy and optionally virial
     */
    PotentialResults compute(AtomicSystem& system,
                            NeighborList& neighbors,
                            bool compute_forces = true,
                            bool compute_virial = true) {
        return static_cast<Derived*>(this)->compute_impl(
            system, neighbors, compute_forces, compute_virial);
    }

protected:
    // Default implementations (can be overridden by derived classes)
    void bind_to_impl(AtomicSystem& /*system*/, NeighborList& /*neighbors*/) {
        // Default: do nothing
    }
};

/**
 * @brief Abstract base class for runtime polymorphism
 *
 * Use this when you need to store potentials in containers or
 * switch between potentials at runtime.
 */
class Potential {
public:
    virtual ~Potential() = default;

    virtual Scalar cutoff() const = 0;
    virtual void bind_to(AtomicSystem& system, NeighborList& neighbors) = 0;
    virtual PotentialResults compute(AtomicSystem& system,
                                     NeighborList& neighbors,
                                     bool compute_forces = true,
                                     bool compute_virial = true) = 0;
};

/**
 * @brief Wrapper to use CRTP potentials with virtual interface
 */
template<typename CRTPPotential>
class PotentialWrapper : public Potential {
public:
    template<typename... Args>
    explicit PotentialWrapper(Args&&... args)
        : potential_(std::forward<Args>(args)...) {}

    Scalar cutoff() const override {
        return potential_.cutoff();
    }

    void bind_to(AtomicSystem& system, NeighborList& neighbors) override {
        potential_.bind_to(system, neighbors);
    }

    PotentialResults compute(AtomicSystem& system,
                            NeighborList& neighbors,
                            bool compute_forces = true,
                            bool compute_virial = true) override {
        return potential_.compute(system, neighbors, compute_forces, compute_virial);
    }

    CRTPPotential& get() { return potential_; }
    const CRTPPotential& get() const { return potential_; }

private:
    CRTPPotential potential_;
};

} // namespace atomistica
