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

#include <cmath>

#include "../include/atomistica/config.hpp"
#include "../include/atomistica/core/atomic_system.hpp"

namespace atomistica {

// Boltzmann constant in eV/K
constexpr double kB_eV = 8.617333262e-5;

// 1 amu*(Å/fs)^2 = AMU_AFSQ_PER_EV eV  →  sqrt of this is the eV/A/fs time scale
constexpr double AMU_AFSQ_PER_EV_CONST = 103.636;

// sqrt(AMU_AFSQ_PER_EV): multiply internal time by this to get femtoseconds
constexpr double SQRT_AMU_AFSQ_PER_EV = 10.180;  // sqrt(103.636)

// -----------------------------------------------------------------------
// Group / freeze helpers
// -----------------------------------------------------------------------

// Returns true if atom i is frozen (group <= 0).
// If no "group" property exists, all atoms are considered mobile.
inline bool is_frozen(const AtomicSystem& sys, size_t i) {
    if (!sys.properties().has("group")) return false;
    return sys.properties().get<ArrayXi>("group")[static_cast<int>(i)] <= 0;
}

// -----------------------------------------------------------------------
// Kinematic helpers
// -----------------------------------------------------------------------

// Total kinetic energy (eV), optionally skipping frozen atoms.
inline double kinetic_energy(const AtomicSystem& sys, bool skip_frozen = false) {
    double ekin = 0.0;
    size_t n = sys.num_atoms();
    for (size_t i = 0; i < n; ++i) {
        if (skip_frozen && is_frozen(sys, i)) continue;
        double m = sys.mass(i);
        Vec3 v = sys.velocity(i);
        ekin += 0.5 * m * v.squaredNorm();
    }
    return ekin;
}

// Number of degrees of freedom (3N - 3, subtracting COM motion).
inline int n_dof(const AtomicSystem& sys, bool skip_frozen = false) {
    int n_mobile = 0;
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        if (!skip_frozen || !is_frozen(sys, i)) ++n_mobile;
    return std::max(1, 3 * n_mobile - 3);
}

// Instantaneous temperature (K) from kinetic energy.
inline double temperature_K(const AtomicSystem& sys, bool skip_frozen = false) {
    double ekin = kinetic_energy(sys, skip_frozen);
    int dof = n_dof(sys, skip_frozen);
    return 2.0 * ekin / (dof * kB_eV);
}

// Kinetic contribution to the virial tensor: W_kin = sum_i m_i v_i ⊗ v_i
inline Mat3 kinetic_virial(const AtomicSystem& sys) {
    Mat3 W = Mat3::Zero();
    for (size_t i = 0; i < sys.num_atoms(); ++i) {
        double m = sys.mass(i);
        Vec3 v = sys.velocity(i);
        W += m * v * v.transpose();
    }
    return W;
}

// Scalar pressure (eV/Å³) from potential virial and kinetic virial.
// virial = ctx.results.virial (convention: stored as -W_pot, so pressure = +virial + W_kin)
inline double pressure(const Mat3& pot_virial, const AtomicSystem& sys) {
    double vol = sys.volume();
    if (vol <= 0.0) return 0.0;
    Mat3 W_kin = kinetic_virial(sys);
    return (pot_virial + W_kin).trace() / (3.0 * vol);
}

// Maximum force magnitude (eV/Å).
inline double fmax(const AtomicSystem& sys) {
    double fm2 = 0.0;
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        fm2 = std::max(fm2, sys.forces().col(i).matrix().squaredNorm());
    return std::sqrt(fm2);
}

// Wrap all positions into the primary cell.
inline void wrap_positions(AtomicSystem& sys) {
    for (size_t i = 0; i < sys.num_atoms(); ++i)
        sys.set_position(i, sys.wrap_position(sys.position(i)));
}

} // namespace atomistica
