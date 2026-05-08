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
#include <random>

#include "../include/atomistica/config.hpp"
#include "../include/atomistica/core/atomic_system.hpp"
#include "../include/atomistica/core/neighbor_list.hpp"
#include "../include/atomistica/integrators/thermostats.hpp"  // for kB_eV_K

namespace atomistica {

// Peters (DPD) thermostat.
// Reference: E.A.J.F. Peters, Europhys. Lett. 66, 311 (2004).
//
// Works in the C++ natural unit system (eV_A):
//   - masses in amu
//   - velocities in sqrt(eV/amu)
//   - dt in sqrt(amu*Å²/eV)
//   - gamma in amu / sqrt(amu*Å²/eV) = sqrt(amu*eV)/Å  (makes weight dimensionless)
//
// For eV_A_fs input: convert gamma → gamma / sqrt(103.636) before constructing.
class PetersT {
public:
    PetersT(double T, double gamma, double cutoff, unsigned int seed = 12345)
        : T_(T), gamma_(gamma), cutoff_(cutoff), rng_(seed) {}

    void set_temperature(double T) { T_ = T; }
    double temperature() const { return T_; }

    // Apply one thermostat step.
    // dt is in C++ natural time units (sqrt(amu*Å²/eV)).
    void apply(AtomicSystem& system, const NeighborList& nl, double dt) {
        std::normal_distribution<double> normal(0.0, 1.0);
        size_t nat = system.num_atoms();

        for (size_t i = 0; i < nat; ++i) {
            auto [begin, end] = nl.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                size_t j = it->index;
                if (j <= i) continue;  // each pair once

                // Displacement vector with PBC image
                Vec3 shift;
                shift << static_cast<Scalar>(it->cell_shift[0]),
                         static_cast<Scalar>(it->cell_shift[1]),
                         static_cast<Scalar>(it->cell_shift[2]);
                Vec3 dr = system.position(j) + system.cell() * shift
                          - system.position(i);

                Scalar r2 = dr.squaredNorm();
                if (r2 >= cutoff_ * cutoff_) continue;

                Scalar r = std::sqrt(r2);
                Vec3 rhat = dr / r;

                Scalar m_i = system.mass(i);
                Scalar m_j = system.mass(j);
                Scalar r_muij = 1.0 / m_i + 1.0 / m_j;
                Scalar mu = 1.0 / r_muij;

                // Linear decay kernel: w(r) = 1 - r/cutoff
                Scalar kernel = 1.0 - r / cutoff_;

                Scalar w = static_cast<Scalar>(dt) * r_muij
                           * static_cast<Scalar>(gamma_) * kernel;

                Scalar a = mu * (1.0 - std::exp(-w));
                Scalar b_var = kB_eV_K * static_cast<Scalar>(T_) * mu
                               * (1.0 - std::exp(-2.0 * w));
                Scalar b = (b_var > 0.0 ? std::sqrt(b_var) : 0.0)
                           * static_cast<Scalar>(normal(rng_));

                Vec3 vij = system.velocity(i) - system.velocity(j);
                Vec3 dmom = (-a * vij.dot(rhat) + b) * rhat;

                system.set_velocity(i, system.velocity(i) + dmom / m_i);
                system.set_velocity(j, system.velocity(j) - dmom / m_j);
            }
        }
    }

private:
    double T_;
    double gamma_;
    double cutoff_;
    std::mt19937 rng_;
};

} // namespace atomistica
