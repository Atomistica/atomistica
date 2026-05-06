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
#include <stdexcept>

#include "../../config.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../potential_base.hpp"

namespace atomistica {

// ============================================================================
// Helper: element pair matching (Z=0 means wildcard = any element)
// ============================================================================

inline bool pair_matches(int Za, int Zb, int Z1, int Z2) {
    bool fwd = (Z1 == 0 || Za == Z1) && (Z2 == 0 || Zb == Z2);
    bool rev = (Z1 == 0 || Zb == Z1) && (Z2 == 0 || Za == Z2);
    return fwd || rev;
}

// ============================================================================
// BornMayer potential
// ============================================================================

/**
 * @brief Born-Mayer repulsive potential
 *
 * V(r) = A * exp(-r/rho) - shift
 * where shift = A * exp(-cutoff/rho) ensures V(cutoff) = 0.
 *
 * Reference: Born & Mayer (1932)
 */
class BornMayer : public PotentialBase<BornMayer> {
public:
    Scalar A = 1.0;
    Scalar rho = 1.0;
    Scalar cutoff_radius = 1.0;
    int Z1 = 0;  // 0 = wildcard
    int Z2 = 0;

    BornMayer() = default;

    BornMayer(Scalar A_, Scalar rho_, Scalar cutoff_, int Z1_ = 0, int Z2_ = 0)
        : A(A_), rho(rho_), cutoff_radius(cutoff_), Z1(Z1_), Z2(Z2_) {}

    Scalar cutoff() const { return cutoff_radius; }

    void bind_to(AtomicSystem&, NeighborList& nl) {
        nl.set_cutoff(std::max(nl.cutoff(), cutoff_radius));
    }

    PotentialResults compute(AtomicSystem& system, NeighborList& neighbors,
                             bool compute_forces = true,
                             bool compute_virial = true)
    {
        PotentialResults results;
        const Scalar shift = A * std::exp(-cutoff_radius / rho);
        const Mat3& cell = system.cell();

        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            int Zi = system.atomic_numbers()(i);
            auto [begin, end] = neighbors.neighbors(i);

            for (auto it = begin; it != end; ++it) {
                std::size_t j = it->index;
                if (j <= i) continue;

                int Zj = system.atomic_numbers()(j);
                if (!pair_matches(Zi, Zj, Z1, Z2)) continue;

                Vec3 dr = system.position(j) - system.position(i);
                dr += cell.col(0) * it->cell_shift[0];
                dr += cell.col(1) * it->cell_shift[1];
                dr += cell.col(2) * it->cell_shift[2];

                Scalar r = dr.norm();
                if (r >= cutoff_radius || r < 1e-10) continue;

                Scalar exp_r = A * std::exp(-r / rho);
                Scalar V = exp_r - shift;
                results.energy += V;

                if (compute_forces || compute_virial) {
                    Scalar dV_dr = -(A / rho) * std::exp(-r / rho);
                    Vec3 f = -(dV_dr / r) * dr;

                    if (compute_forces) {
                        system.forces().col(i) -= f.array();
                        system.forces().col(j) += f.array();
                    }

                    if (compute_virial) {
                        results.virial += dr * f.transpose();
                    }
                }
            }
        }

        return results;
    }
};

// ============================================================================
// Harmonic potential
// ============================================================================

/**
 * @brief Harmonic spring potential
 *
 * V(r) = 0.5 * k * (r - r0)^2 - offset
 * If shift=true: offset = 0.5 * k * (cutoff - r0)^2 (so V(cutoff)=0)
 * If shift=false: offset = 0
 */
class Harmonic : public PotentialBase<Harmonic> {
public:
    Scalar k = 1.0;
    Scalar r0 = 1.0;
    Scalar cutoff_radius = 1.5;
    bool shift = false;
    int Z1 = 0;
    int Z2 = 0;

    Harmonic() = default;

    Harmonic(Scalar k_, Scalar r0_, Scalar cutoff_, bool shift_ = false,
             int Z1_ = 0, int Z2_ = 0)
        : k(k_), r0(r0_), cutoff_radius(cutoff_), shift(shift_), Z1(Z1_), Z2(Z2_) {}

    Scalar cutoff() const { return cutoff_radius; }

    void bind_to(AtomicSystem&, NeighborList& nl) {
        nl.set_cutoff(std::max(nl.cutoff(), cutoff_radius));
    }

    PotentialResults compute(AtomicSystem& system, NeighborList& neighbors,
                             bool compute_forces = true,
                             bool compute_virial = true)
    {
        PotentialResults results;
        Scalar offset = shift ? 0.5 * k * (cutoff_radius - r0) * (cutoff_radius - r0) : 0.0;
        const Mat3& cell = system.cell();

        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            int Zi = system.atomic_numbers()(i);
            auto [begin, end] = neighbors.neighbors(i);

            for (auto it = begin; it != end; ++it) {
                std::size_t j = it->index;
                if (j <= i) continue;

                int Zj = system.atomic_numbers()(j);
                if (!pair_matches(Zi, Zj, Z1, Z2)) continue;

                Vec3 dr = system.position(j) - system.position(i);
                dr += cell.col(0) * it->cell_shift[0];
                dr += cell.col(1) * it->cell_shift[1];
                dr += cell.col(2) * it->cell_shift[2];

                Scalar r = dr.norm();
                if (r >= cutoff_radius || r < 1e-10) continue;

                Scalar dr_r0 = r - r0;
                Scalar V = 0.5 * k * dr_r0 * dr_r0 - offset;
                results.energy += V;

                if (compute_forces || compute_virial) {
                    Scalar dV_dr = k * dr_r0;
                    Vec3 f = -(dV_dr / r) * dr;

                    if (compute_forces) {
                        system.forces().col(i) -= f.array();
                        system.forces().col(j) += f.array();
                    }

                    if (compute_virial) {
                        results.virial += dr * f.transpose();
                    }
                }
            }
        }

        return results;
    }
};

// ============================================================================
// DoubleHarmonic potential
// ============================================================================

/**
 * @brief Double-harmonic potential (two-well spring)
 *
 * For r < rm = (r1 + r2) / 2:  V(r) = 0.5 * k1 * (r - r1)^2
 * For r >= rm:                  V(r) = 0.5 * k2 * (r - r2)^2
 *
 * Typically r1 < r2 so there is a well at r1 and another at r2.
 */
class DoubleHarmonic : public PotentialBase<DoubleHarmonic> {
public:
    Scalar k1 = 1.0;
    Scalar r1 = 1.0;
    Scalar k2 = 1.0;
    Scalar r2 = 1.2;
    Scalar cutoff_radius = 1.5;
    int Z1 = 0;
    int Z2 = 0;

    DoubleHarmonic() = default;

    DoubleHarmonic(Scalar k1_, Scalar r1_, Scalar k2_, Scalar r2_,
                   Scalar cutoff_, int Z1_ = 0, int Z2_ = 0)
        : k1(k1_), r1(r1_), k2(k2_), r2(r2_), cutoff_radius(cutoff_),
          Z1(Z1_), Z2(Z2_) {}

    Scalar cutoff() const { return cutoff_radius; }

    void bind_to(AtomicSystem&, NeighborList& nl) {
        nl.set_cutoff(std::max(nl.cutoff(), cutoff_radius));
    }

    PotentialResults compute(AtomicSystem& system, NeighborList& neighbors,
                             bool compute_forces = true,
                             bool compute_virial = true)
    {
        PotentialResults results;
        Scalar rm = 0.5 * (r1 + r2);
        const Mat3& cell = system.cell();

        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            int Zi = system.atomic_numbers()(i);
            auto [begin, end] = neighbors.neighbors(i);

            for (auto it = begin; it != end; ++it) {
                std::size_t j = it->index;
                if (j <= i) continue;

                int Zj = system.atomic_numbers()(j);
                if (!pair_matches(Zi, Zj, Z1, Z2)) continue;

                Vec3 dr = system.position(j) - system.position(i);
                dr += cell.col(0) * it->cell_shift[0];
                dr += cell.col(1) * it->cell_shift[1];
                dr += cell.col(2) * it->cell_shift[2];

                Scalar r = dr.norm();
                if (r >= cutoff_radius || r < 1e-10) continue;

                Scalar V, dV_dr;
                if (r < rm) {
                    Scalar d = r - r1;
                    V = 0.5 * k1 * d * d;
                    dV_dr = k1 * d;
                } else {
                    Scalar d = r - r2;
                    V = 0.5 * k2 * d * d;
                    dV_dr = k2 * d;
                }

                results.energy += V;

                if (compute_forces || compute_virial) {
                    Vec3 f = -(dV_dr / r) * dr;

                    if (compute_forces) {
                        system.forces().col(i) -= f.array();
                        system.forces().col(j) += f.array();
                    }

                    if (compute_virial) {
                        results.virial += dr * f.transpose();
                    }
                }
            }
        }

        return results;
    }
};

// ============================================================================
// R6 (r^-6) potential
// ============================================================================

/**
 * @brief r^-6 (London dispersion) potential
 *
 * V(r) = A / (r0 + r)^6
 *
 * Can be used as a dispersion correction. A < 0 gives attraction.
 */
class R6 : public PotentialBase<R6> {
public:
    Scalar A = 1.0;
    Scalar r0 = 0.0;
    Scalar cutoff_radius = 5.0;
    int Z1 = 0;
    int Z2 = 0;

    R6() = default;

    R6(Scalar A_, Scalar r0_, Scalar cutoff_, int Z1_ = 0, int Z2_ = 0)
        : A(A_), r0(r0_), cutoff_radius(cutoff_), Z1(Z1_), Z2(Z2_) {}

    Scalar cutoff() const { return cutoff_radius; }

    void bind_to(AtomicSystem&, NeighborList& nl) {
        nl.set_cutoff(std::max(nl.cutoff(), cutoff_radius));
    }

    PotentialResults compute(AtomicSystem& system, NeighborList& neighbors,
                             bool compute_forces = true,
                             bool compute_virial = true)
    {
        PotentialResults results;
        const Mat3& cell = system.cell();

        for (std::size_t i = 0; i < system.num_atoms(); ++i) {
            int Zi = system.atomic_numbers()(i);
            auto [begin, end] = neighbors.neighbors(i);

            for (auto it = begin; it != end; ++it) {
                std::size_t j = it->index;
                if (j <= i) continue;

                int Zj = system.atomic_numbers()(j);
                if (!pair_matches(Zi, Zj, Z1, Z2)) continue;

                Vec3 dr = system.position(j) - system.position(i);
                dr += cell.col(0) * it->cell_shift[0];
                dr += cell.col(1) * it->cell_shift[1];
                dr += cell.col(2) * it->cell_shift[2];

                Scalar r = dr.norm();
                if (r >= cutoff_radius || r < 1e-10) continue;

                Scalar rp = r0 + r;
                Scalar rp6 = rp * rp * rp * rp * rp * rp;
                Scalar V = A / rp6;
                results.energy += V;

                if (compute_forces || compute_virial) {
                    Scalar dV_dr = -6.0 * A / (rp6 * rp);
                    Vec3 f = -(dV_dr / r) * dr;

                    if (compute_forces) {
                        system.forces().col(i) -= f.array();
                        system.forces().col(j) += f.array();
                    }

                    if (compute_virial) {
                        results.virial += dr * f.transpose();
                    }
                }
            }
        }

        return results;
    }
};

} // namespace atomistica
