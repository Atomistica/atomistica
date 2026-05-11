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

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "../../include/atomistica/config.hpp"
#include "../../include/atomistica/core/atomic_system.hpp"
#include "../../include/atomistica/core/neighbor_list.hpp"
#include "../../include/atomistica/potentials/potential_base.hpp"
#include "../../include/atomistica/potentials/coulomb/coulomb.hpp"
#include "../../include/atomistica/potentials/coulomb/pme.hpp"

#include "../config.hpp"
#include "../registry.hpp"
#include "../coulomb_solver.hpp"

namespace atomistica {

// k_e = 1/(4πε₀) in eV·Å/e²
static constexpr double COULOMB_K = 14.3996447794;

// ---------------------------------------------------------------------------
// Helper: compute dr vector for a neighbor (handles PBC image shift)
// ---------------------------------------------------------------------------
inline Vec3 neighbor_dr(const AtomicSystem& sys, size_t i, size_t j,
                         const std::array<int,3>& cell_shift) {
    Vec3 shift(static_cast<Scalar>(cell_shift[0]),
               static_cast<Scalar>(cell_shift[1]),
               static_cast<Scalar>(cell_shift[2]));
    return sys.position(j) + sys.cell() * shift - sys.position(i);
}

// ---------------------------------------------------------------------------
// DirectCoulombWrapper  (O(N²), no NL needed)
// ---------------------------------------------------------------------------
class DirectCoulombWrapper : public CoulombSolver {
    DirectCoulomb solver_;
    double        k_eff_;

public:
    explicit DirectCoulombWrapper(const Config& cfg) {
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_.set_epsilon_r(eps_r);
        k_eff_ = COULOMB_K / eps_r;
    }

    Scalar cutoff() const override { return solver_.cutoff(); }
    void bind_to(AtomicSystem&, NeighborList&) override {}

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                              bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            solver_.set_charges(std::vector<Scalar>(q.data(), q.data() + q.size()));
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    // phi_i = k_eff * sum_{j≠i} q_j / r_ij  (minimum-image PBC)
    void compute_potential(AtomicSystem& sys, NeighborList&,
                           const ArrayX& q, ArrayX& phi) override {
        int n = static_cast<int>(sys.num_atoms());
        phi.resize(n);
        phi.setZero();
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (i == j) continue;
                Vec3 dr = sys.minimum_image(sys.position(j) - sys.position(i));
                double r = dr.norm();
                if (r < 1e-10) continue;
                phi[i] += q[j] / r;
            }
            phi[i] *= static_cast<Scalar>(k_eff_);
        }
    }
};
REGISTER_COULOMB("DirectCoulomb", DirectCoulombWrapper)

// ---------------------------------------------------------------------------
// CutoffCoulombWrapper
// ---------------------------------------------------------------------------
class CutoffCoulombWrapper : public CoulombSolver {
    CutoffCoulomb solver_;
    double        k_eff_;

public:
    explicit CutoffCoulombWrapper(const Config& cfg) {
        double rc    = cfg.get_or<double>("cutoff", 10.0);
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_.set_cutoff(rc);
        solver_.set_epsilon_r(eps_r);
        k_eff_ = COULOMB_K / eps_r;
    }

    Scalar cutoff() const override { return solver_.cutoff(); }
    void bind_to(AtomicSystem&, NeighborList&) override {}

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                              bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            solver_.set_charges(std::vector<Scalar>(q.data(), q.data() + q.size()));
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    // phi_i = k_eff * sum_{j, r_ij<rc} q_j / r_ij  (using NL)
    void compute_potential(AtomicSystem& sys, NeighborList& nl,
                           const ArrayX& q, ArrayX& phi) override {
        int n = static_cast<int>(sys.num_atoms());
        phi.resize(n);
        phi.setZero();
        double rc2 = solver_.cutoff() * solver_.cutoff();

        for (size_t i = 0; i < sys.num_atoms(); ++i) {
            auto [beg, end] = nl.neighbors(i);
            for (auto it = beg; it != end; ++it) {
                size_t j = it->index;
                Vec3 dr = neighbor_dr(sys, i, j, it->cell_shift);
                double r2 = dr.squaredNorm();
                if (r2 >= rc2 || r2 < 1e-20) continue;
                phi[i] += static_cast<Scalar>(q[j] / std::sqrt(r2));
            }
            phi[i] *= static_cast<Scalar>(k_eff_);
        }
    }
};
REGISTER_COULOMB("CutoffCoulomb", CutoffCoulombWrapper)

// ---------------------------------------------------------------------------
// WolfCoulombWrapper  (Damped Shifted Force method)
// ---------------------------------------------------------------------------
class WolfCoulombWrapper : public CoulombSolver {
    WolfCoulomb solver_;
    double k_eff_;
    double alpha_;
    double rc_;
    double shift_pot_;  // erfc(alpha*rc)/rc — precomputed in bind_to

public:
    explicit WolfCoulombWrapper(const Config& cfg) {
        rc_    = cfg.get_or<double>("cutoff", 10.0);
        alpha_ = cfg.get_or<double>("alpha", 0.2);
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_.set_cutoff(rc_);
        solver_.set_alpha(alpha_);
        solver_.set_epsilon_r(eps_r);
        k_eff_ = COULOMB_K / eps_r;
        // Actual alpha may be auto-selected by library; retrieve it.
        alpha_    = solver_.alpha();
        shift_pot_ = std::erfc(alpha_ * rc_) / rc_;
    }

    void bind_to(AtomicSystem&, NeighborList&) override {
        // Refresh in case alpha was auto-selected at construction.
        alpha_    = solver_.alpha();
        shift_pot_ = std::erfc(alpha_ * rc_) / rc_;
    }

    Scalar cutoff() const override { return solver_.cutoff(); }

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                              bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            solver_.set_charges(std::vector<Scalar>(q.data(), q.data() + q.size()));
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    // DSF electrostatic potential (no self-term; suitable for QEq SCF).
    // phi_i = k_eff * sum_{j≠i, r<rc} q_j * [erfc(alpha*r)/r - erfc(alpha*rc)/rc]
    void compute_potential(AtomicSystem& sys, NeighborList& nl,
                           const ArrayX& q, ArrayX& phi) override {
        int n = static_cast<int>(sys.num_atoms());
        phi.resize(n);
        phi.setZero();
        double rc2 = rc_ * rc_;

        for (size_t i = 0; i < sys.num_atoms(); ++i) {
            auto [beg, end] = nl.neighbors(i);
            for (auto it = beg; it != end; ++it) {
                size_t j = it->index;
                Vec3 dr = neighbor_dr(sys, i, j, it->cell_shift);
                double r2 = dr.squaredNorm();
                if (r2 >= rc2 || r2 < 1e-20) continue;
                double r = std::sqrt(r2);
                double contrib = std::erfc(alpha_ * r) / r - shift_pot_;
                phi[i] += static_cast<Scalar>(q[j] * contrib);
            }
            phi[i] *= static_cast<Scalar>(k_eff_);
        }
    }
};
REGISTER_COULOMB("WolfCoulomb", WolfCoulombWrapper)

// ---------------------------------------------------------------------------
// PMECoulombWrapper  (Particle Mesh Ewald)
// compute_potential() uses a CutoffCoulomb approximation for the SCF loop
// since the library does not expose a separate potential calculation for PME.
// ---------------------------------------------------------------------------
class PMECoulombWrapper : public CoulombSolver {
    PMECoulomb solver_;
    double     rc_;
    double     k_eff_;
    double     alpha_;
    bool       potential_warning_printed_ = false;

public:
    explicit PMECoulombWrapper(const Config& cfg) {
        rc_    = cfg.get_or<double>("cutoff", 10.0);
        int gx = cfg.get_or<int>("grid_x", 32);
        int gy = cfg.get_or<int>("grid_y", 32);
        int gz = cfg.get_or<int>("grid_z", 32);
        int ord = cfg.get_or<int>("order", 4);
        alpha_ = cfg.get_or<double>("alpha", 0.0);
        double eps_r = cfg.get_or<double>("epsilon_r", 1.0);
        solver_ = PMECoulomb(rc_, gx, gy, gz, ord, alpha_);
        k_eff_  = COULOMB_K / eps_r;
    }

    Scalar cutoff() const override { return solver_.cutoff(); }
    void bind_to(AtomicSystem&, NeighborList&) override {}

    PotentialResults compute(AtomicSystem& sys, NeighborList& nl,
                              bool cf = true, bool cv = true) override {
        if (sys.properties().has("charges")) {
            const ArrayX& q = sys.properties().get<ArrayX>("charges");
            solver_.set_charges(std::vector<Scalar>(q.data(), q.data() + q.size()));
        }
        return solver_.compute(sys, nl, cf, cv);
    }

    // Approximation: use real-space erfc sum only (omits k-space correction).
    // This is exact for ChargeEquilibration when rc is large enough relative
    // to the system; the k-space error decays exponentially with alpha*rc.
    void compute_potential(AtomicSystem& sys, NeighborList& nl,
                           const ArrayX& q, ArrayX& phi) override {
        if (!potential_warning_printed_) {
            std::fprintf(stderr,
                "PMECoulomb::compute_potential: using real-space-only "
                "approximation for SCF (k-space correction omitted).\n");
            potential_warning_printed_ = true;
        }
        int n = static_cast<int>(sys.num_atoms());
        phi.resize(n);
        phi.setZero();
        double alpha = (alpha_ > 0.0) ? alpha_
                                       : std::sqrt(12.0 * std::log(10.0)) / rc_;
        double rc2 = rc_ * rc_;
        double shift = std::erfc(alpha * rc_) / rc_;

        for (size_t i = 0; i < sys.num_atoms(); ++i) {
            auto [beg, end] = nl.neighbors(i);
            for (auto it = beg; it != end; ++it) {
                size_t j = it->index;
                Vec3 dr = neighbor_dr(sys, i, j, it->cell_shift);
                double r2 = dr.squaredNorm();
                if (r2 >= rc2 || r2 < 1e-20) continue;
                double r = std::sqrt(r2);
                phi[i] += static_cast<Scalar>(q[j] * (std::erfc(alpha*r)/r - shift));
            }
            phi[i] *= static_cast<Scalar>(k_eff_);
        }
    }
};
REGISTER_COULOMB("PMECoulomb", PMECoulombWrapper)

} // namespace atomistica
