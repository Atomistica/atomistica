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
#include <vector>

#include "../../config.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../potential_base.hpp"

namespace atomistica {

// Physical constants in eV and Angstrom units
// Coulomb constant k_e = 1/(4*pi*eps_0) in eV*Angstrom/e^2
// = 14.3996 eV*Angstrom
constexpr Scalar COULOMB_CONST = 14.3996447794;  // eV*Angstrom/e^2

/**
 * @brief Direct Coulomb summation for non-periodic systems
 *
 * Computes electrostatic energy and forces using direct pairwise summation:
 *   E = k_e * sum_{i<j} q_i * q_j / r_ij
 *
 * This is O(N^2) and should only be used for small, non-periodic systems.
 *
 * Reference: Standard electrostatics
 */
class DirectCoulomb : public PotentialBase<DirectCoulomb> {
public:
    /**
     * @brief Construct with optional dielectric constant
     * @param epsilon_r Relative dielectric constant (default: 1.0 = vacuum)
     */
    explicit DirectCoulomb(Scalar epsilon_r = 1.0)
        : epsilon_r_(epsilon_r)
        , k_eff_(COULOMB_CONST / epsilon_r) {}

    /**
     * @brief Set relative dielectric constant
     */
    void set_epsilon_r(Scalar epsilon_r) {
        epsilon_r_ = epsilon_r;
        k_eff_ = COULOMB_CONST / epsilon_r;
    }

    Scalar epsilon_r() const { return epsilon_r_; }

    /**
     * @brief Set charges for all atoms
     * @param charges Vector of charges (in units of e)
     */
    void set_charges(const std::vector<Scalar>& charges) {
        charges_ = charges;
    }

    /**
     * @brief Get charges
     */
    const std::vector<Scalar>& charges() const { return charges_; }

    // CRTP implementation
    Scalar cutoff_impl() const {
        // Direct Coulomb has infinite cutoff (no truncation)
        return std::numeric_limits<Scalar>::infinity();
    }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& /*neighbors*/,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    Scalar epsilon_r_ = 1.0;
    Scalar k_eff_;
    std::vector<Scalar> charges_;
};

/**
 * @brief Cutoff Coulomb with hard truncation
 *
 * WARNING: This method introduces discontinuities at the cutoff!
 * Use WolfCoulomb for better behavior.
 *
 * E = k_e * sum_{i<j, r_ij < r_c} q_i * q_j / r_ij
 */
class CutoffCoulomb : public PotentialBase<CutoffCoulomb> {
public:
    /**
     * @brief Construct with cutoff and optional dielectric constant
     * @param cutoff Cutoff radius in Angstrom
     * @param epsilon_r Relative dielectric constant
     */
    explicit CutoffCoulomb(Scalar cutoff = 10.0, Scalar epsilon_r = 1.0)
        : cutoff_(cutoff)
        , cutoff_sq_(cutoff * cutoff)
        , epsilon_r_(epsilon_r)
        , k_eff_(COULOMB_CONST / epsilon_r) {}

    void set_cutoff(Scalar cutoff) {
        cutoff_ = cutoff;
        cutoff_sq_ = cutoff * cutoff;
    }

    void set_epsilon_r(Scalar epsilon_r) {
        epsilon_r_ = epsilon_r;
        k_eff_ = COULOMB_CONST / epsilon_r;
    }

    Scalar epsilon_r() const { return epsilon_r_; }

    void set_charges(const std::vector<Scalar>& charges) {
        charges_ = charges;
    }

    const std::vector<Scalar>& charges() const { return charges_; }

    // CRTP implementation
    Scalar cutoff_impl() const { return cutoff_; }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    Scalar cutoff_ = 10.0;
    Scalar cutoff_sq_;
    Scalar epsilon_r_ = 1.0;
    Scalar k_eff_;
    std::vector<Scalar> charges_;
};

/**
 * @brief Wolf summation for Coulomb interactions
 *
 * Implements the Wolf summation method which provides a smooth, charge-neutral
 * truncation of the Coulomb interaction. Much better than hard cutoff.
 *
 * E = k_e * sum_{i<j, r_ij < r_c} q_i * q_j * [erfc(alpha*r_ij)/r_ij
 *                                             - erfc(alpha*r_c)/r_c
 *                                             + (erfc(alpha*r_c)/r_c^2
 *                                               + 2*alpha/sqrt(pi)*exp(-alpha^2*r_c^2)/r_c)
 *                                               * (r_ij - r_c)]
 *     - k_e * sum_i q_i^2 * [erfc(alpha*r_c)/(2*r_c) + alpha/sqrt(pi)]
 *
 * Reference: Wolf, Keblinski, Phillpot, Eggebrecht, J. Chem. Phys. 110, 8254 (1999)
 *
 * The damped shifted force (DSF) variant is used, which ensures both energy
 * and force are continuous at the cutoff.
 *
 * Reference: Fennell, Gezelter, J. Chem. Phys. 124, 234104 (2006)
 */
class WolfCoulomb : public PotentialBase<WolfCoulomb> {
public:
    /**
     * @brief Construct Wolf Coulomb solver
     * @param cutoff Cutoff radius in Angstrom
     * @param alpha Damping parameter (1/Angstrom). If <= 0, auto-computed.
     * @param epsilon_r Relative dielectric constant
     */
    explicit WolfCoulomb(Scalar cutoff = 10.0, Scalar alpha = 0.0, Scalar epsilon_r = 1.0)
        : cutoff_(cutoff)
        , cutoff_sq_(cutoff * cutoff)
        , epsilon_r_(epsilon_r)
        , k_eff_(COULOMB_CONST / epsilon_r) {
        set_alpha(alpha);
    }

    void set_cutoff(Scalar cutoff) {
        cutoff_ = cutoff;
        cutoff_sq_ = cutoff * cutoff;
        precompute();
    }

    /**
     * @brief Set damping parameter
     * @param alpha Damping parameter. If <= 0, auto-computed as sqrt(ln(10)*12)/cutoff
     */
    void set_alpha(Scalar alpha) {
        if (alpha <= 0) {
            // Auto-compute alpha for good convergence
            // This gives erfc(alpha*cutoff) ~ 1e-6
            alpha_ = std::sqrt(12.0 * std::log(10.0)) / cutoff_;
        } else {
            alpha_ = alpha;
        }
        precompute();
    }

    void set_epsilon_r(Scalar epsilon_r) {
        epsilon_r_ = epsilon_r;
        k_eff_ = COULOMB_CONST / epsilon_r;
    }

    Scalar cutoff() const { return cutoff_; }
    Scalar alpha() const { return alpha_; }
    Scalar epsilon_r() const { return epsilon_r_; }

    void set_charges(const std::vector<Scalar>& charges) {
        charges_ = charges;
    }

    const std::vector<Scalar>& charges() const { return charges_; }

    // CRTP implementation
    Scalar cutoff_impl() const { return cutoff_; }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    void precompute();

    Scalar cutoff_ = 10.0;
    Scalar cutoff_sq_;
    Scalar alpha_ = 0.0;
    Scalar epsilon_r_ = 1.0;
    Scalar k_eff_;
    std::vector<Scalar> charges_;

    // Precomputed quantities at cutoff
    Scalar erfc_alpha_rc_;       // erfc(alpha * r_c)
    Scalar erfc_alpha_rc_over_rc_;  // erfc(alpha * r_c) / r_c
    Scalar shift_potential_;     // Potential shift for DSF
    Scalar shift_force_;         // Force shift for DSF
    Scalar self_energy_factor_;  // Self-energy per q^2
};

// ============================================================================
// Implementation
// ============================================================================

inline PotentialResults DirectCoulomb::compute_impl(
    AtomicSystem& system,
    NeighborList& /*neighbors*/,
    bool compute_forces,
    bool compute_virial)
{
    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();

    if (charges_.size() != num_atoms) {
        throw std::runtime_error("DirectCoulomb: charges array size mismatch");
    }

    // Direct O(N^2) summation - all pairs
    for (std::size_t i = 0; i < num_atoms; ++i) {
        Scalar qi = charges_[i];
        if (std::abs(qi) < 1e-15) continue;

        Vec3 ri = system.position(i).matrix();

        for (std::size_t j = i + 1; j < num_atoms; ++j) {
            Scalar qj = charges_[j];
            if (std::abs(qj) < 1e-15) continue;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;

            // For non-periodic systems, apply minimum image if PBC is set
            if (system.pbc()[0] || system.pbc()[1] || system.pbc()[2]) {
                dr = system.minimum_image(dr);
            }

            Scalar r = dr.norm();
            if (r < 1e-10) continue;  // Skip overlapping atoms

            Scalar inv_r = 1.0 / r;

            // Energy: k_eff * qi * qj / r
            Scalar pair_energy = k_eff_ * qi * qj * inv_r;
            results.energy += pair_energy;

            if (compute_forces || compute_virial) {
                // Force: -dE/dr * r_hat = k_eff * qi * qj / r^2 * r_hat
                // Force on i is -k_eff * qi * qj / r^3 * dr (dr points from i to j)
                Scalar force_mag = -k_eff_ * qi * qj * inv_r * inv_r * inv_r;
                Vec3 force = force_mag * dr;

                if (compute_forces) {
                    // Force on i from j
                    system.forces().col(i) += force.array();
                    // Newton's third law
                    system.forces().col(j) -= force.array();
                }

                if (compute_virial) {
                    // Virial contribution: r_ij * F_i (consistent with BOPKernel sign)
                    results.virial += dr * force.transpose();
                }
            }
        }
    }

    return results;
}

inline PotentialResults CutoffCoulomb::compute_impl(
    AtomicSystem& system,
    NeighborList& neighbors,
    bool compute_forces,
    bool compute_virial)
{
    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();
    const Mat3& cell = system.cell();

    if (charges_.size() != num_atoms) {
        throw std::runtime_error("CutoffCoulomb: charges array size mismatch");
    }

    // Use neighbor list
    for (std::size_t i = 0; i < num_atoms; ++i) {
        Scalar qi = charges_[i];
        if (std::abs(qi) < 1e-15) continue;

        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;

            Scalar qj = charges_[j];
            if (std::abs(qj) < 1e-15) continue;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq >= cutoff_sq_ || r_sq < 1e-20) continue;

            Scalar r = std::sqrt(r_sq);
            Scalar inv_r = 1.0 / r;

            // Energy: k_eff * qi * qj / r (halved for full neighbor list)
            Scalar pair_energy = 0.5 * k_eff_ * qi * qj * inv_r;
            results.energy += pair_energy;

            if (compute_forces || compute_virial) {
                // Force on i: -dV/dr_i = k_eff * qi * qj / r^3 * (rj - ri)
                // force_over_r = k*qi*qj/r^3 (note: positive sign here,
                // consistent with WolfCoulomb gradient convention so that
                // forces(i) -= force gives F_i = -force = k*qi*qj/r^3 * dr)
                Scalar force_over_r = k_eff_ * qi * qj * inv_r * inv_r * inv_r;
                Vec3 force = force_over_r * dr;

                if (compute_forces) {
                    system.forces().col(i) -= force.array();
                }

                if (compute_virial) {
                    results.virial -= 0.5 * dr * force.transpose();
                }
            }
        }
    }

    return results;
}

inline void WolfCoulomb::precompute() {
    // Constants
    const Scalar sqrt_pi = std::sqrt(M_PI);
    const Scalar two_alpha_sqrt_pi = 2.0 * alpha_ / sqrt_pi;

    // erfc at cutoff
    erfc_alpha_rc_ = std::erfc(alpha_ * cutoff_);
    erfc_alpha_rc_over_rc_ = erfc_alpha_rc_ / cutoff_;

    // Exponential at cutoff
    Scalar exp_alpha_rc_sq = std::exp(-alpha_ * alpha_ * cutoff_ * cutoff_);

    // DSF shift terms (ensures energy and force are zero at cutoff)
    // Potential shift: erfc(alpha*rc)/rc
    shift_potential_ = erfc_alpha_rc_over_rc_;

    // Force shift: d/dr[erfc(alpha*r)/r] at r = rc
    // = -erfc(alpha*rc)/rc^2 - 2*alpha/sqrt(pi)*exp(-alpha^2*rc^2)/rc
    shift_force_ = erfc_alpha_rc_over_rc_ / cutoff_
                 + two_alpha_sqrt_pi * exp_alpha_rc_sq / cutoff_;

    // Self-energy factor: erfc(alpha*rc)/(2*rc) + alpha/sqrt(pi)
    self_energy_factor_ = erfc_alpha_rc_over_rc_ / 2.0 + alpha_ / sqrt_pi;
}

inline PotentialResults WolfCoulomb::compute_impl(
    AtomicSystem& system,
    NeighborList& neighbors,
    bool compute_forces,
    bool compute_virial)
{
    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();
    const Mat3& cell = system.cell();

    if (charges_.size() != num_atoms) {
        throw std::runtime_error("WolfCoulomb: charges array size mismatch");
    }

    const Scalar sqrt_pi = std::sqrt(M_PI);
    const Scalar two_alpha_sqrt_pi = 2.0 * alpha_ / sqrt_pi;

    // Self-energy contribution
    Scalar q_sum_sq = 0.0;
    for (std::size_t i = 0; i < num_atoms; ++i) {
        q_sum_sq += charges_[i] * charges_[i];
    }
    results.energy -= k_eff_ * self_energy_factor_ * q_sum_sq;

    // Pair interactions
    for (std::size_t i = 0; i < num_atoms; ++i) {
        Scalar qi = charges_[i];
        if (std::abs(qi) < 1e-15) continue;

        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;

            Scalar qj = charges_[j];
            if (std::abs(qj) < 1e-15) continue;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq >= cutoff_sq_ || r_sq < 1e-20) continue;

            Scalar r = std::sqrt(r_sq);
            Scalar inv_r = 1.0 / r;

            // erfc and exponential at this distance
            Scalar erfc_alpha_r = std::erfc(alpha_ * r);
            Scalar exp_alpha_r_sq = std::exp(-alpha_ * alpha_ * r_sq);

            // DSF energy: erfc(alpha*r)/r - shift_potential + shift_force * (r - rc)
            Scalar phi_dsf = erfc_alpha_r * inv_r
                           - shift_potential_
                           + shift_force_ * (r - cutoff_);

            // Energy (halved for full neighbor list)
            results.energy += 0.5 * k_eff_ * qi * qj * phi_dsf;

            if (compute_forces || compute_virial) {
                // Force derivative: d(phi_dsf)/dr
                // d/dr[erfc(alpha*r)/r] = -erfc(alpha*r)/r^2 - 2*alpha/sqrt(pi)*exp(-alpha^2*r^2)/r
                Scalar dphi_dr = -erfc_alpha_r * inv_r * inv_r
                               - two_alpha_sqrt_pi * exp_alpha_r_sq * inv_r
                               + shift_force_;

                // Force magnitude: -k_eff * qi * qj * dphi_dr / r
                // Force direction: -r_hat = -dr/r
                Scalar force_over_r = -k_eff_ * qi * qj * dphi_dr * inv_r;
                Vec3 force = force_over_r * dr;

                if (compute_forces) {
                    // Full neighbor list: only add to atom i
                    system.forces().col(i) -= force.array();
                }

                if (compute_virial) {
                    results.virial -= 0.5 * dr * force.transpose();
                }
            }
        }
    }

    return results;
}

} // namespace atomistica
