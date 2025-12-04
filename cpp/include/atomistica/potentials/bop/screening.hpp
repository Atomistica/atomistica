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
#include <vector>

#include "../../config.hpp"

namespace atomistica {

/**
 * @brief Screening parameters for a pair type
 *
 * Screening function based on:
 * - Baskes et al., Modelling Simul. Mater. Sci. Eng. 2, 505 (1994)
 * - Pastewka et al., Phys. Rev. B 78, 161402(R) (2008)
 */
struct ScreeningParams {
    // Screening cutoff parameters
    Scalar Cmin = 1.0;          // Inner screening parameter (fully unscreened if C > Cmax)
    Scalar Cmax = 3.0;          // Outer screening parameter (fully screened if C < Cmin)
    Scalar dC = 2.0;            // Cmax - Cmin (precomputed)

    // Cutoff regions for screened potential
    Scalar cut_in_l = 0.0;      // Inner cutoff start (below: unscreened, fc=1)
    Scalar cut_in_h = 0.0;      // Inner cutoff end
    Scalar cut_out_l = 0.0;     // Outer cutoff start (screening region)
    Scalar cut_out_h = 0.0;     // Outer cutoff end (above: fc=0)
    Scalar cut_bo_l = 0.0;      // Bond-order cutoff start
    Scalar cut_bo_h = 0.0;      // Bond-order cutoff end

    // Precomputed
    Scalar max_cut_sq = 0.0;    // Maximum cutoff squared
    Scalar C_dr_cut = 0.0;      // Screening neighbor cutoff factor

    // Thresholds
    Scalar screening_threshold = std::log(1e-6);
    Scalar dot_threshold = 1e-10;

    void precompute() {
        dC = Cmax - Cmin;
        max_cut_sq = cut_out_h * cut_out_h;
        // C_dr_cut determines how far to look for screening atoms
        // Based on geometry, atoms within C_dr_cut * r_ij^2 can screen
        C_dr_cut = 1.0 + 2.0 / (Cmin + 1e-10);
    }
};

/**
 * @brief Result of screening function evaluation
 */
struct ScreeningResult {
    Scalar S = 1.0;             // Screening function value [0, 1]
    Scalar dS_drij = 0.0;       // Derivative w.r.t. r_ij magnitude
    bool fully_screened = false; // Bond is completely screened
};

/**
 * @brief Data for a screening neighbor (atom k that screens bond i-j)
 */
struct ScreeningNeighbor {
    std::size_t k;              // Atom index
    Scalar dS_drik = 0.0;       // dS/dr_ik (scaled by 1/r_ik)
    Scalar dS_drjk = 0.0;       // dS/dr_jk (scaled by 1/r_jk)
};

/**
 * @brief Compute the screening function for bond i-j
 *
 * The screening function S reduces the interaction when intermediate atoms
 * lie between atoms i and j. This implements the algorithm from:
 * Pastewka et al., Phys. Rev. B 78, 161402(R) (2008)
 *
 * The screening parameter C for atom k is:
 *   C = (2*(x_ik + x_jk) - (x_ik - x_jk)^2 - 1) / (1 - (x_ik - x_jk)^2)
 * where x_ik = r_ik^2 / r_ij^2, x_jk = r_jk^2 / r_ij^2
 *
 * - C > Cmax: atom k does not screen (no contribution)
 * - C < Cmin: bond is fully screened (S = 0)
 * - Cmin < C < Cmax: partial screening
 *
 * @param params Screening parameters
 * @param rij_vec Distance vector r_j - r_i
 * @param rij_sq Squared distance |r_ij|^2
 * @param neighbors List of potential screening atoms with their r_ik vectors
 * @param screening_neighbors Output: list of atoms that contribute to screening
 * @return Screening result with S value and derivatives
 */
template<typename NeighborIterator>
ScreeningResult compute_screening(
    const ScreeningParams& params,
    const Vec3& rij_vec,
    Scalar rij_sq,
    NeighborIterator nb_begin,
    NeighborIterator nb_end,
    const std::function<Vec3(const typename std::iterator_traits<NeighborIterator>::value_type&)>& get_rik,
    std::vector<ScreeningNeighbor>& screening_neighbors)
{
    ScreeningResult result;
    result.S = 0.0;  // Using polynomial screening: S = sum of (Cmax-C)^2/(C-Cmin)^2
    result.dS_drij = 0.0;
    result.fully_screened = false;

    screening_neighbors.clear();

    const Scalar C_dr_cut_rij_sq = params.C_dr_cut * rij_sq;

    for (auto it = nb_begin; it != nb_end; ++it) {
        Vec3 rik_vec = get_rik(*it);
        Scalar rik_sq = rik_vec.squaredNorm();

        // Skip atoms too far away to screen
        if (rik_sq >= C_dr_cut_rij_sq) continue;

        // Compute r_jk = r_ik - r_ij
        Vec3 rjk_vec = rik_vec - rij_vec;
        Scalar rjk_sq = rjk_vec.squaredNorm();

        // Compute dot products for geometry check
        Scalar dot_ij_ik = rij_vec.dot(rik_vec);
        Scalar dot_ij_jk = rij_vec.dot(rjk_vec);

        // Atom k must be "between" i and j (geometrically)
        // dot_ij_ik > 0: k is on the j-side of i
        // dot_ij_jk < 0: k is on the i-side of j
        if (dot_ij_ik <= params.dot_threshold || dot_ij_jk >= -params.dot_threshold) {
            continue;
        }

        // Compute screening parameter C
        Scalar xik = rik_sq / rij_sq;
        Scalar xjk = rjk_sq / rij_sq;

        Scalar xik_m_xjk = xik - xjk;
        Scalar xik_p_xjk = xik + xjk;

        Scalar denom = 1.0 - xik_m_xjk * xik_m_xjk;
        if (std::abs(denom) < 1e-15) continue;  // Degenerate case

        Scalar fac = 1.0 / denom;
        Scalar C = (2.0 * xik_p_xjk - xik_m_xjk * xik_m_xjk - 1.0) * fac;

        if (C <= params.Cmin) {
            // Fully screened by this atom
            result.fully_screened = true;
            result.S = params.screening_threshold;  // Very negative = fully screened
            return result;
        }

        if (C < params.Cmax) {
            // Partial screening - accumulate contribution
            Scalar Cmax_C = params.Cmax - C;
            Scalar C_Cmin = C - params.Cmin;

            // Polynomial screening: S -= (Cmax-C)^2 / (C-Cmin)^2
            Scalar ratio = Cmax_C / C_Cmin;
            result.S -= ratio * ratio;

            // Derivatives of C with respect to distances
            // dC/d(xik) and dC/d(xjk) then chain rule to distances
            Scalar dCdxik = 4.0 * xik * fac * (1.0 + (C - 1.0) * xik_m_xjk);
            Scalar dCdxjk = 4.0 * xjk * fac * (1.0 - (C - 1.0) * xik_m_xjk);
            Scalar dCdrij_sq = -(dCdxik + dCdxjk);  // d/d(rij^2)

            // dS/dC = 2 * (Cmax-C) * dC / (C-Cmin)^3
            Scalar dSdC = 2.0 * Cmax_C * params.dC / (C_Cmin * C_Cmin * C_Cmin);

            // dS/dr_ij (via chain rule through rij^2)
            result.dS_drij += dSdC * dCdrij_sq;

            // Store screening neighbor with derivatives
            ScreeningNeighbor sn;
            sn.k = it->index;
            // dS/dr_ik = dS/dC * dC/d(xik) * d(xik)/d(rik^2) * d(rik^2)/d(rik)
            //          = dSdC * dCdxik * (1/rij^2) * 2*rik
            // We store dS/d(rik) / rik = dSdC * dCdxik * 2 / rij^2
            sn.dS_drik = dSdC * dCdxik * 2.0 / rij_sq;
            sn.dS_drjk = dSdC * dCdxjk * 2.0 / rij_sq;

            screening_neighbors.push_back(sn);
        }
        // If C >= Cmax, atom k doesn't contribute to screening
    }

    // Check if accumulated screening exceeds threshold
    if (result.S < params.screening_threshold) {
        result.fully_screened = true;
    }

    // Convert dS_drij from d/d(rij^2) to d/d(rij)
    // d(rij^2)/d(rij) = 2*rij, so dS/d(rij) = dS/d(rij^2) * 2 * rij
    // We need to divide by rij to get the per-unit-distance derivative
    // that will be multiplied by unit vector later
    Scalar rij = std::sqrt(rij_sq);
    if (rij > 1e-10) {
        result.dS_drij *= 2.0 * rij;
    }

    return result;
}

/**
 * @brief Simpler screening computation without storing neighbors (for cutoff only)
 */
inline ScreeningResult compute_screening_simple(
    const ScreeningParams& params,
    const Vec3& rij_vec,
    Scalar rij_sq,
    const std::vector<Vec3>& rik_vectors)
{
    ScreeningResult result;
    result.S = 0.0;
    result.dS_drij = 0.0;
    result.fully_screened = false;

    const Scalar C_dr_cut_rij_sq = params.C_dr_cut * rij_sq;

    for (const auto& rik_vec : rik_vectors) {
        Scalar rik_sq = rik_vec.squaredNorm();
        if (rik_sq >= C_dr_cut_rij_sq) continue;

        Vec3 rjk_vec = rik_vec - rij_vec;
        Scalar rjk_sq = rjk_vec.squaredNorm();

        Scalar dot_ij_ik = rij_vec.dot(rik_vec);
        Scalar dot_ij_jk = rij_vec.dot(rjk_vec);

        if (dot_ij_ik <= params.dot_threshold || dot_ij_jk >= -params.dot_threshold) {
            continue;
        }

        Scalar xik = rik_sq / rij_sq;
        Scalar xjk = rjk_sq / rij_sq;
        Scalar xik_m_xjk = xik - xjk;
        Scalar xik_p_xjk = xik + xjk;

        Scalar denom = 1.0 - xik_m_xjk * xik_m_xjk;
        if (std::abs(denom) < 1e-15) continue;

        Scalar fac = 1.0 / denom;
        Scalar C = (2.0 * xik_p_xjk - xik_m_xjk * xik_m_xjk - 1.0) * fac;

        if (C <= params.Cmin) {
            result.fully_screened = true;
            result.S = params.screening_threshold;
            return result;
        }

        if (C < params.Cmax) {
            Scalar Cmax_C = params.Cmax - C;
            Scalar C_Cmin = C - params.Cmin;
            Scalar ratio = Cmax_C / C_Cmin;
            result.S -= ratio * ratio;
        }
    }

    if (result.S < params.screening_threshold) {
        result.fully_screened = true;
    }

    return result;
}

} // namespace atomistica
