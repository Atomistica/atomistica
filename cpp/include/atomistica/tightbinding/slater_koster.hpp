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

#include "../config.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Slater-Koster transformation functions
 *
 * These functions transform tabulated SK integrals (sss, sps, pps, etc.)
 * to Cartesian orbital basis using direction cosines.
 *
 * Orbital indices:
 *   1: s
 *   2: px, 3: py, 4: pz
 *   5: dxy, 6: dyz, 7: dzx, 8: dx2-y2, 9: d3z2-r2
 *
 * SK integral indices in array:
 *   0: dds (d-d sigma)
 *   1: ddp (d-d pi)
 *   2: ddd (d-d delta)
 *   3: pds (p-d sigma)
 *   4: pdp (p-d pi)
 *   5: pps (p-p sigma)
 *   6: ppp (p-p pi)
 *   7: sds (s-d sigma)
 *   8: sps (s-p sigma)
 *   9: sss (s-s sigma)
 */

constexpr Scalar SQRT3 = 1.7320508075688772935;

/**
 * @brief Transform SK integrals to Cartesian matrix element
 *
 * @param a First orbital index (1-9)
 * @param b Second orbital index (1-9)
 * @param c Direction cosines [l, m, n] from atom i to j
 * @param sk SK integrals array (dds, ddp, ddd, pds, pdp, pps, ppp, sds, sps, sss)
 * @return Transformed matrix element H_ab or S_ab
 */
inline Scalar transform_orb(int a, int b, const Vec3& c,
                            const std::array<Scalar, NUM_SK_INTEGRALS>& sk) {
    // Extract direction cosines
    const Scalar l = c[0];
    const Scalar m = c[1];
    const Scalar n = c[2];

    // Precompute powers
    const Scalar ll = l * l;
    const Scalar mm = m * m;
    const Scalar nn = n * n;

    // Extract SK integrals
    const Scalar dds = sk[0];
    const Scalar ddp = sk[1];
    const Scalar ddd = sk[2];
    const Scalar pds = sk[3];
    const Scalar pdp = sk[4];
    const Scalar pps = sk[5];
    const Scalar ppp = sk[6];
    const Scalar sds = sk[7];
    const Scalar sps = sk[8];
    const Scalar sss = sk[9];

    // Handle parity correction for swapped orbitals
    // When a > b, apply parity factor (-1)^(l_a + l_b)
    bool swapped = false;
    if (a > b) {
        std::swap(a, b);
        swapped = true;
    }

    Scalar result = 0.0;

    // s-s interaction (a=1, b=1)
    if (a == 1 && b == 1) {
        result = sss;
    }
    // s-p interactions (a=1, b=2,3,4)
    else if (a == 1 && b == 2) {  // s-px
        result = l * sps;
    }
    else if (a == 1 && b == 3) {  // s-py
        result = m * sps;
    }
    else if (a == 1 && b == 4) {  // s-pz
        result = n * sps;
    }
    // s-d interactions (a=1, b=5,6,7,8,9)
    else if (a == 1 && b == 5) {  // s-dxy
        result = SQRT3 * l * m * sds;
    }
    else if (a == 1 && b == 6) {  // s-dyz
        result = SQRT3 * m * n * sds;
    }
    else if (a == 1 && b == 7) {  // s-dzx
        result = SQRT3 * n * l * sds;
    }
    else if (a == 1 && b == 8) {  // s-dx2-y2
        result = 0.5 * SQRT3 * (ll - mm) * sds;
    }
    else if (a == 1 && b == 9) {  // s-d3z2-r2
        result = (nn - 0.5 * (ll + mm)) * sds;
    }
    // p-p interactions (a=2,3,4, b=2,3,4)
    else if (a == 2 && b == 2) {  // px-px
        result = ll * pps + (1.0 - ll) * ppp;
    }
    else if (a == 2 && b == 3) {  // px-py
        result = l * m * (pps - ppp);
    }
    else if (a == 2 && b == 4) {  // px-pz
        result = l * n * (pps - ppp);
    }
    else if (a == 3 && b == 3) {  // py-py
        result = mm * pps + (1.0 - mm) * ppp;
    }
    else if (a == 3 && b == 4) {  // py-pz
        result = m * n * (pps - ppp);
    }
    else if (a == 4 && b == 4) {  // pz-pz
        result = nn * pps + (1.0 - nn) * ppp;
    }
    // p-d interactions (a=2,3,4, b=5,6,7,8,9)
    else if (a == 2 && b == 5) {  // px-dxy
        result = SQRT3 * ll * m * pds + m * (1.0 - 2.0 * ll) * pdp;
    }
    else if (a == 2 && b == 6) {  // px-dyz
        result = SQRT3 * l * m * n * pds - 2.0 * l * m * n * pdp;
    }
    else if (a == 2 && b == 7) {  // px-dzx
        result = SQRT3 * ll * n * pds + n * (1.0 - 2.0 * ll) * pdp;
    }
    else if (a == 2 && b == 8) {  // px-dx2-y2
        result = 0.5 * SQRT3 * l * (ll - mm) * pds + l * (1.0 - ll + mm) * pdp;
    }
    else if (a == 2 && b == 9) {  // px-d3z2-r2
        result = l * (nn - 0.5 * (ll + mm)) * pds - SQRT3 * l * nn * pdp;
    }
    else if (a == 3 && b == 5) {  // py-dxy
        result = SQRT3 * mm * l * pds + l * (1.0 - 2.0 * mm) * pdp;
    }
    else if (a == 3 && b == 6) {  // py-dyz
        result = SQRT3 * mm * n * pds + n * (1.0 - 2.0 * mm) * pdp;
    }
    else if (a == 3 && b == 7) {  // py-dzx
        result = SQRT3 * l * m * n * pds - 2.0 * l * m * n * pdp;
    }
    else if (a == 3 && b == 8) {  // py-dx2-y2
        result = 0.5 * SQRT3 * m * (ll - mm) * pds - m * (1.0 + ll - mm) * pdp;
    }
    else if (a == 3 && b == 9) {  // py-d3z2-r2
        result = m * (nn - 0.5 * (ll + mm)) * pds - SQRT3 * m * nn * pdp;
    }
    else if (a == 4 && b == 5) {  // pz-dxy
        result = SQRT3 * l * m * n * pds - 2.0 * l * m * n * pdp;
    }
    else if (a == 4 && b == 6) {  // pz-dyz
        result = SQRT3 * nn * m * pds + m * (1.0 - 2.0 * nn) * pdp;
    }
    else if (a == 4 && b == 7) {  // pz-dzx
        result = SQRT3 * nn * l * pds + l * (1.0 - 2.0 * nn) * pdp;
    }
    else if (a == 4 && b == 8) {  // pz-dx2-y2
        result = 0.5 * SQRT3 * n * (ll - mm) * pds - n * (ll - mm) * pdp;
    }
    else if (a == 4 && b == 9) {  // pz-d3z2-r2
        result = n * (nn - 0.5 * (ll + mm)) * pds + SQRT3 * n * (ll + mm) * pdp;
    }
    // d-d interactions (a=5,6,7,8,9, b=5,6,7,8,9)
    else if (a == 5 && b == 5) {  // dxy-dxy
        result = 3.0 * ll * mm * dds + (ll + mm - 4.0 * ll * mm) * ddp + (nn + ll * mm) * ddd;
    }
    else if (a == 5 && b == 6) {  // dxy-dyz
        result = 3.0 * l * mm * n * dds + l * n * (1.0 - 4.0 * mm) * ddp + l * n * (mm - 1.0) * ddd;
    }
    else if (a == 5 && b == 7) {  // dxy-dzx
        result = 3.0 * ll * m * n * dds + m * n * (1.0 - 4.0 * ll) * ddp + m * n * (ll - 1.0) * ddd;
    }
    else if (a == 5 && b == 8) {  // dxy-dx2-y2
        result = 1.5 * l * m * (ll - mm) * dds + 2.0 * l * m * (mm - ll) * ddp + 0.5 * l * m * (ll - mm) * ddd;
    }
    else if (a == 5 && b == 9) {  // dxy-d3z2-r2
        result = SQRT3 * l * m * (nn - 0.5 * (ll + mm)) * dds - 2.0 * SQRT3 * l * m * nn * ddp
                 + 0.5 * SQRT3 * l * m * (1.0 + nn) * ddd;
    }
    else if (a == 6 && b == 6) {  // dyz-dyz
        result = 3.0 * mm * nn * dds + (mm + nn - 4.0 * mm * nn) * ddp + (ll + mm * nn) * ddd;
    }
    else if (a == 6 && b == 7) {  // dyz-dzx
        result = 3.0 * l * m * nn * dds + l * m * (1.0 - 4.0 * nn) * ddp + l * m * (nn - 1.0) * ddd;
    }
    else if (a == 6 && b == 8) {  // dyz-dx2-y2
        result = 1.5 * m * n * (ll - mm) * dds - m * n * (1.0 + 2.0 * (ll - mm)) * ddp
                 + m * n * (1.0 + 0.5 * (ll - mm)) * ddd;
    }
    else if (a == 6 && b == 9) {  // dyz-d3z2-r2
        result = SQRT3 * m * n * (nn - 0.5 * (ll + mm)) * dds + SQRT3 * m * n * (ll + mm - nn) * ddp
                 - 0.5 * SQRT3 * m * n * (ll + mm) * ddd;
    }
    else if (a == 7 && b == 7) {  // dzx-dzx
        result = 3.0 * ll * nn * dds + (ll + nn - 4.0 * ll * nn) * ddp + (mm + ll * nn) * ddd;
    }
    else if (a == 7 && b == 8) {  // dzx-dx2-y2
        result = 1.5 * n * l * (ll - mm) * dds + n * l * (1.0 - 2.0 * (ll - mm)) * ddp
                 - n * l * (1.0 - 0.5 * (ll - mm)) * ddd;
    }
    else if (a == 7 && b == 9) {  // dzx-d3z2-r2
        result = SQRT3 * l * n * (nn - 0.5 * (ll + mm)) * dds + SQRT3 * l * n * (ll + mm - nn) * ddp
                 - 0.5 * SQRT3 * l * n * (ll + mm) * ddd;
    }
    else if (a == 8 && b == 8) {  // dx2-y2 - dx2-y2
        Scalar lm2 = ll - mm;
        result = 0.75 * lm2 * lm2 * dds + (ll + mm - lm2 * lm2) * ddp + (nn + 0.25 * lm2 * lm2) * ddd;
    }
    else if (a == 8 && b == 9) {  // dx2-y2 - d3z2-r2
        result = 0.5 * SQRT3 * (ll - mm) * (nn - 0.5 * (ll + mm)) * dds + SQRT3 * nn * (mm - ll) * ddp
                 + 0.25 * SQRT3 * (1.0 + nn) * (ll - mm) * ddd;
    }
    else if (a == 9 && b == 9) {  // d3z2-r2 - d3z2-r2
        Scalar nnh = nn - 0.5 * (ll + mm);
        result = nnh * nnh * dds + 3.0 * nn * (ll + mm) * ddp + 0.75 * (ll + mm) * (ll + mm) * ddd;
    }

    // Apply parity factor for swapped orbitals
    // Factor is (-1)^(l_a + l_b) where l is angular momentum
    if (swapped) {
        int la = ORBITAL_L[a - 1];
        int lb = ORBITAL_L[b - 1];
        if ((la + lb) % 2 == 1) {
            result = -result;
        }
    }

    return result;
}

/**
 * @brief Compute derivatives of direction cosines
 *
 * @param c Direction cosines [l, m, n]
 * @param r Distance
 * @return Array of direction cosine derivatives [dl/dx, dl/dy, dl/dz, dm/dx, ...]
 */
inline std::array<Scalar, 9> compute_dc_derivatives(const Vec3& c, Scalar r) {
    const Scalar l = c[0];
    const Scalar m = c[1];
    const Scalar n = c[2];
    const Scalar r_inv = 1.0 / r;

    // Derivatives of direction cosines: d(c_i)/d(x_j) = (delta_ij - c_i * c_j) / r
    std::array<Scalar, 9> dc;
    dc[0] = (1.0 - l * l) * r_inv;  // dl/dx
    dc[1] = -l * m * r_inv;         // dl/dy
    dc[2] = -l * n * r_inv;         // dl/dz
    dc[3] = -m * l * r_inv;         // dm/dx
    dc[4] = (1.0 - m * m) * r_inv;  // dm/dy
    dc[5] = -m * n * r_inv;         // dm/dz
    dc[6] = -n * l * r_inv;         // dn/dx
    dc[7] = -n * m * r_inv;         // dn/dy
    dc[8] = (1.0 - n * n) * r_inv;  // dn/dz

    return dc;
}

/**
 * @brief Compute spatial derivatives of SK-transformed matrix element
 *
 * This implements the chain rule: dH/dr_i = dH/dSK * dSK/dr * c_i + H_geom * dc/dr
 *
 * @param a First orbital index (1-9)
 * @param b Second orbital index (1-9)
 * @param c Direction cosines [l, m, n]
 * @param r Distance
 * @param sk SK integrals
 * @param dsk SK integral derivatives (dSK/dr)
 * @return Gradient of matrix element [dH/dx, dH/dy, dH/dz]
 */
inline Vec3 transform_orb_derivative(int a, int b, const Vec3& c, Scalar r,
                                     const std::array<Scalar, NUM_SK_INTEGRALS>& sk,
                                     const std::array<Scalar, NUM_SK_INTEGRALS>& dsk) {
    Vec3 gradient = Vec3::Zero();

    // Extract direction cosines
    const Scalar l = c[0];
    const Scalar m = c[1];
    const Scalar n = c[2];

    // Precompute powers
    const Scalar ll = l * l;
    const Scalar mm = m * m;
    const Scalar nn = n * n;

    // Get direction cosine derivatives
    auto dc = compute_dc_derivatives(c, r);
    const Scalar li_x = dc[0], li_y = dc[1], li_z = dc[2];
    const Scalar mi_x = dc[3], mi_y = dc[4], mi_z = dc[5];
    const Scalar ni_x = dc[6], ni_y = dc[7], ni_z = dc[8];

    // Handle parity correction
    bool swapped = false;
    if (a > b) {
        std::swap(a, b);
        swapped = true;
    }

    // Extract SK integrals and derivatives
    const Scalar dds = sk[0], ddp = sk[1], ddd = sk[2];
    const Scalar pds = sk[3], pdp = sk[4];
    const Scalar pps = sk[5], ppp = sk[6];
    const Scalar sds = sk[7], sps = sk[8], sss = sk[9];

    const Scalar d_dds = dsk[0], d_ddp = dsk[1], d_ddd = dsk[2];
    const Scalar d_pds = dsk[3], d_pdp = dsk[4];
    const Scalar d_pps = dsk[5], d_ppp = dsk[6];
    const Scalar d_sds = dsk[7], d_sps = dsk[8], d_sss = dsk[9];

    // For each orbital pair, compute both radial and geometric contributions
    // The full derivative is: dH/dx_k = (dH/dr) * c_k + H_geometric_deriv

    // Note: This is a simplified implementation. The full mdiff function in Fortran
    // is about 300 lines. Here we implement the most common cases.

    // s-s interaction
    if (a == 1 && b == 1) {
        gradient[0] = d_sss * l;
        gradient[1] = d_sss * m;
        gradient[2] = d_sss * n;
    }
    // s-p interactions
    else if (a == 1 && b == 2) {  // s-px
        gradient[0] = d_sps * l * l + sps * li_x;
        gradient[1] = d_sps * l * m + sps * li_y;
        gradient[2] = d_sps * l * n + sps * li_z;
    }
    else if (a == 1 && b == 3) {  // s-py
        gradient[0] = d_sps * m * l + sps * mi_x;
        gradient[1] = d_sps * m * m + sps * mi_y;
        gradient[2] = d_sps * m * n + sps * mi_z;
    }
    else if (a == 1 && b == 4) {  // s-pz
        gradient[0] = d_sps * n * l + sps * ni_x;
        gradient[1] = d_sps * n * m + sps * ni_y;
        gradient[2] = d_sps * n * n + sps * ni_z;
    }
    // p-p interactions
    else if (a == 2 && b == 2) {  // px-px
        Scalar geom = ll * pps + (1.0 - ll) * ppp;
        gradient[0] = (ll * d_pps + (1.0 - ll) * d_ppp) * l + 2.0 * l * li_x * (pps - ppp);
        gradient[1] = (ll * d_pps + (1.0 - ll) * d_ppp) * m + 2.0 * l * li_y * (pps - ppp);
        gradient[2] = (ll * d_pps + (1.0 - ll) * d_ppp) * n + 2.0 * l * li_z * (pps - ppp);
    }
    else if (a == 2 && b == 3) {  // px-py
        Scalar diff = pps - ppp;
        Scalar d_diff = d_pps - d_ppp;
        gradient[0] = d_diff * l * m * l + diff * (li_x * m + l * mi_x);
        gradient[1] = d_diff * l * m * m + diff * (li_y * m + l * mi_y);
        gradient[2] = d_diff * l * m * n + diff * (li_z * m + l * mi_z);
    }
    else if (a == 2 && b == 4) {  // px-pz
        Scalar diff = pps - ppp;
        Scalar d_diff = d_pps - d_ppp;
        gradient[0] = d_diff * l * n * l + diff * (li_x * n + l * ni_x);
        gradient[1] = d_diff * l * n * m + diff * (li_y * n + l * ni_y);
        gradient[2] = d_diff * l * n * n + diff * (li_z * n + l * ni_z);
    }
    else if (a == 3 && b == 3) {  // py-py
        gradient[0] = (mm * d_pps + (1.0 - mm) * d_ppp) * l + 2.0 * m * mi_x * (pps - ppp);
        gradient[1] = (mm * d_pps + (1.0 - mm) * d_ppp) * m + 2.0 * m * mi_y * (pps - ppp);
        gradient[2] = (mm * d_pps + (1.0 - mm) * d_ppp) * n + 2.0 * m * mi_z * (pps - ppp);
    }
    else if (a == 3 && b == 4) {  // py-pz
        Scalar diff = pps - ppp;
        Scalar d_diff = d_pps - d_ppp;
        gradient[0] = d_diff * m * n * l + diff * (mi_x * n + m * ni_x);
        gradient[1] = d_diff * m * n * m + diff * (mi_y * n + m * ni_y);
        gradient[2] = d_diff * m * n * n + diff * (mi_z * n + m * ni_z);
    }
    else if (a == 4 && b == 4) {  // pz-pz
        gradient[0] = (nn * d_pps + (1.0 - nn) * d_ppp) * l + 2.0 * n * ni_x * (pps - ppp);
        gradient[1] = (nn * d_pps + (1.0 - nn) * d_ppp) * m + 2.0 * n * ni_y * (pps - ppp);
        gradient[2] = (nn * d_pps + (1.0 - nn) * d_ppp) * n + 2.0 * n * ni_z * (pps - ppp);
    }
    // s-d interactions
    else if (a == 1 && b == 5) {  // s-dxy
        Scalar geom = SQRT3 * l * m;
        gradient[0] = SQRT3 * d_sds * l * m * l + SQRT3 * sds * (li_x * m + l * mi_x);
        gradient[1] = SQRT3 * d_sds * l * m * m + SQRT3 * sds * (li_y * m + l * mi_y);
        gradient[2] = SQRT3 * d_sds * l * m * n + SQRT3 * sds * (li_z * m + l * mi_z);
    }
    else if (a == 1 && b == 6) {  // s-dyz
        gradient[0] = SQRT3 * d_sds * m * n * l + SQRT3 * sds * (mi_x * n + m * ni_x);
        gradient[1] = SQRT3 * d_sds * m * n * m + SQRT3 * sds * (mi_y * n + m * ni_y);
        gradient[2] = SQRT3 * d_sds * m * n * n + SQRT3 * sds * (mi_z * n + m * ni_z);
    }
    else if (a == 1 && b == 7) {  // s-dzx
        gradient[0] = SQRT3 * d_sds * n * l * l + SQRT3 * sds * (ni_x * l + n * li_x);
        gradient[1] = SQRT3 * d_sds * n * l * m + SQRT3 * sds * (ni_y * l + n * li_y);
        gradient[2] = SQRT3 * d_sds * n * l * n + SQRT3 * sds * (ni_z * l + n * li_z);
    }
    else if (a == 1 && b == 8) {  // s-dx2-y2
        Scalar lm = ll - mm;
        gradient[0] = 0.5 * SQRT3 * d_sds * lm * l + SQRT3 * sds * (l * li_x - m * mi_x);
        gradient[1] = 0.5 * SQRT3 * d_sds * lm * m + SQRT3 * sds * (l * li_y - m * mi_y);
        gradient[2] = 0.5 * SQRT3 * d_sds * lm * n + SQRT3 * sds * (l * li_z - m * mi_z);
    }
    else if (a == 1 && b == 9) {  // s-d3z2-r2
        Scalar nnh = nn - 0.5 * (ll + mm);
        gradient[0] = d_sds * nnh * l + sds * (2.0 * n * ni_x - l * li_x - m * mi_x);
        gradient[1] = d_sds * nnh * m + sds * (2.0 * n * ni_y - l * li_y - m * mi_y);
        gradient[2] = d_sds * nnh * n + sds * (2.0 * n * ni_z - l * li_z - m * mi_z);
    }
    // For higher orbital combinations (p-d, d-d), we use numerical differentiation
    // in the actual implementation. Here we provide placeholder zeros.
    else {
        // For a complete implementation, these would need the full mdiff formulas
        // from the Fortran code. For now, use a simple finite difference approach
        // in the Hamiltonian class.
        gradient.setZero();
    }

    // Apply parity factor for swapped orbitals
    if (swapped) {
        int la = ORBITAL_L[a - 1];
        int lb = ORBITAL_L[b - 1];
        if ((la + lb) % 2 == 1) {
            gradient = -gradient;
        }
    }

    return gradient;
}

/**
 * @brief Get the required SK integral indices for a given orbital configuration
 *
 * @param no1 Number of orbitals on atom 1
 * @param no2 Number of orbitals on atom 2
 * @return Vector of SK integral indices needed
 */
inline std::vector<int> get_required_integrals(int no1, int no2) {
    std::vector<int> result;

    // Determine the combined orbital configuration
    int l_max1 = (no1 == 1) ? 0 : ((no1 == 4) ? 1 : 2);
    int l_max2 = (no2 == 1) ? 0 : ((no2 == 4) ? 1 : 2);

    // s-s always needed if both have s
    if (l_max1 >= 0 && l_max2 >= 0) result.push_back(9);  // sss

    // s-p needed if one has s and other has p
    if ((l_max1 >= 0 && l_max2 >= 1) || (l_max1 >= 1 && l_max2 >= 0)) result.push_back(8);  // sps

    // p-p needed if both have p
    if (l_max1 >= 1 && l_max2 >= 1) {
        result.push_back(5);  // pps
        result.push_back(6);  // ppp
    }

    // s-d needed if one has s and other has d
    if ((l_max1 >= 0 && l_max2 >= 2) || (l_max1 >= 2 && l_max2 >= 0)) result.push_back(7);  // sds

    // p-d needed if one has p and other has d
    if ((l_max1 >= 1 && l_max2 >= 2) || (l_max1 >= 2 && l_max2 >= 1)) {
        result.push_back(3);  // pds
        result.push_back(4);  // pdp
    }

    // d-d needed if both have d
    if (l_max1 >= 2 && l_max2 >= 2) {
        result.push_back(0);  // dds
        result.push_back(1);  // ddp
        result.push_back(2);  // ddd
    }

    return result;
}

/**
 * @brief Map orbital index from reduced to full basis
 *
 * For elements that don't have all orbitals, this maps the condensed
 * orbital index to the absolute orbital index.
 *
 * @param no Number of orbitals (1, 4, 5, 6, 8, or 9)
 * @param a0 Input orbital index (1-based, within reduced basis)
 * @return Absolute orbital index (1-9)
 */
inline int get_absolute_orbital(int no, int a0) {
    int a = a0;
    if (no == 5) a = a + 4;  // d orbitals only (5->9)
    if (no == 8 && a > 0) a = a + 1;  // pd orbitals (skip s)
    if (no == 6 && a > 1) a = a + 3;  // sd orbitals
    return a;
}

} // namespace tb
} // namespace atomistica
