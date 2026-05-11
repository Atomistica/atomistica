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
 * This implements the complete Slater-Koster derivative transformation,
 * following the mdiff function from dense_forces.f90. The derivative has
 * two contributions:
 * - d: radial part where SK integrals are differentiated
 * - g: geometric part where direction cosines are differentiated
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

    // Handle parity correction for swapped orbitals
    bool swapped = false;
    if (a > b) {
        std::swap(a, b);
        swapped = true;
    }

    // Extract SK integrals
    const Scalar dds = sk[0], ddp = sk[1], ddd = sk[2];
    const Scalar pds = sk[3], pdp = sk[4];
    const Scalar pps = sk[5], ppp = sk[6];
    const Scalar sds = sk[7], sps = sk[8], sss = sk[9];

    // Loop over x, y, z components
    for (int i = 0; i < 3; ++i) {
        // Compute radial derivatives: d(sk)/dr * c[i]
        const Scalar ddsi = dsk[0] * c[i];
        const Scalar ddpi = dsk[1] * c[i];
        const Scalar dddi = dsk[2] * c[i];
        const Scalar pdsi = dsk[3] * c[i];
        const Scalar pdpi = dsk[4] * c[i];
        const Scalar ppsi = dsk[5] * c[i];
        const Scalar pppi = dsk[6] * c[i];
        const Scalar sdsi = dsk[7] * c[i];
        const Scalar spsi = dsk[8] * c[i];
        const Scalar sssi = dsk[9] * c[i];

        // Direction cosine derivatives: d(c_j)/d(x_i) = (delta_ij - c_i * c_j) / r
        const Scalar li = ((i == 0 ? 1.0 : 0.0) - l * c[i]) / r;
        const Scalar mi = ((i == 1 ? 1.0 : 0.0) - m * c[i]) / r;
        const Scalar ni = ((i == 2 ? 1.0 : 0.0) - n * c[i]) / r;
        const Scalar lli = 2.0 * l * li;
        const Scalar mmi = 2.0 * m * mi;
        const Scalar nni = 2.0 * n * ni;

        Scalar d = 0.0;  // Radial contribution
        Scalar g = 0.0;  // Geometric contribution

        // Select transformation rule based on orbital pair (a <= b)
        if (a == 1) {  // s orbital
            if (b == 1) {  // s-s
                d = sssi;
                g = 0.0;
            } else if (b == 2) {  // s-px
                d = l * spsi;
                g = li * sps;
            } else if (b == 3) {  // s-py
                d = m * spsi;
                g = mi * sps;
            } else if (b == 4) {  // s-pz
                d = n * spsi;
                g = ni * sps;
            } else if (b == 5) {  // s-dxy
                d = SQRT3 * l * m * sdsi;
                g = SQRT3 * (li * m + l * mi) * sds;
            } else if (b == 6) {  // s-dyz
                d = SQRT3 * m * n * sdsi;
                g = SQRT3 * (mi * n + m * ni) * sds;
            } else if (b == 7) {  // s-dzx
                d = SQRT3 * n * l * sdsi;
                g = SQRT3 * (ni * l + n * li) * sds;
            } else if (b == 8) {  // s-dx2-y2
                d = 0.5 * SQRT3 * (ll - mm) * sdsi;
                g = 0.5 * SQRT3 * (lli - mmi) * sds;
            } else if (b == 9) {  // s-d3z2-r2
                d = (nn - 0.5 * (ll + mm)) * sdsi;
                g = (nni - 0.5 * (lli + mmi)) * sds;
            }
        } else if (a == 2) {  // px orbital
            if (b == 2) {  // px-px
                d = ll * ppsi + (1.0 - ll) * pppi;
                g = lli * pps + (-lli) * ppp;
            } else if (b == 3) {  // px-py
                d = l * m * ppsi - l * m * pppi;
                g = (li * m + l * mi) * pps - (li * m + l * mi) * ppp;
            } else if (b == 4) {  // px-pz
                d = l * n * ppsi - l * n * pppi;
                g = (li * n + l * ni) * pps - (li * n + l * ni) * ppp;
            } else if (b == 5) {  // px-dxy
                d = SQRT3 * ll * m * pdsi + m * (1.0 - 2.0 * ll) * pdpi;
                g = SQRT3 * (lli * m + ll * mi) * pds + (mi * (1.0 - 2.0 * ll) + m * (-2.0 * lli)) * pdp;
            } else if (b == 6) {  // px-dyz
                d = SQRT3 * l * m * n * pdsi - 2.0 * l * m * n * pdpi;
                g = SQRT3 * (li * m * n + l * mi * n + l * m * ni) * pds - 2.0 * (li * m * n + l * mi * n + l * m * ni) * pdp;
            } else if (b == 7) {  // px-dzx
                d = SQRT3 * ll * n * pdsi + n * (1.0 - 2.0 * ll) * pdpi;
                g = SQRT3 * (lli * n + ll * ni) * pds + (ni * (1.0 - 2.0 * ll) + n * (-2.0 * lli)) * pdp;
            } else if (b == 8) {  // px-dx2-y2
                d = 0.5 * SQRT3 * l * (ll - mm) * pdsi + l * (1.0 - ll + mm) * pdpi;
                g = 0.5 * SQRT3 * (li * (ll - mm) + l * (lli - mmi)) * pds + (li * (1.0 - ll + mm) + l * (-lli + mmi)) * pdp;
            } else if (b == 9) {  // px-d3z2-r2
                d = l * (nn - 0.5 * (ll + mm)) * pdsi - SQRT3 * l * nn * pdpi;
                g = (li * (nn - 0.5 * (ll + mm)) + l * (nni - 0.5 * (lli + mmi))) * pds - SQRT3 * (li * nn + l * nni) * pdp;
            }
        } else if (a == 3) {  // py orbital
            if (b == 3) {  // py-py
                d = mm * ppsi + (1.0 - mm) * pppi;
                g = mmi * pps + (-mmi) * ppp;
            } else if (b == 4) {  // py-pz
                d = m * n * ppsi - m * n * pppi;
                g = (mi * n + m * ni) * pps - (mi * n + m * ni) * ppp;
            } else if (b == 5) {  // py-dxy
                d = SQRT3 * mm * l * pdsi + l * (1.0 - 2.0 * mm) * pdpi;
                g = SQRT3 * (mmi * l + mm * li) * pds + (li * (1.0 - 2.0 * mm) + l * (-2.0 * mmi)) * pdp;
            } else if (b == 6) {  // py-dyz
                d = SQRT3 * mm * n * pdsi + n * (1.0 - 2.0 * mm) * pdpi;
                g = SQRT3 * (mmi * n + mm * ni) * pds + (ni * (1.0 - 2.0 * mm) + n * (-2.0 * mmi)) * pdp;
            } else if (b == 7) {  // py-dzx
                d = SQRT3 * m * n * l * pdsi - 2.0 * m * n * l * pdpi;
                g = SQRT3 * (mi * n * l + m * ni * l + m * n * li) * pds - 2.0 * (mi * n * l + m * ni * l + m * n * li) * pdp;
            } else if (b == 8) {  // py-dx2-y2
                d = 0.5 * SQRT3 * m * (ll - mm) * pdsi - m * (1.0 + ll - mm) * pdpi;
                g = 0.5 * SQRT3 * (mi * (ll - mm) + m * (lli - mmi)) * pds - (mi * (1.0 + ll - mm) + m * (lli - mmi)) * pdp;
            } else if (b == 9) {  // py-d3z2-r2
                d = m * (nn - 0.5 * (ll + mm)) * pdsi - SQRT3 * m * nn * pdpi;
                g = (mi * (nn - 0.5 * (ll + mm)) + m * (nni - 0.5 * (lli + mmi))) * pds - SQRT3 * (mi * nn + m * nni) * pdp;
            }
        } else if (a == 4) {  // pz orbital
            if (b == 4) {  // pz-pz
                d = nn * ppsi + (1.0 - nn) * pppi;
                g = nni * pps + (-nni) * ppp;
            } else if (b == 5) {  // pz-dxy
                d = SQRT3 * l * m * n * pdsi - 2.0 * m * n * l * pdpi;
                g = SQRT3 * (li * m * n + l * mi * n + l * m * ni) * pds - 2.0 * (mi * n * l + m * ni * l + m * n * li) * pdp;
            } else if (b == 6) {  // pz-dyz
                d = SQRT3 * nn * m * pdsi + m * (1.0 - 2.0 * nn) * pdpi;
                g = SQRT3 * (nni * m + nn * mi) * pds + (mi * (1.0 - 2.0 * nn) + m * (-2.0 * nni)) * pdp;
            } else if (b == 7) {  // pz-dzx
                d = SQRT3 * nn * l * pdsi + l * (1.0 - 2.0 * nn) * pdpi;
                g = SQRT3 * (nni * l + nn * li) * pds + (li * (1.0 - 2.0 * nn) + l * (-2.0 * nni)) * pdp;
            } else if (b == 8) {  // pz-dx2-y2
                d = 0.5 * SQRT3 * n * (ll - mm) * pdsi - n * (ll - mm) * pdpi;
                g = 0.5 * SQRT3 * (ni * (ll - mm) + n * (lli - mmi)) * pds - (ni * (ll - mm) + n * (lli - mmi)) * pdp;
            } else if (b == 9) {  // pz-d3z2-r2
                d = n * (nn - 0.5 * (ll + mm)) * pdsi + SQRT3 * n * (ll + mm) * pdpi;
                g = (ni * (nn - 0.5 * (ll + mm)) + n * (nni - 0.5 * (lli + mmi))) * pds + SQRT3 * (ni * (ll + mm) + n * (lli + mmi)) * pdp;
            }
        } else if (a == 5) {  // dxy orbital
            if (b == 5) {  // dxy-dxy
                d = 3.0 * ll * mm * ddsi + (ll + mm - 4.0 * ll * mm) * ddpi + (nn + ll * mm) * dddi;
                g = 3.0 * (lli * mm + ll * mmi) * dds + (lli + mmi - 4.0 * (lli * mm + ll * mmi)) * ddp + (nni + (lli * mm + ll * mmi)) * ddd;
            } else if (b == 6) {  // dxy-dyz
                d = 3.0 * l * mm * n * ddsi + l * n * (1.0 - 4.0 * mm) * ddpi + l * n * (mm - 1.0) * dddi;
                g = 3.0 * (li * mm * n + l * mmi * n + l * mm * ni) * dds
                    + (li * n * (1.0 - 4.0 * mm) + l * ni * (1.0 - 4.0 * mm) + l * n * (-4.0 * mmi)) * ddp
                    + (li * n * (mm - 1.0) + l * ni * (mm - 1.0) + l * n * mmi) * ddd;
            } else if (b == 7) {  // dxy-dzx
                d = 3.0 * ll * m * n * ddsi + m * n * (1.0 - 4.0 * ll) * ddpi + m * n * (ll - 1.0) * dddi;
                g = 3.0 * (lli * m * n + ll * mi * n + ll * m * ni) * dds
                    + (mi * n * (1.0 - 4.0 * ll) + m * ni * (1.0 - 4.0 * ll) + m * n * (-4.0 * lli)) * ddp
                    + (mi * n * (ll - 1.0) + m * ni * (ll - 1.0) + m * n * lli) * ddd;
            } else if (b == 8) {  // dxy-dx2-y2
                d = 1.5 * l * m * (ll - mm) * ddsi + 2.0 * l * m * (mm - ll) * ddpi + 0.5 * l * m * (ll - mm) * dddi;
                g = 1.5 * (li * m * (ll - mm) + l * mi * (ll - mm) + l * m * (lli - mmi)) * dds
                    + 2.0 * (li * m * (mm - ll) + l * mi * (mm - ll) + l * m * (mmi - lli)) * ddp
                    + 0.5 * (li * m * (ll - mm) + l * mi * (ll - mm) + l * m * (lli - mmi)) * ddd;
            } else if (b == 9) {  // dxy-d3z2-r2
                d = SQRT3 * l * m * (nn - 0.5 * (ll + mm)) * ddsi - 2.0 * SQRT3 * l * m * nn * ddpi + 0.5 * SQRT3 * l * m * (1.0 + nn) * dddi;
                g = SQRT3 * (li * m * (nn - 0.5 * (ll + mm)) + l * mi * (nn - 0.5 * (ll + mm)) + l * m * (nni - 0.5 * (lli + mmi))) * dds
                    - 2.0 * SQRT3 * (li * m * nn + l * mi * nn + l * m * nni) * ddp
                    + 0.5 * SQRT3 * (li * m * (1.0 + nn) + l * mi * (1.0 + nn) + l * m * nni) * ddd;
            }
        } else if (a == 6) {  // dyz orbital
            if (b == 6) {  // dyz-dyz
                d = 3.0 * mm * nn * ddsi + (mm + nn - 4.0 * mm * nn) * ddpi + (ll + mm * nn) * dddi;
                g = 3.0 * (mmi * nn + mm * nni) * dds + (mmi + nni - 4.0 * (mmi * nn + mm * nni)) * ddp + (lli + mmi * nn + mm * nni) * ddd;
            } else if (b == 7) {  // dyz-dzx
                d = 3.0 * m * nn * l * ddsi + m * l * (1.0 - 4.0 * nn) * ddpi + m * l * (nn - 1.0) * dddi;
                g = 3.0 * (mi * nn * l + m * nni * l + m * nn * li) * dds
                    + (mi * l * (1.0 - 4.0 * nn) + m * li * (1.0 - 4.0 * nn) + m * l * (-4.0 * nni)) * ddp
                    + (mi * l * (nn - 1.0) + m * li * (nn - 1.0) + m * l * nni) * ddd;
            } else if (b == 8) {  // dyz-dx2-y2
                d = 1.5 * m * n * (ll - mm) * ddsi - m * n * (1.0 + 2.0 * (ll - mm)) * ddpi + m * n * (1.0 + 0.5 * (ll - mm)) * dddi;
                g = 1.5 * (mi * n * (ll - mm) + m * ni * (ll - mm) + m * n * (lli - mmi)) * dds
                    - (mi * n * (1.0 + 2.0 * (ll - mm)) + m * ni * (1.0 + 2.0 * (ll - mm)) + m * n * (2.0 * lli - 2.0 * mmi)) * ddp
                    + (mi * n * (1.0 + 0.5 * (ll - mm)) + m * ni * (1.0 + 0.5 * (ll - mm)) + m * n * (0.5 * (lli - mmi))) * ddd;
            } else if (b == 9) {  // dyz-d3z2-r2
                d = SQRT3 * m * n * (nn - 0.5 * (ll + mm)) * ddsi + SQRT3 * m * n * (ll + mm - nn) * ddpi - 0.5 * SQRT3 * m * n * (ll + mm) * dddi;
                g = SQRT3 * (mi * n * (nn - 0.5 * (ll + mm)) + m * ni * (nn - 0.5 * (ll + mm)) + m * n * (nni - 0.5 * (lli + mmi))) * dds
                    + SQRT3 * (mi * n * (ll + mm - nn) + m * ni * (ll + mm - nn) + m * n * (lli + mmi - nni)) * ddp
                    - 0.5 * SQRT3 * (mi * n * (ll + mm) + m * ni * (ll + mm) + m * n * (lli + mmi)) * ddd;
            }
        } else if (a == 7) {  // dzx orbital
            if (b == 7) {  // dzx-dzx
                d = 3.0 * nn * ll * ddsi + (nn + ll - 4.0 * nn * ll) * ddpi + (mm + nn * ll) * dddi;
                g = 3.0 * (nni * ll + nn * lli) * dds + (nni + lli - 4.0 * (nni * ll + nn * lli)) * ddp + (mmi + nni * ll + nn * lli) * ddd;
            } else if (b == 8) {  // dzx-dx2-y2
                d = 1.5 * n * l * (ll - mm) * ddsi + n * l * (1.0 - 2.0 * (ll - mm)) * ddpi - n * l * (1.0 - 0.5 * (ll - mm)) * dddi;
                g = 1.5 * (ni * l * (ll - mm) + n * li * (ll - mm) + n * l * (lli - mmi)) * dds
                    + (ni * l * (1.0 - 2.0 * (ll - mm)) + n * li * (1.0 - 2.0 * (ll - mm)) + n * l * (-2.0 * (lli - mmi))) * ddp
                    - (ni * l * (1.0 - 0.5 * (ll - mm)) + n * li * (1.0 - 0.5 * (ll - mm)) + n * l * (-0.5 * (lli - mmi))) * ddd;
            } else if (b == 9) {  // dzx-d3z2-r2
                d = SQRT3 * l * n * (nn - 0.5 * (ll + mm)) * ddsi + SQRT3 * l * n * (ll + mm - nn) * ddpi - 0.5 * SQRT3 * l * n * (ll + mm) * dddi;
                g = SQRT3 * (li * n * (nn - 0.5 * (ll + mm)) + l * ni * (nn - 0.5 * (ll + mm)) + l * n * (nni - 0.5 * (lli + mmi))) * dds
                    + SQRT3 * (li * n * (ll + mm - nn) + l * ni * (ll + mm - nn) + l * n * (lli + mmi - nni)) * ddp
                    - 0.5 * SQRT3 * (li * n * (ll + mm) + l * ni * (ll + mm) + l * n * (lli + mmi)) * ddd;
            }
        } else if (a == 8) {  // dx2-y2 orbital
            if (b == 8) {  // dx2-y2 - dx2-y2
                Scalar lm2 = ll - mm;
                d = 0.75 * lm2 * lm2 * ddsi + (ll + mm - lm2 * lm2) * ddpi + (nn + 0.25 * lm2 * lm2) * dddi;
                g = 0.75 * 2.0 * lm2 * (lli - mmi) * dds + (lli + mmi - 2.0 * lm2 * (lli - mmi)) * ddp + (nni + 0.25 * 2.0 * lm2 * (lli - mmi)) * ddd;
            } else if (b == 9) {  // dx2-y2 - d3z2-r2
                d = 0.5 * SQRT3 * (ll - mm) * (nn - 0.5 * (ll + mm)) * ddsi + SQRT3 * nn * (mm - ll) * ddpi + 0.25 * SQRT3 * (1.0 + nn) * (ll - mm) * dddi;
                g = 0.5 * SQRT3 * ((lli - mmi) * (nn - 0.5 * (ll + mm)) + (ll - mm) * (nni - 0.5 * (lli + mmi))) * dds
                    + SQRT3 * (nni * (mm - ll) + nn * (mmi - lli)) * ddp
                    + 0.25 * SQRT3 * (nni * (ll - mm) + (1.0 + nn) * (lli - mmi)) * ddd;
            }
        } else if (a == 9) {  // d3z2-r2 orbital
            if (b == 9) {  // d3z2-r2 - d3z2-r2
                Scalar nnh = nn - 0.5 * (ll + mm);
                d = nnh * nnh * ddsi + 3.0 * nn * (ll + mm) * ddpi + 0.75 * (ll + mm) * (ll + mm) * dddi;
                g = 2.0 * nnh * (nni - 0.5 * (lli + mmi)) * dds + 3.0 * (nni * (ll + mm) + nn * (lli + mmi)) * ddp + 0.75 * 2.0 * (ll + mm) * (lli + mmi) * ddd;
            }
        }

        gradient[i] = d + g;
    }

    // Apply parity factor for swapped orbitals: (-1)^(l_a + l_b)
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
