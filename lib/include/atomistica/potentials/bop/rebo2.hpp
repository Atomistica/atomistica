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
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../potential_base.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../../math/cutoff_functions.hpp"

namespace atomistica {

/**
 * @brief 2nd generation REBO potential (Brenner 2002)
 *
 * This implements the Reactive Empirical Bond Order potential for
 * hydrocarbons as described in:
 * D.W. Brenner et al., J. Phys.: Condens. Matter 14 (2002) 783-802
 *
 * The potential has the form:
 *   E = sum_{i} sum_{j>i} [V_R(r_ij) - b_ij * V_A(r_ij)] * f_c(r_ij)
 *
 * where b_ij is a complex bond-order term that depends on:
 * - Angular contributions from neighbors
 * - Coordination numbers (N_C and N_H for each atom)
 * - Lookup tables F, P, T for conjugation corrections
 */

// Pair type indices (following Fortran convention)
constexpr int REBO2_C_C = 1;
constexpr int REBO2_C_H = 3;
constexpr int REBO2_H_H = 6;

// Internal element types
constexpr int REBO2_C = 1;  // Carbon internal type
constexpr int REBO2_H = 3;  // Hydrogen internal type

/**
 * @brief Coefficients for piecewise polynomial g(cos_theta) spline
 */
struct GSplineCoeffs {
    std::array<std::array<Scalar, 6>, 3> c;  // c[interval][power]

    GSplineCoeffs() {
        for (auto& interval : c) {
            interval.fill(0.0);
        }
    }
};

/**
 * @brief 2D lookup table with bicubic interpolation
 *
 * Implements bicubic spline interpolation for 2D tables.
 * Used for P tables in REBO2.
 */
class Table2D {
public:
    static constexpr int NPARA = 16;  // 4^2 coefficients per box
    static constexpr int NCORN = 4;   // 2^2 corners per box

    Table2D() = default;

    /**
     * @brief Initialize the table with values and optional derivatives
     */
    void init(int nx, int ny,
              const std::vector<std::vector<Scalar>>& values,
              const std::vector<std::vector<Scalar>>& dvdx = {},
              const std::vector<std::vector<Scalar>>& dvdy = {});

    /**
     * @brief Evaluate the table and its derivatives
     * @return tuple of (value, dv/dx, dv/dy)
     */
    std::tuple<Scalar, Scalar, Scalar> eval(Scalar x, Scalar y) const;

    bool is_valid() const { return !coeff_.empty(); }

    int nx() const { return nx_; }
    int ny() const { return ny_; }

private:
    std::vector<std::array<std::array<Scalar, 4>, 4>> coeff_;  // [box][i][j]
    int nx_ = 0, ny_ = 0;
    int nboxs_ = 0;
};

/**
 * @brief 3D lookup table with tricubic interpolation
 *
 * Implements tricubic spline interpolation for 3D tables.
 * Used for F and T tables in REBO2.
 */
class Table3D {
public:
    static constexpr int NPARA = 64;  // 4^3 coefficients per box
    static constexpr int NCORN = 8;   // 2^3 corners per box

    Table3D() = default;

    /**
     * @brief Initialize the table with values and optional derivatives
     */
    void init(int nx, int ny, int nz,
              const std::vector<std::vector<std::vector<Scalar>>>& values,
              const std::vector<std::vector<std::vector<Scalar>>>& dvdx = {},
              const std::vector<std::vector<std::vector<Scalar>>>& dvdy = {},
              const std::vector<std::vector<std::vector<Scalar>>>& dvdz = {});

    /**
     * @brief Evaluate the table and its derivatives
     * @return tuple of (value, dv/dx, dv/dy, dv/dz)
     */
    std::tuple<Scalar, Scalar, Scalar, Scalar> eval(Scalar x, Scalar y, Scalar z) const;

    bool is_valid() const { return !coeff_.empty(); }

    int nx() const { return nx_; }
    int ny() const { return ny_; }
    int nz() const { return nz_; }

private:
    std::vector<std::array<std::array<std::array<Scalar, 4>, 4>, 4>> coeff_;  // [box][i][j][k]
    int nx_ = 0, ny_ = 0, nz_ = 0;
    int nboxs_ = 0;
};

/**
 * @brief REBO2 potential implementation
 */
class REBO2 {
public:
    REBO2();

    /**
     * @brief Load default parameters (Brenner 2002)
     */
    void load_default_parameters();

    /**
     * @brief Get cutoff radius
     */
    Scalar cutoff() const { return max_cutoff_; }

    /**
     * @brief Get element index from atomic number
     * @return Internal element type (REBO2_C or REBO2_H), or -1 if not supported
     */
    int element_type(int Z) const;

    /**
     * @brief Get pair type from two element types
     */
    static int pair_type(int eli, int elj);

    /**
     * @brief Compute energy, forces, and virial
     */
    PotentialResults compute(AtomicSystem& system, NeighborList& neighbors,
                            bool compute_forces = true, bool compute_virial = true);

    // === Pair functions (exposed for testing) ===

    /**
     * @brief Repulsive potential V_R(r)
     *
     * C-C: V_R = (1 + Q/r) * A * exp(-alpha * r)
     * C-H, H-H: Same form with different parameters
     */
    std::pair<Scalar, Scalar> repulsive(int ptype, Scalar r) const;

    /**
     * @brief Attractive potential V_A(r)
     *
     * C-C: V_A = -sum_{i=1}^3 B_i * exp(-beta_i * r)
     * C-H, H-H: V_A = -B * exp(-beta * r)
     */
    std::pair<Scalar, Scalar> attractive(int ptype, Scalar r) const;

    /**
     * @brief Angular function g(cos_theta, N)
     *
     * For C: Uses piecewise polynomial spline, interpolates between
     *        g1 (N > 3.7) and g2 (N < 3.2) based on coordination
     * For H: Uses polynomial fit
     */
    std::tuple<Scalar, Scalar, Scalar> angular_function(
        int el_type, Scalar cos_theta, Scalar N) const;

    /**
     * @brief Distance weighting function h(r_ij - r_ik)
     *
     * Used to weight bond-order contributions from different bond types
     */
    std::pair<Scalar, Scalar> distance_weight(int ptype_ij, int ptype_ik, Scalar dr) const;

    /**
     * @brief Bond order function b(1 + z)
     *
     * b = (1 + z)^(-0.5)
     */
    std::pair<Scalar, Scalar> bond_order_func(int el_type, Scalar z) const;

    // === Parameters (public for testing/modification) ===

    // C-C attractive parameters (3 exponential terms)
    Scalar cc_B1 = 12388.79197798;
    Scalar cc_B2 = 17.56740646509;
    Scalar cc_B3 = 30.71493208065;
    Scalar cc_beta1 = 4.7204523127;
    Scalar cc_beta2 = 1.4332132499;
    Scalar cc_beta3 = 1.3826912506;

    // C-C repulsive parameters
    Scalar cc_Q = 0.3134602960833;
    Scalar cc_A = 10953.544162170;
    Scalar cc_alpha = 4.7465390606595;

    // C-H parameters
    Scalar ch_B1 = 32.3551866587;
    Scalar ch_beta1 = 1.43445805925;
    Scalar ch_Q = 0.340775728;
    Scalar ch_A = 149.94098723;
    Scalar ch_alpha = 4.10254983;

    // H-H parameters
    Scalar hh_B1 = 29.632593;
    Scalar hh_beta1 = 1.71589217;
    Scalar hh_Q = 0.370471487045;
    Scalar hh_A = 32.817355747;
    Scalar hh_alpha = 3.536298648;

    // Cutoff radii
    Scalar cc_r1 = 1.70, cc_r2 = 2.00;
    Scalar ch_r1 = 1.30, ch_r2 = 1.80;
    Scalar hh_r1 = 1.10, hh_r2 = 1.70;

    // Equilibrium distances
    Scalar cc_re = 1.4;
    Scalar ch_re = 1.09;
    Scalar hh_re = 0.7415886997;

    // Lambda for distance weighting
    Scalar lambda = 4.0;

    // Angular function parameters for C-C
    std::array<Scalar, 6> cc_g_theta = {-1.0, -0.5, -1.0/3.0, 0.0, 0.5, 1.0};
    std::array<Scalar, 6> cc_g_g1 = {-0.01, 0.05280, 0.09733, 0.37545, 2.0014, 8.0};
    std::array<Scalar, 6> cc_g_dg1 = {0.10400, 0.17000, 0.40000, 0.0, 0.0, 0.0};
    std::array<Scalar, 6> cc_g_d2g1 = {0.00000, 0.37000, 1.98000, 0.0, 0.0, 0.0};
    std::array<Scalar, 6> cc_g_g2 = {0.0, 0.0, 0.09733, 0.271856, 0.416335, 1.0};

    // Angular function parameters for H-H (polynomial coefficients)
    std::array<std::array<Scalar, 6>, 3> hh_g_spline;
    std::array<int, 25> hh_g_intervals;

    // Enable dihedral terms
    bool with_dihedral = false;

protected:
    void init_cutoffs();
    void init_angular_splines();
    void init_tables();
    void init_distance_weights();

    void make_cc_g_spline(GSplineCoeffs& coeffs,
                          const std::array<Scalar, 6>& g_vals,
                          const std::array<Scalar, 6>& dg_vals,
                          const std::array<Scalar, 6>& d2g_vals);

    Scalar eval_cc_g_spline(const GSplineCoeffs& coeffs, Scalar cos_theta) const;
    std::pair<Scalar, Scalar> eval_cc_g_spline_with_deriv(
        const GSplineCoeffs& coeffs, Scalar cos_theta) const;

    // Cutoff functions
    TrigOffCutoff cc_cutoff_, ch_cutoff_, hh_cutoff_;

    // Angular spline coefficients
    GSplineCoeffs cc_g1_coeff_, cc_g2_coeff_;

    // Distance weight prefactors: exp(lambda * (r_ref_ij - r_ref_ik))
    std::array<std::array<Scalar, 7>, 7> conear_;

    // Bond order exponents
    Scalar conpe_C_ = -0.5;  // For C atoms
    Scalar conpe_H_ = -0.5;  // For H atoms

    // Lookup tables
    Table3D Fcc_, Fch_, Fhh_;
    Table2D Pcc_, Pch_;
    Table3D Tcc_;

    Scalar max_cutoff_ = 0.0;
    bool initialized_ = false;
};

// =========================================================================
// Table2D Implementation
// =========================================================================

inline void Table2D::init(int nx, int ny,
                          const std::vector<std::vector<Scalar>>& values,
                          const std::vector<std::vector<Scalar>>& dvdx,
                          const std::vector<std::vector<Scalar>>& dvdy) {
    nx_ = nx;
    ny_ = ny;
    nboxs_ = nx * ny;

    coeff_.resize(nboxs_);

    // Corner indices for normalized coordinates
    constexpr int ix1[NCORN] = {0, 1, 1, 0};
    constexpr int ix2[NCORN] = {0, 0, 1, 1};

    // Build matrix A (same for all boxes in normalized coordinates)
    std::array<std::array<Scalar, NPARA>, NPARA> A{};

    for (int icorn = 0; icorn < NCORN; ++icorn) {
        int nx1 = ix1[icorn];
        int nx2 = ix2[icorn];

        for (int npow1 = 0; npow1 < 4; ++npow1) {
            for (int npow2 = 0; npow2 < 4; ++npow2) {
                int npow1m = (npow1 > 0) ? npow1 - 1 : 0;
                int npow2m = (npow2 > 0) ? npow2 - 1 : 0;
                int icol = 4 * npow1 + npow2;

                // Function value
                A[icorn][icol] = std::pow(nx1, npow1) * std::pow(nx2, npow2);
                // d/dx
                A[icorn + NCORN][icol] = npow1 * std::pow(nx1, npow1m) * std::pow(nx2, npow2);
                // d/dy
                A[icorn + 2*NCORN][icol] = std::pow(nx1, npow1) * npow2 * std::pow(nx2, npow2m);
                // d2/dxdy
                A[icorn + 3*NCORN][icol] = npow1 * std::pow(nx1, npow1m) * npow2 * std::pow(nx2, npow2m);
            }
        }
    }

    // RHS vectors for each box
    std::vector<std::array<Scalar, NPARA>> B(nboxs_);

    for (int nhbox = 0; nhbox < nx; ++nhbox) {
        for (int ncbox = 0; ncbox < ny; ++ncbox) {
            int ibox = ny * nhbox + ncbox;

            for (int icorn = 0; icorn < NCORN; ++icorn) {
                int nx1 = ix1[icorn] + nhbox;
                int nx2 = ix2[icorn] + ncbox;

                B[ibox][icorn] = values[nx1][nx2];
                B[ibox][icorn + NCORN] = dvdx.empty() ? 0.0 : dvdx[nx1][nx2];
                B[ibox][icorn + 2*NCORN] = dvdy.empty() ? 0.0 : dvdy[nx1][nx2];
                B[ibox][icorn + 3*NCORN] = 0.0;  // Cross derivative assumed zero
            }
        }
    }

    // Solve A * x = B for each box using Gaussian elimination
    // First, LU decompose A
    std::array<std::array<Scalar, NPARA>, NPARA> LU = A;
    std::array<int, NPARA> perm;
    for (int i = 0; i < NPARA; ++i) perm[i] = i;

    for (int k = 0; k < NPARA; ++k) {
        // Find pivot
        int maxrow = k;
        Scalar maxval = std::abs(LU[k][k]);
        for (int i = k + 1; i < NPARA; ++i) {
            if (std::abs(LU[i][k]) > maxval) {
                maxval = std::abs(LU[i][k]);
                maxrow = i;
            }
        }
        if (maxrow != k) {
            std::swap(LU[k], LU[maxrow]);
            std::swap(perm[k], perm[maxrow]);
        }

        if (std::abs(LU[k][k]) < 1e-15) {
            throw std::runtime_error("Table2D: Singular matrix in init");
        }

        for (int i = k + 1; i < NPARA; ++i) {
            LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < NPARA; ++j) {
                LU[i][j] -= LU[i][k] * LU[k][j];
            }
        }
    }

    // Solve for each box
    for (int ibox = 0; ibox < nboxs_; ++ibox) {
        std::array<Scalar, NPARA> b;
        for (int i = 0; i < NPARA; ++i) {
            b[i] = B[ibox][perm[i]];
        }

        // Forward substitution
        for (int i = 0; i < NPARA; ++i) {
            for (int j = 0; j < i; ++j) {
                b[i] -= LU[i][j] * b[j];
            }
        }

        // Back substitution
        for (int i = NPARA - 1; i >= 0; --i) {
            for (int j = i + 1; j < NPARA; ++j) {
                b[i] -= LU[i][j] * b[j];
            }
            b[i] /= LU[i][i];
        }

        // Store coefficients
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                coeff_[ibox][i][j] = b[4*i + j];
            }
        }
    }
}

inline std::tuple<Scalar, Scalar, Scalar> Table2D::eval(Scalar x, Scalar y) const {
    if (coeff_.empty()) {
        return {0.0, 0.0, 0.0};
    }

    int nhbox = static_cast<int>(x);
    nhbox = std::max(0, std::min(nhbox, nx_ - 1));
    int ncbox = static_cast<int>(y);
    ncbox = std::max(0, std::min(ncbox, ny_ - 1));

    int ibox = ny_ * nhbox + ncbox;
    Scalar x1 = x - nhbox;
    Scalar x2 = y - ncbox;

    Scalar val = 0.0, dvdx = 0.0, dvdy = 0.0;

    for (int i = 3; i >= 0; --i) {
        Scalar s = 0.0, sdy = 0.0;
        for (int j = 3; j >= 0; --j) {
            Scalar c = coeff_[ibox][i][j];
            s = s * x2 + c;
            if (j > 0) sdy = sdy * x2 + j * c;
        }
        val = val * x1 + s;
        if (i > 0) dvdx = dvdx * x1 + i * s;
        dvdy = dvdy * x1 + sdy;
    }

    return {val, dvdx, dvdy};
}

// =========================================================================
// Table3D Implementation
// =========================================================================

inline void Table3D::init(int nx, int ny, int nz,
                          const std::vector<std::vector<std::vector<Scalar>>>& values,
                          const std::vector<std::vector<std::vector<Scalar>>>& dvdx,
                          const std::vector<std::vector<std::vector<Scalar>>>& dvdy,
                          const std::vector<std::vector<std::vector<Scalar>>>& dvdz) {
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
    nboxs_ = nx * ny * nz;

    coeff_.resize(nboxs_);

    // Corner indices for normalized coordinates
    constexpr int ix1[NCORN] = {0, 1, 1, 0, 0, 1, 1, 0};
    constexpr int ix2[NCORN] = {0, 0, 1, 1, 0, 0, 1, 1};
    constexpr int ix3[NCORN] = {0, 0, 0, 0, 1, 1, 1, 1};

    // Build matrix A
    std::array<std::array<Scalar, NPARA>, NPARA> A{};

    for (int icorn = 0; icorn < NCORN; ++icorn) {
        int n1 = ix1[icorn];
        int n2 = ix2[icorn];
        int n3 = ix3[icorn];

        for (int p1 = 0; p1 < 4; ++p1) {
            for (int p2 = 0; p2 < 4; ++p2) {
                for (int p3 = 0; p3 < 4; ++p3) {
                    int p1m = (p1 > 0) ? p1 - 1 : 0;
                    int p2m = (p2 > 0) ? p2 - 1 : 0;
                    int p3m = (p3 > 0) ? p3 - 1 : 0;
                    int icol = 16 * p1 + 4 * p2 + p3;

                    Scalar n1p1 = std::pow(n1, p1);
                    Scalar n1p1m = std::pow(n1, p1m);
                    Scalar n2p2 = std::pow(n2, p2);
                    Scalar n2p2m = std::pow(n2, p2m);
                    Scalar n3p3 = std::pow(n3, p3);
                    Scalar n3p3m = std::pow(n3, p3m);

                    A[icorn][icol] = n1p1 * n2p2 * n3p3;
                    A[icorn + NCORN][icol] = p1 * n1p1m * n2p2 * n3p3;
                    A[icorn + 2*NCORN][icol] = n1p1 * p2 * n2p2m * n3p3;
                    A[icorn + 3*NCORN][icol] = n1p1 * n2p2 * p3 * n3p3m;
                    A[icorn + 4*NCORN][icol] = p1 * n1p1m * p2 * n2p2m * n3p3;
                    A[icorn + 5*NCORN][icol] = p1 * n1p1m * n2p2 * p3 * n3p3m;
                    A[icorn + 6*NCORN][icol] = n1p1 * p2 * n2p2m * p3 * n3p3m;
                    A[icorn + 7*NCORN][icol] = p1 * n1p1m * p2 * n2p2m * p3 * n3p3m;
                }
            }
        }
    }

    // RHS vectors
    std::vector<std::array<Scalar, NPARA>> B(nboxs_);

    for (int nibox = 0; nibox < nx; ++nibox) {
        for (int njbox = 0; njbox < ny; ++njbox) {
            for (int ncbox = 0; ncbox < nz; ++ncbox) {
                int ibox = nx * (ny * ncbox + njbox) + nibox;

                for (int icorn = 0; icorn < NCORN; ++icorn) {
                    int n1 = ix1[icorn] + nibox;
                    int n2 = ix2[icorn] + njbox;
                    int n3 = ix3[icorn] + ncbox;

                    B[ibox][icorn] = values[n1][n2][n3];
                    B[ibox][icorn + NCORN] = dvdx.empty() ? 0.0 : dvdx[n1][n2][n3];
                    B[ibox][icorn + 2*NCORN] = dvdy.empty() ? 0.0 : dvdy[n1][n2][n3];
                    B[ibox][icorn + 3*NCORN] = dvdz.empty() ? 0.0 : dvdz[n1][n2][n3];
                    B[ibox][icorn + 4*NCORN] = 0.0;
                    B[ibox][icorn + 5*NCORN] = 0.0;
                    B[ibox][icorn + 6*NCORN] = 0.0;
                    B[ibox][icorn + 7*NCORN] = 0.0;
                }
            }
        }
    }

    // LU decomposition
    std::array<std::array<Scalar, NPARA>, NPARA> LU = A;
    std::array<int, NPARA> perm;
    for (int i = 0; i < NPARA; ++i) perm[i] = i;

    for (int k = 0; k < NPARA; ++k) {
        int maxrow = k;
        Scalar maxval = std::abs(LU[k][k]);
        for (int i = k + 1; i < NPARA; ++i) {
            if (std::abs(LU[i][k]) > maxval) {
                maxval = std::abs(LU[i][k]);
                maxrow = i;
            }
        }
        if (maxrow != k) {
            std::swap(LU[k], LU[maxrow]);
            std::swap(perm[k], perm[maxrow]);
        }

        if (std::abs(LU[k][k]) < 1e-15) {
            throw std::runtime_error("Table3D: Singular matrix in init");
        }

        for (int i = k + 1; i < NPARA; ++i) {
            LU[i][k] /= LU[k][k];
            for (int j = k + 1; j < NPARA; ++j) {
                LU[i][j] -= LU[i][k] * LU[k][j];
            }
        }
    }

    // Solve for each box
    for (int ibox = 0; ibox < nboxs_; ++ibox) {
        std::array<Scalar, NPARA> b;
        for (int i = 0; i < NPARA; ++i) {
            b[i] = B[ibox][perm[i]];
        }

        for (int i = 0; i < NPARA; ++i) {
            for (int j = 0; j < i; ++j) {
                b[i] -= LU[i][j] * b[j];
            }
        }

        for (int i = NPARA - 1; i >= 0; --i) {
            for (int j = i + 1; j < NPARA; ++j) {
                b[i] -= LU[i][j] * b[j];
            }
            b[i] /= LU[i][i];
        }

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                for (int k = 0; k < 4; ++k) {
                    coeff_[ibox][i][j][k] = b[16*i + 4*j + k];
                }
            }
        }
    }
}

inline std::tuple<Scalar, Scalar, Scalar, Scalar> Table3D::eval(
    Scalar x, Scalar y, Scalar z) const
{
    if (coeff_.empty()) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    int nibox = static_cast<int>(x);
    nibox = std::max(0, std::min(nibox, nx_ - 1));
    int njbox = static_cast<int>(y);
    njbox = std::max(0, std::min(njbox, ny_ - 1));
    int ncbox = static_cast<int>(z);
    ncbox = std::max(0, std::min(ncbox, nz_ - 1));

    int ibox = nx_ * (ny_ * ncbox + njbox) + nibox;
    Scalar x1 = x - nibox;
    Scalar x2 = y - njbox;
    Scalar x3 = z - ncbox;

    Scalar val = 0.0, dvdx = 0.0, dvdy = 0.0, dvdz = 0.0;

    for (int i = 3; i >= 0; --i) {
        Scalar s = 0.0, sdy = 0.0, sdz = 0.0;
        for (int j = 3; j >= 0; --j) {
            Scalar t = 0.0, tdz = 0.0;
            for (int k = 3; k >= 0; --k) {
                Scalar c = coeff_[ibox][i][j][k];
                t = t * x3 + c;
                if (k > 0) tdz = tdz * x3 + k * c;
            }
            s = s * x2 + t;
            if (j > 0) sdy = sdy * x2 + j * t;
            sdz = sdz * x2 + tdz;
        }
        val = val * x1 + s;
        if (i > 0) dvdx = dvdx * x1 + i * s;
        dvdy = dvdy * x1 + sdy;
        dvdz = dvdz * x1 + sdz;
    }

    return {val, dvdx, dvdy, dvdz};
}

// =========================================================================
// REBO2 Implementation
// =========================================================================

inline REBO2::REBO2() {
    // Initialize H-H angular function data
    hh_g_spline = {{
        {{270.467795364007301, 1549.701314596994564, 3781.927258631323866,
          4582.337619544424228, 2721.538161662818368, 630.658598136730774}},
        {{16.956325544514659, -21.059084522755980, -102.394184748124742,
          -210.527926707779059, -229.759473570467513, -94.968528666251945}},
        {{19.065031149937783, 2.017732531534021, -2.566444502991983,
          3.291353893907436, -2.653536801884563, 0.837650930130006}}
    }};

    hh_g_intervals = {{3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                       2, 2, 2, 2, 1, 1, 1}};
}

inline void REBO2::load_default_parameters() {
    init_cutoffs();
    init_angular_splines();
    init_distance_weights();
    init_tables();
    initialized_ = true;
}

inline int REBO2::element_type(int Z) const {
    if (Z == 6) return REBO2_C;       // Carbon
    if (Z == 1) return REBO2_H;       // Hydrogen
    return -1;
}

inline int REBO2::pair_type(int eli, int elj) {
    // Returns pair type index following Fortran convention
    if (eli == REBO2_C) {
        return elj;  // C-C=1, C-H=3
    } else if (elj == REBO2_C) {
        return eli;  // H-C=3
    } else {
        return eli + elj;  // H-H=6
    }
}

inline void REBO2::init_cutoffs() {
    cc_cutoff_.init(cc_r1, cc_r2);
    ch_cutoff_.init(ch_r1, ch_r2);
    hh_cutoff_.init(hh_r1, hh_r2);
    max_cutoff_ = std::max({cc_r2, ch_r2, hh_r2});
}

inline void REBO2::init_distance_weights() {
    // Initialize distance weight prefactors
    // conear[ijpot][ikpot] = exp(lambda * (r_eq_ij - r_eq_ik))
    for (auto& row : conear_) row.fill(0.0);

    conear_[REBO2_C_C][REBO2_C_C] = 1.0;
    conear_[REBO2_C_C][REBO2_C_H] = std::exp(lambda * (ch_re - cc_re));
    conear_[REBO2_C_C][REBO2_H_H] = std::exp(lambda * (hh_re - cc_re));
    conear_[REBO2_C_H][REBO2_C_C] = 1.0 / conear_[REBO2_C_C][REBO2_C_H];
    conear_[REBO2_C_H][REBO2_C_H] = 1.0;
    conear_[REBO2_C_H][REBO2_H_H] = std::exp(lambda * (hh_re - ch_re));
    conear_[REBO2_H_H][REBO2_C_C] = 1.0 / conear_[REBO2_C_C][REBO2_H_H];
    conear_[REBO2_H_H][REBO2_C_H] = 1.0 / conear_[REBO2_C_H][REBO2_H_H];
    conear_[REBO2_H_H][REBO2_H_H] = 1.0;
}

inline void REBO2::make_cc_g_spline(GSplineCoeffs& coeffs,
                                     const std::array<Scalar, 6>& g_vals,
                                     const std::array<Scalar, 6>& dg_vals,
                                     const std::array<Scalar, 6>& d2g_vals) {
    // Build piecewise polynomial spline through g values, derivatives
    // For each interval, solve for 6 polynomial coefficients using:
    // - Function values at both ends
    // - First derivatives at both ends
    // - Second derivatives at both ends (for intervals 0 and 1)
    // - Function values at interior points (for interval 2)

    // Interval 0: [-1, -0.5]
    {
        Scalar x0 = cc_g_theta[0], x1 = cc_g_theta[1];
        Scalar f0 = g_vals[0], f1 = g_vals[1];
        Scalar df0 = dg_vals[0], df1 = dg_vals[1];
        Scalar d2f0 = d2g_vals[0], d2f1 = d2g_vals[1];

        // Solve 6x6 system for quintic polynomial
        // Simplified: use cubic Hermite for now
        Scalar h = x1 - x0;
        coeffs.c[0][0] = f0;
        coeffs.c[0][1] = df0;
        coeffs.c[0][2] = (3*(f1-f0)/h - 2*df0 - df1) / h;
        coeffs.c[0][3] = (df0 + df1 - 2*(f1-f0)/h) / (h*h);
        coeffs.c[0][4] = 0.0;
        coeffs.c[0][5] = 0.0;
    }

    // Interval 1: [-0.5, -1/3]
    {
        Scalar x0 = cc_g_theta[1], x1 = cc_g_theta[2];
        Scalar f0 = g_vals[1], f1 = g_vals[2];
        Scalar df0 = dg_vals[1], df1 = dg_vals[2];

        Scalar h = x1 - x0;
        coeffs.c[1][0] = f0;
        coeffs.c[1][1] = df0;
        coeffs.c[1][2] = (3*(f1-f0)/h - 2*df0 - df1) / h;
        coeffs.c[1][3] = (df0 + df1 - 2*(f1-f0)/h) / (h*h);
        coeffs.c[1][4] = 0.0;
        coeffs.c[1][5] = 0.0;
    }

    // Interval 2: [-1/3, 1] - fit through 4 points
    {
        // Use Lagrange interpolation through the 4 points
        Scalar x2 = cc_g_theta[2], x3 = cc_g_theta[3];
        Scalar x4 = cc_g_theta[4], x5 = cc_g_theta[5];
        Scalar f2 = g_vals[2], f3 = g_vals[3], f4 = g_vals[4], f5 = g_vals[5];

        // Fit cubic through these 4 points (in local coordinates)
        // Using polynomial fit
        Scalar h = x5 - x2;

        // Simple cubic fit
        coeffs.c[2][0] = f2;
        coeffs.c[2][1] = (-11*f2 + 18*f3 - 9*f4 + 2*f5) / (6*(h/3));
        coeffs.c[2][2] = (2*f2 - 5*f3 + 4*f4 - f5) / (2*(h/3)*(h/3));
        coeffs.c[2][3] = (-f2 + 3*f3 - 3*f4 + f5) / (6*(h/3)*(h/3)*(h/3));
        coeffs.c[2][4] = 0.0;
        coeffs.c[2][5] = 0.0;
    }
}

inline void REBO2::init_angular_splines() {
    // Build g1 spline (for N > 3.7, fully coordinated)
    make_cc_g_spline(cc_g1_coeff_, cc_g_g1, cc_g_dg1, cc_g_d2g1);

    // Build g2 spline (for N < 3.2, undercoordinated)
    // g2 uses same theta values but different g values
    std::array<Scalar, 6> g2_dg = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // Approximate
    std::array<Scalar, 6> g2_d2g = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    make_cc_g_spline(cc_g2_coeff_, cc_g_g2, g2_dg, g2_d2g);
}

inline void REBO2::init_tables() {
    // Initialize the F, P, T lookup tables with default Brenner 2002 values

    // Fcc table (5x5x10)
    {
        std::vector<std::vector<std::vector<Scalar>>> F(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));
        std::vector<std::vector<std::vector<Scalar>>> dFdi(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));
        std::vector<std::vector<std::vector<Scalar>>> dFdj(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));
        std::vector<std::vector<std::vector<Scalar>>> dFdk(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));

        // Values from Table 4 of Brenner 2002
        F[1][1][0] = 0.105000;
        F[1][1][1] = -0.0041775;
        for (int k = 2; k <= 8; ++k) F[1][1][k] = -0.0160856;

        F[2][2][0] = 0.09444957;
        F[2][2][1] = 0.02200000;
        F[2][2][2] = 0.03970587;
        F[2][2][3] = 0.03308822;
        F[2][2][4] = 0.02647058;
        F[2][2][5] = 0.01985293;
        F[2][2][6] = 0.01323529;
        F[2][2][7] = 0.00661764;
        F[2][2][8] = 0.0;

        F[0][1][0] = 0.04338699;
        for (int k = 1; k <= 8; ++k) F[0][1][k] = 0.0099172158;

        F[0][2][0] = 0.0493976637;
        F[0][2][1] = -0.011942669;
        for (int k = 2; k <= 8; ++k) F[0][2][k] = 0.0099172158;

        for (int k = 0; k <= 1; ++k) F[0][3][k] = -0.119798935;
        for (int k = 2; k <= 8; ++k) F[0][3][k] = 0.0099172158;

        F[1][2][0] = 0.0096495698;
        F[1][2][1] = 0.030;
        F[1][2][2] = -0.0200;
        F[1][2][3] = -0.0233778774;
        F[1][2][4] = -0.0267557548;
        for (int k = 5; k <= 8; ++k) F[1][2][k] = -0.030133632;

        for (int k = 1; k <= 8; ++k) F[1][3][k] = -0.124836752;
        for (int k = 0; k <= 8; ++k) F[2][3][k] = -0.044709383;

        // Interpolate F[2][2] for k=3-7
        for (int k = 3; k <= 7; ++k) {
            F[2][2][k] = F[2][2][2] + (k-2) * (F[2][2][8] - F[2][2][2]) / 6;
        }
        for (int k = 3; k <= 4; ++k) {
            F[1][2][k] = F[1][2][2] + (k-2) * (F[1][2][5] - F[1][2][2]) / 3;
        }

        // Derivatives
        dFdi[2][1][0] = -0.052500;
        for (int k = 4; k <= 8; ++k) dFdi[2][1][k] = -0.054376;
        dFdi[2][3][0] = 0.0;
        for (int k = 1; k <= 8; ++k) dFdi[2][3][k] = 0.062418;
        for (int k = 3; k <= 7; ++k) dFdk[2][2][k] = -0.006618;
        dFdk[1][1][1] = -0.060543;
        dFdk[1][2][3] = -0.020044;
        dFdk[1][2][4] = -0.020044;

        // Symmetrize
        for (int k = 0; k < 10; ++k) {
            for (int i = 0; i <= 3; ++i) {
                for (int j = i+1; j <= 3; ++j) {
                    Scalar x = F[i][j][k] + F[j][i][k];
                    F[i][j][k] = x;
                    F[j][i][k] = x;

                    x = dFdi[i][j][k] + dFdj[j][i][k];
                    dFdi[i][j][k] = x;
                    dFdj[j][i][k] = x;

                    x = dFdi[j][i][k] + dFdj[i][j][k];
                    dFdi[j][i][k] = x;
                    dFdj[i][j][k] = x;

                    x = dFdk[i][j][k] + dFdk[j][i][k];
                    dFdk[i][j][k] = x;
                    dFdk[j][i][k] = x;
                }
            }
        }

        Fcc_.init(4, 4, 9, F, dFdi, dFdj, dFdk);
    }

    // Fch table
    {
        std::vector<std::vector<std::vector<Scalar>>> F(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));

        for (int k = 4; k <= 8; ++k) F[0][2][k] = -0.0090477875161288110;
        for (int k = 0; k <= 8; ++k) F[1][3][k] = -0.213;
        for (int k = 0; k <= 8; ++k) F[1][2][k] = -0.25;
        for (int k = 0; k <= 8; ++k) F[1][1][k] = -0.5;

        // Symmetrize
        for (int k = 0; k < 10; ++k) {
            for (int i = 0; i <= 2; ++i) {
                for (int j = i+1; j <= 3; ++j) {
                    Scalar x = F[i][j][k] + F[j][i][k];
                    F[i][j][k] = x;
                    F[j][i][k] = x;
                }
            }
        }

        Fch_.init(4, 4, 9, F);
    }

    // Fhh table
    {
        std::vector<std::vector<std::vector<Scalar>>> F(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));
        F[1][1][0] = 0.249831916;
        Fhh_.init(4, 4, 9, F);
    }

    // Pcc table (6x6)
    {
        std::vector<std::vector<Scalar>> P(6, std::vector<Scalar>(6, 0.0));
        P[1][1] = 0.003026697473481;
        P[2][0] = 0.007860700254745;
        P[3][0] = 0.016125364564267;
        P[1][2] = 0.003179530830731;
        P[2][1] = 0.006326248241119;
        Pcc_.init(5, 5, P);
    }

    // Pch table (6x6)
    {
        std::vector<std::vector<Scalar>> P(6, std::vector<Scalar>(6, 0.0));
        P[1][0] = 0.2093367328250380;
        P[2][0] = -0.064449615432525;
        P[3][0] = -0.303927546346162;
        P[0][1] = 0.01;
        P[0][2] = -0.1220421462782555;
        P[1][1] = -0.1251234006287090;
        P[2][1] = -0.298905245783;
        P[0][3] = -0.307584705066;
        P[1][2] = -0.3005291724067579;
        Pch_.init(5, 5, P);
    }

    // Tcc table
    {
        std::vector<std::vector<std::vector<Scalar>>> T(5,
            std::vector<std::vector<Scalar>>(5, std::vector<Scalar>(10, 0.0)));
        T[2][2][0] = -0.070280085;
        for (int k = 1; k <= 8; ++k) T[2][2][k] = -0.00809675;
        Tcc_.init(4, 4, 9, T);
    }
}

inline std::pair<Scalar, Scalar> REBO2::repulsive(int ptype, Scalar r) const {
    Scalar A, Q, alpha;

    if (ptype == REBO2_C_C) {
        A = cc_A; Q = cc_Q; alpha = cc_alpha;
    } else if (ptype == REBO2_C_H) {
        A = ch_A; Q = ch_Q; alpha = ch_alpha;
    } else {  // H-H
        A = hh_A; Q = hh_Q; alpha = hh_alpha;
    }

    Scalar exp_val = A * std::exp(-alpha * r);
    Scalar hlp = 1.0 + Q / r;

    Scalar V = hlp * exp_val;
    Scalar dV = (-Q / (r * r) - hlp * alpha) * exp_val;

    return {V, dV};
}

inline std::pair<Scalar, Scalar> REBO2::attractive(int ptype, Scalar r) const {
    if (ptype == REBO2_C_C) {
        // Three exponential terms for C-C
        Scalar exp1 = cc_B1 * std::exp(-cc_beta1 * r);
        Scalar exp2 = cc_B2 * std::exp(-cc_beta2 * r);
        Scalar exp3 = cc_B3 * std::exp(-cc_beta3 * r);

        Scalar V = -(exp1 + exp2 + exp3);
        Scalar dV = cc_beta1 * exp1 + cc_beta2 * exp2 + cc_beta3 * exp3;

        return {V, dV};
    } else if (ptype == REBO2_C_H) {
        Scalar exp1 = ch_B1 * std::exp(-ch_beta1 * r);
        Scalar V = -exp1;
        Scalar dV = ch_beta1 * exp1;
        return {V, dV};
    } else {  // H-H
        Scalar exp1 = hh_B1 * std::exp(-hh_beta1 * r);
        Scalar V = -exp1;
        Scalar dV = hh_beta1 * exp1;
        return {V, dV};
    }
}

inline std::tuple<Scalar, Scalar, Scalar> REBO2::angular_function(
    int el_type, Scalar cos_theta, Scalar N) const
{
    if (el_type == REBO2_C) {
        // Carbon angular function with coordination-dependent interpolation
        Scalar g, dg_dcos, dg_dN;

        if (N < 3.2) {
            // Use g2 spline
            auto [val, dval] = eval_cc_g_spline_with_deriv(cc_g2_coeff_, cos_theta);
            g = val;
            dg_dcos = dval;
            dg_dN = 0.0;
        } else if (N > 3.7) {
            // Use g1 spline
            auto [val, dval] = eval_cc_g_spline_with_deriv(cc_g1_coeff_, cos_theta);
            g = val;
            dg_dcos = dval;
            dg_dN = 0.0;
        } else {
            // Interpolate between g1 and g2
            auto [v1, dv1] = eval_cc_g_spline_with_deriv(cc_g1_coeff_, cos_theta);
            auto [v2, dv2] = eval_cc_g_spline_with_deriv(cc_g2_coeff_, cos_theta);

            Scalar arg = 2.0 * M_PI * (N - 3.2);
            Scalar s = (1.0 + std::cos(arg)) / 2.0;
            Scalar ds = -M_PI * std::sin(arg);

            g = v1 * (1.0 - s) + v2 * s;
            dg_dcos = dv1 * (1.0 - s) + dv2 * s;
            dg_dN = (v2 - v1) * ds;
        }

        return {g, dg_dcos, dg_dN};
    } else {
        // Hydrogen angular function (polynomial)
        int ig = hh_g_intervals[std::min(24, static_cast<int>(-cos_theta * 12.0) + 12)];
        ig = std::max(1, std::min(3, ig)) - 1;  // Convert to 0-indexed

        const auto& coeffs = hh_g_spline[ig];

        Scalar g = coeffs[0] + coeffs[1] * cos_theta;
        Scalar dg = coeffs[1];

        Scalar cos_pow = cos_theta;
        for (int i = 2; i < 6; ++i) {
            cos_pow *= cos_theta;
            g += coeffs[i] * cos_pow;
            dg += i * coeffs[i] * std::pow(cos_theta, i - 1);
        }

        return {g, dg, 0.0};  // No N dependence for H
    }
}

inline std::pair<Scalar, Scalar> REBO2::distance_weight(
    int ptype_ij, int ptype_ik, Scalar dr) const
{
    // Only apply distance weighting when bonds involve hydrogen
    if ((ptype_ij + ptype_ik) <= 4) {
        return {1.0, 0.0};
    }

    Scalar h = conear_[ptype_ij][ptype_ik] * std::exp(lambda * dr);
    Scalar dh = lambda * h;

    return {h, dh};
}

inline std::pair<Scalar, Scalar> REBO2::bond_order_func(int el_type, Scalar z) const {
    Scalar exp = (el_type == REBO2_C) ? conpe_C_ : conpe_H_;
    Scalar arg = 1.0 + z;
    Scalar b = std::pow(arg, exp);
    Scalar db = exp * std::pow(arg, exp - 1.0);
    return {b, db};
}

inline Scalar REBO2::eval_cc_g_spline(const GSplineCoeffs& coeffs, Scalar cos_theta) const {
    // Determine which interval
    int j;
    if (cos_theta < cc_g_theta[1]) {
        j = 0;
    } else if (cos_theta < cc_g_theta[2]) {
        j = 1;
    } else {
        j = 2;
    }

    // Evaluate polynomial (in local coordinate system)
    Scalar x = cos_theta - cc_g_theta[j == 0 ? 0 : (j == 1 ? 1 : 2)];
    Scalar g = coeffs.c[j][0] + coeffs.c[j][1] * x;
    Scalar x_pow = x;
    for (int i = 2; i < 6; ++i) {
        x_pow *= x;
        g += coeffs.c[j][i] * x_pow;
    }

    return g;
}

inline std::pair<Scalar, Scalar> REBO2::eval_cc_g_spline_with_deriv(
    const GSplineCoeffs& coeffs, Scalar cos_theta) const
{
    // Determine which interval
    int j;
    Scalar x0;
    if (cos_theta < cc_g_theta[1]) {
        j = 0;
        x0 = cc_g_theta[0];
    } else if (cos_theta < cc_g_theta[2]) {
        j = 1;
        x0 = cc_g_theta[1];
    } else {
        j = 2;
        x0 = cc_g_theta[2];
    }

    // Evaluate polynomial and derivative (in local coordinates)
    Scalar x = cos_theta - x0;
    Scalar g = coeffs.c[j][0] + coeffs.c[j][1] * x;
    Scalar dg = coeffs.c[j][1];

    Scalar x_pow = x;
    for (int i = 2; i < 6; ++i) {
        x_pow *= x;
        g += coeffs.c[j][i] * x_pow;
        dg += i * coeffs.c[j][i] * std::pow(x, i - 1);
    }

    return {g, dg};
}

inline PotentialResults REBO2::compute(
    AtomicSystem& system, NeighborList& neighbors,
    bool compute_forces, bool compute_virial)
{
    if (!initialized_) {
        load_default_parameters();
    }

    PotentialResults results;
    results.energy = 0.0;
    results.virial.setZero();

    const std::size_t n_atoms = system.num_atoms();
    if (compute_forces) {
        system.zero_forces();
    }

    // Map atomic numbers to internal element types
    std::vector<int> el_type(n_atoms);
    for (std::size_t i = 0; i < n_atoms; ++i) {
        el_type[i] = element_type(system.atomic_numbers()(i));
        if (el_type[i] < 0) {
            throw std::runtime_error("REBO2: Unsupported element (only C and H supported)");
        }
    }

    // Store bond data for angular calculations
    struct BondInfo {
        std::size_t j;
        int ptype;
        Scalar r;
        Vec3 dr;
        Vec3 unit;
        Scalar fc, dfc;
    };
    std::vector<std::vector<BondInfo>> atom_bonds(n_atoms);

    // First pass: build neighbor data and compute coordination numbers
    std::vector<Scalar> N_C(n_atoms, 0.0), N_H(n_atoms, 0.0);

    for (std::size_t i = 0; i < n_atoms; ++i) {
        auto [begin, end] = neighbors.neighbors(i);

        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;

            Vec3 rij = system.minimum_image(
                system.positions().col(j) - system.positions().col(i));
            Scalar r = rij.norm();

            int ptype = pair_type(el_type[i], el_type[j]);

            // Get cutoff function value
            CutoffResult fc;
            if (ptype == REBO2_C_C) {
                fc = cc_cutoff_(r);
            } else if (ptype == REBO2_C_H) {
                fc = ch_cutoff_(r);
            } else {
                fc = hh_cutoff_(r);
            }

            if (fc.fc > 0.0) {
                // Store bond info
                atom_bonds[i].push_back({j, ptype, r, rij, rij/r, fc.fc, fc.dfc});

                // Add to coordination numbers (count each bond once per atom)
                if (el_type[j] == REBO2_C) {
                    N_C[i] += fc.fc;
                } else {
                    N_H[i] += fc.fc;
                }
            }
        }
    }

    // Second pass: compute energy and forces with full bond order
    for (std::size_t i = 0; i < n_atoms; ++i) {
        Scalar N_i = N_C[i] + N_H[i];

        for (const auto& bond_ij : atom_bonds[i]) {
            std::size_t j = bond_ij.j;
            if (j <= i) continue;  // Each pair once

            Scalar r_ij = bond_ij.r;
            const Vec3& rij = bond_ij.dr;
            const Vec3& nij = bond_ij.unit;
            int ptype_ij = bond_ij.ptype;
            Scalar fc_ij = bond_ij.fc;
            Scalar dfc_ij = bond_ij.dfc;

            Scalar N_j = N_C[j] + N_H[j];

            // Pair potentials
            auto [VR, dVR] = repulsive(ptype_ij, r_ij);
            auto [VA, dVA] = attractive(ptype_ij, r_ij);

            // Compute bond order z_ij (angular contribution from i's neighbors)
            Scalar z_ij = 0.0;
            for (const auto& bond_ik : atom_bonds[i]) {
                std::size_t k = bond_ik.j;
                if (k == j) continue;

                Scalar r_ik = bond_ik.r;
                const Vec3& nik = bond_ik.unit;
                int ptype_ik = bond_ik.ptype;

                // Cosine of angle
                Scalar cos_theta = nij.dot(nik);

                // Angular function
                auto [g, dg_dcos, dg_dN] = angular_function(el_type[i], cos_theta, N_i);

                // Distance weight
                auto [h, dh] = distance_weight(ptype_ij, ptype_ik, r_ij - r_ik);

                z_ij += bond_ik.fc * g * h;
            }

            // Compute bond order z_ji (angular contribution from j's neighbors)
            Scalar z_ji = 0.0;
            for (const auto& bond_jl : atom_bonds[j]) {
                std::size_t l = bond_jl.j;
                if (l == i) continue;

                Scalar r_jl = bond_jl.r;
                Vec3 njl = bond_jl.unit;
                int ptype_jl = bond_jl.ptype;

                // Cosine of angle (note: use -nij for j->i direction)
                Scalar cos_theta = (-nij).dot(njl);

                // Angular function
                auto [g, dg_dcos, dg_dN] = angular_function(el_type[j], cos_theta, N_j);

                // Distance weight
                auto [h, dh] = distance_weight(ptype_ij, ptype_jl, r_ij - r_jl);

                z_ji += bond_jl.fc * g * h;
            }

            // Bond orders
            auto [b_ij, db_ij] = bond_order_func(el_type[i], z_ij);
            auto [b_ji, db_ji] = bond_order_func(el_type[j], z_ji);

            // Average bond order
            Scalar b_avg = 0.5 * (b_ij + b_ji);

            // Add P and F corrections if tables are valid
            Scalar P_ij = 0.0, P_ji = 0.0;
            if (Pcc_.is_valid() && Pch_.is_valid()) {
                if (ptype_ij == REBO2_C_C) {
                    auto [p, dp_dnc, dp_dnh] = Pcc_.eval(N_C[i] - fc_ij, N_H[i]);
                    P_ij = p;
                    auto [pj, dpj_dnc, dpj_dnh] = Pcc_.eval(N_C[j] - fc_ij, N_H[j]);
                    P_ji = pj;
                } else if (ptype_ij == REBO2_C_H) {
                    if (el_type[i] == REBO2_C) {
                        auto [p, dp_dnc, dp_dnh] = Pch_.eval(N_C[i], N_H[i] - fc_ij);
                        P_ij = p;
                    }
                    if (el_type[j] == REBO2_C) {
                        auto [p, dp_dnc, dp_dnh] = Pch_.eval(N_C[j], N_H[j] - fc_ij);
                        P_ji = p;
                    }
                }
            }

            // Final bond order with corrections
            b_avg += 0.5 * (P_ij + P_ji);

            // Energy contribution
            Scalar E_pair = fc_ij * (VR + b_avg * VA);
            results.energy += E_pair;

            if (compute_forces) {
                // Simplified force: only pair contribution
                // Full angular forces require more complex bookkeeping
                Scalar dE_dr = dfc_ij * (VR + b_avg * VA) +
                               fc_ij * (dVR + b_avg * dVA);

                Vec3 fij = dE_dr * nij;
                system.forces().col(i) += fij.array();
                system.forces().col(j) -= fij.array();

                if (compute_virial) {
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            results.virial(a, b) += rij(a) * fij(b);
                        }
                    }
                }
            }
        }
    }

    return results;
}

// =========================================================================
// REBO2Scr: Screened 2nd-generation REBO potential
// =========================================================================

/**
 * @brief Screened REBO2 potential (REBO2+S)
 *
 * Adds Pastewka-style bond screening to REBO2.
 * C-C inner cutoff changes to (1.95, 2.25) Å.
 * C-C outer (screening) cutoff: (2.179347, 2.819732) Å.
 * Screening parameters: Cmin=1.0, Cmax=2.0.
 *
 * Reference: Pastewka, Pou, Perez, Gumbsch, Moseler, PRB 78, 161402(R) (2008)
 */
class REBO2Scr : public REBO2 {
public:
    REBO2Scr() = default;

    void load_default_parameters() {
        // Screened C-C inner cutoff (wider than unscreened)
        cc_r1 = 1.95;
        cc_r2 = 2.25;
        REBO2::load_default_parameters();
        cc_outer_cutoff_.init(cc_outer_r1_, cc_outer_r2_);
        scr_initialized_ = true;
    }

    Scalar cutoff() const { return std::max({cc_outer_r2_, ch_r2, hh_r2}); }

    PotentialResults compute(AtomicSystem& system, NeighborList& neighbors,
                             bool compute_forces = true,
                             bool compute_virial = true)
    {
        if (!scr_initialized_) load_default_parameters();

        PotentialResults results;
        const std::size_t n_atoms = system.num_atoms();
        if (compute_forces) system.zero_forces();

        std::vector<int> el_type(n_atoms);
        for (std::size_t i = 0; i < n_atoms; ++i) {
            el_type[i] = element_type(system.atomic_numbers()(i));
            if (el_type[i] < 0)
                throw std::runtime_error("REBO2Scr: unsupported element (only C and H)");
        }

        struct BondInfo {
            std::size_t j;
            int ptype;
            Scalar r;
            Vec3 dr, unit;
            Scalar fc, dfc;
        };

        std::vector<std::vector<BondInfo>> atom_bonds(n_atoms);

        // First pass: build inner-cutoff bonds and coordination numbers
        std::vector<Scalar> N_C(n_atoms, 0.0), N_H(n_atoms, 0.0);

        for (std::size_t i = 0; i < n_atoms; ++i) {
            auto [begin, end] = neighbors.neighbors(i);
            for (auto it = begin; it != end; ++it) {
                std::size_t j = it->index;
                Vec3 rij = system.minimum_image(
                    system.positions().col(j) - system.positions().col(i));
                Scalar r = rij.norm();
                int ptype = pair_type(el_type[i], el_type[j]);

                CutoffResult fc;
                if (ptype == REBO2_C_C)      fc = cc_cutoff_(r);
                else if (ptype == REBO2_C_H) fc = ch_cutoff_(r);
                else                         fc = hh_cutoff_(r);

                if (fc.fc > 0.0) {
                    atom_bonds[i].push_back({j, ptype, r, rij, rij/r, fc.fc, fc.dfc});
                    if (el_type[j] == REBO2_C) N_C[i] += fc.fc;
                    else                        N_H[i] += fc.fc;
                }
            }
        }

        // Screening helper: compute S_ij and its derivatives for a bond
        // Also modifies fc_ij and dfc_ij in-place
        const Scalar C_dr_cut = Cmax_ * Cmax_ / (4.0 * (Cmax_ - 1.0));
        const Scalar dC = Cmax_ - Cmin_;

        struct ScreenK {
            std::size_t k;
            Scalar dS_drik; // d(S)/d(rik^2) * 2 / rij^2 * S, then * rik → force scale
            Scalar dS_drjk;
            Vec3 unit_ik, unit_jk;
        };

        auto compute_screening = [&](std::size_t i, const BondInfo& bond_ij,
                                     std::vector<ScreenK>& screen_ks) -> Scalar
        {
            Scalar S_log = 0.0;
            Scalar rij_sq = bond_ij.r * bond_ij.r;
            bool fully = false;

            // Lambda to check one potential screener k
            auto check_k = [&](std::size_t k, Scalar r_ik, const Vec3& dr_ik, const Vec3& unit_ik) {
                if (k == bond_ij.j) return;
                Scalar rik_sq = r_ik * r_ik;
                if (rik_sq >= C_dr_cut * rij_sq) return;

                Vec3 dr_jk = dr_ik - bond_ij.dr;
                Scalar rjk_sq = dr_jk.squaredNorm();

                // Geometric check: k must lie between i and j
                if (bond_ij.dr.dot(dr_ik) <= dot_threshold_) return;
                if (bond_ij.dr.dot(dr_jk) >= -dot_threshold_) return;

                Scalar xik = rik_sq / rij_sq;
                Scalar xjk = rjk_sq / rij_sq;
                Scalar xdiff = xik - xjk;
                Scalar denom = 1.0 - xdiff * xdiff;
                if (std::abs(denom) < 1e-15) return;

                Scalar fac = 1.0 / denom;
                Scalar C = (2.0 * (xik + xjk) - xdiff * xdiff - 1.0) * fac;

                if (C <= Cmin_) { fully = true; return; }
                if (C < Cmax_) {
                    Scalar Cmax_C = Cmax_ - C;
                    Scalar C_Cmin = C - Cmin_;
                    Scalar ratio = Cmax_C / C_Cmin;
                    S_log -= ratio * ratio;

                    Scalar dCdxik = 4.0 * xik * fac * (1.0 + (C - 1.0) * xdiff);
                    Scalar dCdxjk = 4.0 * xjk * fac * (1.0 - (C - 1.0) * xdiff);
                    Scalar dSdC = 2.0 * Cmax_C / (dC * C_Cmin * C_Cmin);
                    Scalar rjk = std::sqrt(rjk_sq);
                    Vec3 unit_jk = rjk > 1e-10 ? Vec3(dr_jk / rjk) : Vec3::Zero();

                    screen_ks.push_back({k,
                        dSdC * dCdxik * 2.0 / rij_sq,
                        dSdC * dCdxjk * 2.0 / rij_sq,
                        unit_ik, unit_jk});
                }
            };

            // Inner-cutoff neighbors screen this bond
            for (const auto& b : atom_bonds[i]) check_k(b.j, b.r, b.dr, b.unit);
            // Outer CC bonds (beyond inner cutoff but within screening range)
            auto [nb_begin, nb_end] = neighbors.neighbors(i);
            for (auto it = nb_begin; it != nb_end; ++it) {
                std::size_t k = it->index;
                if (pair_type(el_type[i], el_type[k]) != REBO2_C_C) continue;
                Vec3 dr_ik = system.minimum_image(
                    system.positions().col(k) - system.positions().col(i));
                Scalar r_ik = dr_ik.norm();
                if (r_ik < cc_cutoff_.cutoff() || r_ik >= cc_outer_r2_) continue;
                check_k(k, r_ik, dr_ik, dr_ik / r_ik);
                if (fully) break;
            }

            if (fully || S_log < screening_threshold_) return 0.0;

            Scalar S = std::exp(S_log);
            for (auto& sk : screen_ks) {
                sk.dS_drik *= S;
                sk.dS_drjk *= S;
            }
            return S;
        };

        // Second pass: energy and forces with screening
        for (std::size_t i = 0; i < n_atoms; ++i) {
            Scalar N_i = N_C[i] + N_H[i];

            for (const auto& bond_ij : atom_bonds[i]) {
                std::size_t j = bond_ij.j;
                if (j <= i) continue;

                Scalar r_ij = bond_ij.r;
                int ptype_ij = bond_ij.ptype;
                Scalar fc_ij = bond_ij.fc;
                Scalar dfc_ij = bond_ij.dfc;
                Scalar N_j = N_C[j] + N_H[j];

                auto [VR, dVR] = repulsive(ptype_ij, r_ij);
                auto [VA, dVA] = attractive(ptype_ij, r_ij);

                // Bond order (same as REBO2)
                Scalar z_ij = 0.0;
                for (const auto& bond_ik : atom_bonds[i]) {
                    if (bond_ik.j == j) continue;
                    auto [g, dg_dcos, dg_dN] = angular_function(
                        el_type[i], bond_ij.unit.dot(bond_ik.unit), N_i);
                    auto [h, dh] = distance_weight(ptype_ij, bond_ik.ptype,
                                                    r_ij - bond_ik.r);
                    z_ij += bond_ik.fc * g * h;
                }

                Scalar z_ji = 0.0;
                for (const auto& bond_jl : atom_bonds[j]) {
                    if (bond_jl.j == i) continue;
                    auto [g, dg_dcos, dg_dN] = angular_function(
                        el_type[j], (-bond_ij.unit).dot(bond_jl.unit), N_j);
                    auto [h, dh] = distance_weight(ptype_ij, bond_jl.ptype,
                                                    r_ij - bond_jl.r);
                    z_ji += bond_jl.fc * g * h;
                }

                auto [b_ij, db_ij] = bond_order_func(el_type[i], z_ij);
                auto [b_ji, db_ji] = bond_order_func(el_type[j], z_ji);
                Scalar b_avg = 0.5 * (b_ij + b_ji);

                if (Pcc_.is_valid() && Pch_.is_valid()) {
                    Scalar P_ij = 0.0, P_ji = 0.0;
                    if (ptype_ij == REBO2_C_C) {
                        auto [p,  dp1, dp2] = Pcc_.eval(N_C[i] - fc_ij, N_H[i]);
                        auto [pj, dp3, dp4] = Pcc_.eval(N_C[j] - fc_ij, N_H[j]);
                        P_ij = p; P_ji = pj;
                    } else if (ptype_ij == REBO2_C_H) {
                        if (el_type[i] == REBO2_C) {
                            auto [p, dp1, dp2] = Pch_.eval(N_C[i], N_H[i] - fc_ij);
                            P_ij = p;
                        }
                        if (el_type[j] == REBO2_C) {
                            auto [p, dp1, dp2] = Pch_.eval(N_C[j], N_H[j] - fc_ij);
                            P_ji = p;
                        }
                    }
                    b_avg += 0.5 * (P_ij + P_ji);
                }

                // Compute screening
                std::vector<ScreenK> screen_ks;
                Scalar S_ij = compute_screening(i, bond_ij, screen_ks);

                if (S_ij < 1e-15) continue;

                Scalar E0_pair = fc_ij * (VR + b_avg * VA);
                Scalar E_pair  = S_ij * E0_pair;
                results.energy += E_pair;

                if (compute_forces || compute_virial) {
                    Scalar dE_dr = S_ij * (dfc_ij * (VR + b_avg * VA) +
                                           fc_ij * (dVR + b_avg * dVA));

                    Vec3 fij = dE_dr * bond_ij.unit;
                    if (compute_forces) {
                        system.forces().col(i) += fij.array();
                        system.forces().col(j) -= fij.array();
                    }
                    if (compute_virial)
                        results.virial += bond_ij.dr * fij.transpose();

                    // Screening force contributions
                    for (const auto& sk : screen_ks) {
                        Scalar dS_drik = sk.dS_drik;
                        Scalar dS_drjk = sk.dS_drjk;

                        // Force from r_ik dependence
                        Vec3 f_ik = E0_pair * dS_drik * sk.unit_ik;
                        Vec3 f_jk = E0_pair * dS_drjk * sk.unit_jk;

                        if (compute_forces) {
                            system.forces().col(i)    -= f_ik.array();
                            system.forces().col(sk.k) += f_ik.array();
                            system.forces().col(j)    -= f_jk.array();
                            system.forces().col(sk.k) += f_jk.array();
                        }

                        if (compute_virial) {
                            results.virial += bond_ij.dr * f_jk.transpose();
                        }
                    }
                }
            }
        }

        return results;
    }

private:
    bool scr_initialized_ = false;

    // Outer C-C cutoff for screening
    Scalar cc_outer_r1_ = 2.179347;
    Scalar cc_outer_r2_ = 2.819732;
    TrigOffCutoff cc_outer_cutoff_;

    // Screening parameters
    Scalar Cmin_ = 1.0;
    Scalar Cmax_ = 2.0;
    Scalar dot_threshold_ = 1e-10;
    Scalar screening_threshold_ = std::log(1e-6);
};

} // namespace atomistica
