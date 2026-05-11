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
#include <map>
#include <string>
#include <vector>

#include "bop_base.hpp"

namespace atomistica {

/**
 * @brief Tersoff pair parameters
 */
struct TersoffPairParams : public BOPPairParams {
    // Note: Tersoff uses A, lambda (repulsive), B, mu (attractive) from base
    // Additional mixing parameters
    Scalar chi = 1.0;     // Mixing parameter for heteroatomic bonds

    // Screening parameters (only used when Screening=true)
    ScreeningParams screening;
};

/**
 * @brief Tersoff element parameters
 */
struct TersoffElementParams : public BOPElementParams {
    // Angular parameters (stored per element, mixed for pairs)
    Scalar c = 0.0;
    Scalar d = 1.0;
    Scalar h = 0.0;

    // Precomputed
    Scalar c2 = 0.0;
    Scalar d2 = 0.0;
    Scalar c2_d2 = 0.0;

    void precompute_angular() {
        c2 = c * c;
        d2 = d * d;
        c2_d2 = c2 / d2;
        BOPElementParams::precompute();
    }
};

/**
 * @brief Tersoff potential implementation
 *
 * Standard Tersoff potential as described in:
 * J. Tersoff, Phys. Rev. B 39, 5566 (1989)
 *
 * Energy: E = 0.5 * sum_{ij} fc(r_ij) * [V_R(r_ij) + b_ij * V_A(r_ij)]
 *
 * V_R(r) = A * exp(-lambda * r)
 * V_A(r) = -B * exp(-mu * r)
 * b_ij = (1 + beta^n * zeta_ij^n)^(-1/(2n))
 * zeta_ij = sum_{k != j} fc(r_ik) * g(cos_theta_jik) * exp(mu^3 * (r_ij - r_ik)^3)
 * g(cos_theta) = 1 + c^2/d^2 - c^2/(d^2 + (h - cos_theta)^2)
 *
 * @tparam Screening Enable screening (default: false)
 */
template<bool Screening = false>
class Tersoff : public BOPBase<Tersoff<Screening>, Screening> {
public:
    using Base = BOPBase<Tersoff<Screening>, Screening>;
    friend Base;

    Tersoff() = default;

    /**
     * @brief Set element mapping (atomic number -> internal index)
     */
    void add_element(int Z, const TersoffElementParams& params) {
        int idx = static_cast<int>(element_params_.size());
        element_map_[Z] = idx;
        element_params_.push_back(params);
        element_params_.back().precompute_angular();
        update_pair_count();
    }

    /**
     * @brief Set pair parameters
     */
    void set_pair_params(int Z1, int Z2, const TersoffPairParams& params) {
        int el1 = element_index(Z1);
        int el2 = element_index(Z2);
        if (el1 < 0 || el2 < 0) return;

        int ptype = pair_type(el1, el2);
        ensure_pair_storage(ptype);

        pair_params_[ptype] = params;
        pair_params_[ptype].init_cutoff();

        update_max_cutoff();
    }

    /**
     * @brief Load parameters from built-in database
     *
     * @param name Parameter set name (e.g., "Tersoff_PRB_39_5566_Si_C")
     */
    void load_parameters(const std::string& name);

    // Required interface for CRTP base
    int element_index(int Z) const {
        auto it = element_map_.find(Z);
        return (it != element_map_.end()) ? it->second : -1;
    }

    int pair_type(int eli, int elj) const {
        return pair_index(eli, elj, num_elements());
    }

    int num_elements() const {
        return static_cast<int>(element_params_.size());
    }

    Scalar cutoff_impl() const {
        return max_cutoff_;
    }

    Scalar pair_cutoff(int ptype) const {
        return pair_params_[ptype].r2;
    }

    // For screened version: determines neighbor search radius
    Scalar screened_cutoff(int ptype) const {
        return pair_params_[ptype].screening.cut_out_h;
    }

    CutoffResult cutoff_function(int ptype, Scalar r) const {
        return pair_params_[ptype].cutoff(r);
    }

    /**
     * @brief Get screening parameters for a pair type (only used when Screening=true)
     */
    const ScreeningParams& screening_params(int ptype) const {
        return pair_params_[ptype].screening;
    }

    /**
     * @brief Repulsive potential V_R(r) = A * exp(-lambda * r)
     */
    std::pair<Scalar, Scalar> repulsive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.lambda * r);
        Scalar VR = p.A * exp_val;
        Scalar dVR = -p.lambda * VR;
        return {VR, dVR};
    }

    /**
     * @brief Attractive potential V_A(r) = -B * exp(-mu * r)
     */
    std::pair<Scalar, Scalar> attractive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.mu * r);
        Scalar VA = -p.B * exp_val;
        Scalar dVA = -p.mu * VA;  // = p.B * p.mu * exp_val
        return {VA, dVA};
    }

    /**
     * @brief Angular function g(cos_theta)
     *
     * g(cos) = 1 + c^2/d^2 - c^2/(d^2 + (h - cos)^2)
     * dg/dcos = -2 * c^2 * (h - cos) / (d^2 + (h - cos)^2)^2
     */
    std::pair<Scalar, Scalar> angular_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) const
    {
        // Use parameters from central atom i
        const auto& p = element_params_[eli];

        Scalar h_cos = p.h - cos_theta;
        Scalar h_cos2 = h_cos * h_cos;
        Scalar denom = p.d2 + h_cos2;
        Scalar denom_inv = 1.0 / denom;

        Scalar g = 1.0 + p.c2_d2 - p.c2 * denom_inv;
        // Negative sign is correct: dg/dcos = -2*c^2*(h-cos)/(d^2+(h-cos)^2)^2
        Scalar dg = -2.0 * p.c2 * h_cos * denom_inv * denom_inv;

        return {g, dg};
    }

    /**
     * @brief Distance-dependent function h(r_ij, r_ik)
     *
     * h = exp(mu^3 * (r_ij - r_ik)^3) for Tersoff
     *
     * Note: Some Tersoff variants use h = 1 (no distance modulation)
     */
    std::tuple<Scalar, Scalar, Scalar> distance_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik,
        Scalar r_ij, Scalar r_ik) const
    {
        // Get mu for i-k pair
        const auto& p_ik = pair_params_[ptype_ik];
        Scalar mu3 = p_ik.mu * p_ik.mu * p_ik.mu;

        Scalar dr = r_ij - r_ik;
        Scalar dr3 = dr * dr * dr;

        Scalar h = std::exp(mu3 * dr3);
        Scalar dh_dr = 3.0 * mu3 * dr * dr * h;

        // dh/dr_ik = -dh_dr, dh/dr_ij = +dh_dr
        return {h, -dh_dr, dh_dr};
    }

    /**
     * @brief Bond order function b(z)
     *
     * b(z) = (1 + beta^n * z^n)^(-1/(2n))
     */
    std::pair<Scalar, Scalar> bond_order(int eli, int ptype, Scalar z) const {
        const auto& p = element_params_[eli];

        if (z < 1e-10) {
            // Avoid numerical issues at z=0
            return {1.0, 0.0};
        }

        Scalar beta_n = std::pow(p.beta, p.n);
        Scalar z_n = std::pow(z, p.n);
        Scalar term = 1.0 + beta_n * z_n;
        Scalar exp_val = p.minus_half_over_n;

        Scalar b = std::pow(term, exp_val);
        Scalar db = exp_val * beta_n * p.n * std::pow(z, p.n - 1.0) *
                    std::pow(term, exp_val - 1.0);

        // Apply xi scaling if different from 1
        const auto& p_pair = pair_params_[ptype];
        b *= p_pair.chi;
        db *= p_pair.chi;

        return {b, db};
    }

private:
    void update_pair_count() {
        int n = num_elements();
        int np = atomistica::num_pairs(n);
        if (static_cast<int>(pair_params_.size()) < np) {
            pair_params_.resize(np);
        }
    }

    void ensure_pair_storage(int ptype) {
        if (ptype >= static_cast<int>(pair_params_.size())) {
            pair_params_.resize(ptype + 1);
        }
    }

    void update_max_cutoff() {
        max_cutoff_ = 0.0;
        for (const auto& p : pair_params_) {
            Scalar cut;
            if constexpr (Screening) {
                cut = p.screening.cut_out_h;
            } else {
                cut = p.r2;
            }
            if (cut > max_cutoff_) {
                max_cutoff_ = cut;
            }
        }
    }

    std::map<int, int> element_map_;
    std::vector<TersoffElementParams> element_params_;
    std::vector<TersoffPairParams> pair_params_;
    Scalar max_cutoff_ = 0.0;
};

// ============================================================================
// Built-in parameter sets
// ============================================================================

/**
 * @brief Tersoff Si-C parameters from PRB 39, 5566 (1989)
 *
 * Templated to support both screened and non-screened versions.
 * For screened version, adds screening parameters from Pastewka et al.
 */
template<bool Scr>
inline void load_tersoff_prb_39_5566_si_c(Tersoff<Scr>& pot) {
    // Silicon parameters
    TersoffElementParams si;
    si.beta = 1.1e-6;
    si.n = 0.78734;
    si.c = 100390.0;
    si.d = 16.217;
    si.h = -0.59825;
    pot.add_element(14, si);  // Z=14 for Si

    // Carbon parameters
    TersoffElementParams c;
    c.beta = 1.5724e-7;
    c.n = 0.72751;
    c.c = 38049.0;
    c.d = 4.3484;
    c.h = -0.57058;
    pot.add_element(6, c);  // Z=6 for C

    // Si-Si pair
    TersoffPairParams si_si;
    si_si.A = 1830.8;
    si_si.B = 471.18;
    si_si.lambda = 2.4799;
    si_si.mu = 1.7322;
    si_si.chi = 1.0;

    if constexpr (Scr) {
        // Screened parameters from Fortran tersoff_params.f90
        si_si.r1 = 2.50;
        si_si.r2 = 2.50 * 1.2;  // 3.0
        si_si.screening.cut_in_l = si_si.r1;
        si_si.screening.cut_in_h = si_si.r2;
        si_si.screening.cut_out_l = 3.0;
        si_si.screening.cut_out_h = 3.0 * 2.0;  // 6.0
        si_si.screening.cut_bo_l = 3.0;
        si_si.screening.cut_bo_h = 3.0 * 2.0;
        si_si.screening.Cmin = 1.0;
        si_si.screening.Cmax = 3.0;
        si_si.screening.precompute();
    } else {
        si_si.r1 = 2.7;
        si_si.r2 = 3.0;
    }
    pot.set_pair_params(14, 14, si_si);

    // C-C pair
    TersoffPairParams c_c;
    c_c.A = 1393.6;
    c_c.B = 346.74;
    c_c.lambda = 3.4879;
    c_c.mu = 2.2119;
    c_c.chi = 1.0;

    if constexpr (Scr) {
        c_c.r1 = 2.0;
        c_c.r2 = 2.0 * 1.2;  // 2.4
        c_c.screening.cut_in_l = c_c.r1;
        c_c.screening.cut_in_h = c_c.r2;
        c_c.screening.cut_out_l = 2.0;
        c_c.screening.cut_out_h = 2.0 * 2.0;  // 4.0
        c_c.screening.cut_bo_l = 2.0;
        c_c.screening.cut_bo_h = 2.0 * 2.0;
        c_c.screening.Cmin = 1.0;
        c_c.screening.Cmax = 3.0;
        c_c.screening.precompute();
    } else {
        c_c.r1 = 1.8;
        c_c.r2 = 2.1;
    }
    pot.set_pair_params(6, 6, c_c);

    // Si-C pair (mixed)
    TersoffPairParams si_c;
    si_c.A = std::sqrt(si_si.A * c_c.A);
    si_c.B = std::sqrt(si_si.B * c_c.B);
    si_c.lambda = 0.5 * (si_si.lambda + c_c.lambda);
    si_c.mu = 0.5 * (si_si.mu + c_c.mu);
    si_c.chi = 0.9776;  // Mixing parameter

    if constexpr (Scr) {
        si_c.r1 = std::sqrt(2.50 * 2.0);
        si_c.r2 = si_c.r1 * 1.2;
        si_c.screening.cut_in_l = si_c.r1;
        si_c.screening.cut_in_h = si_c.r2;
        si_c.screening.cut_out_l = std::sqrt(3.0 * 2.0);
        si_c.screening.cut_out_h = si_c.screening.cut_out_l * 2.0;
        si_c.screening.cut_bo_l = si_c.screening.cut_out_l;
        si_c.screening.cut_bo_h = si_c.screening.cut_out_h;
        si_c.screening.Cmin = 1.0;
        si_c.screening.Cmax = 3.0;
        si_c.screening.precompute();
    } else {
        si_c.r1 = std::sqrt(si_si.r1 * c_c.r1);
        si_c.r2 = std::sqrt(si_si.r2 * c_c.r2);
    }
    pot.set_pair_params(14, 6, si_c);
}

// ============================================================================
// Goumri-Said Al-N parameters (Chem. Phys. 302, 135, 2004)
// ============================================================================

/**
 * @brief Goumri-Said Al-N Tersoff parameters
 *
 * S. Goumri-Said et al., Chem. Phys. 302, 135 (2004)
 * Elements: Al (Z=13, index 0), N (Z=7, index 1)
 * Pair ordering: Al-Al=0, Al-N=1, N-N=2
 */
template<bool Scr>
inline void load_goumri_said_chemphys_302_135_al_n(Tersoff<Scr>& pot) {
    // Al element params
    TersoffElementParams al;
    al.beta = 1.094932;
    al.n    = 6.085605;
    al.c    = 0.074836;
    al.d    = 19.569127;
    al.h    = -0.659266;
    pot.add_element(13, al);  // Al -> index 0

    // N element params
    TersoffElementParams n;
    n.beta = 5.2938e-3;
    n.n    = 1.33041;
    n.c    = 2.0312e4;
    n.d    = 20.312;
    n.h    = -0.56239;
    pot.add_element(7, n);  // N -> index 1

    // Al-Al pair (ptype 0)
    TersoffPairParams al_al;
    al_al.A = 746.698; al_al.B = 40.451;
    al_al.lambda = 2.4647; al_al.mu = 0.9683;
    al_al.chi = 1.0;
    al_al.r1 = 3.20; al_al.r2 = 3.60;
    if constexpr (Scr) {
        al_al.screening.cut_in_l  = 3.20; al_al.screening.cut_in_h  = 3.60;
        al_al.screening.cut_out_l = 3.20; al_al.screening.cut_out_h = 7.20;
        al_al.screening.cut_bo_l  = 3.20; al_al.screening.cut_bo_h  = 7.20;
        al_al.screening.Cmin = 1.0; al_al.screening.Cmax = 3.0;
        al_al.screening.precompute();
    }
    pot.set_pair_params(13, 13, al_al);

    // Al-N pair (ptype 1)
    TersoffPairParams al_n;
    al_n.A = 3000.214; al_n.B = 298.81;
    al_n.lambda = 3.53051; al_n.mu = 1.99995;
    al_n.chi = 1.0;
    al_n.r1 = 2.185; al_n.r2 = 2.485;
    if constexpr (Scr) {
        al_n.screening.cut_in_l  = 2.185; al_n.screening.cut_in_h  = 2.485;
        al_n.screening.cut_out_l = 2.185; al_n.screening.cut_out_h = 4.970;
        al_n.screening.cut_bo_l  = 2.185; al_n.screening.cut_bo_h  = 4.970;
        al_n.screening.Cmin = 1.0; al_n.screening.Cmax = 3.0;
        al_n.screening.precompute();
    }
    pot.set_pair_params(13, 7, al_n);

    // N-N pair (ptype 2)
    TersoffPairParams n_n;
    n_n.A = 636.814; n_n.B = 511.76;
    n_n.lambda = 5.43673; n_n.mu = 2.7;
    n_n.chi = 1.0;
    n_n.r1 = 1.60; n_n.r2 = 2.00;
    if constexpr (Scr) {
        n_n.screening.cut_in_l  = 1.60; n_n.screening.cut_in_h  = 2.00;
        n_n.screening.cut_out_l = 1.60; n_n.screening.cut_out_h = 4.00;
        n_n.screening.cut_bo_l  = 1.60; n_n.screening.cut_bo_h  = 4.00;
        n_n.screening.Cmin = 1.0; n_n.screening.Cmax = 3.0;
        n_n.screening.precompute();
    }
    pot.set_pair_params(7, 7, n_n);
}

// ============================================================================
// Matsunaga B-C-N parameters (Jpn. J. Appl. Phys. 39, 48, 2000)
// ============================================================================

/**
 * @brief Matsunaga B-C-N Tersoff parameters
 *
 * K. Matsunaga, C. Fisher, H. Matsubara, Jpn. J. Appl. Phys. 39, 48 (2000)
 * Elements: C (Z=6, index 0), N (Z=7, index 1), B (Z=5, index 2)
 * Pair ordering: C-C=0, C-N=1, C-B=2, N-N=3, N-B=4, B-B=5
 *
 * Mixed pair values computed using geometric mean (A,B,r1,r2) and
 * arithmetic mean (lambda, mu) from homospecies values.
 */
template<bool Scr>
inline void load_matsunaga_fisher_matsubara_b_c_n(Tersoff<Scr>& pot) {
    // C element params
    TersoffElementParams c;
    c.beta = 1.5724e-7; c.n = 7.2751e-1;
    c.c = 3.8049e4; c.d = 4.3484; c.h = -5.7058e-1;
    pot.add_element(6, c);  // C -> index 0

    // N element params
    TersoffElementParams nb;
    nb.beta = 1.0562e-1; nb.n = 12.4498;
    nb.c = 7.9934e4; nb.d = 1.3432e2; nb.h = -0.9973;
    pot.add_element(7, nb);  // N -> index 1

    // B element params
    TersoffElementParams b;
    b.beta = 1.6e-6; b.n = 3.9929;
    b.c = 5.2629e-1; b.d = 1.5870e-3; b.h = 0.5;
    pot.add_element(5, b);  // B -> index 2

    // Helper: set Tersoff pair params
    auto set_pair = [&](int Z1, int Z2, Scalar A, Scalar B_,
                        Scalar lam, Scalar mu_, Scalar chi,
                        Scalar r1, Scalar r2,
                        Scalar or1 = 0.0, Scalar or2 = 0.0) {
        TersoffPairParams p;
        p.A = A; p.B = B_; p.lambda = lam; p.mu = mu_; p.chi = chi;
        p.r1 = r1; p.r2 = r2;
        if constexpr (Scr) {
            p.screening.cut_in_l  = r1;  p.screening.cut_in_h  = r2;
            p.screening.cut_out_l = or1; p.screening.cut_out_h = or2;
            p.screening.cut_bo_l  = or1; p.screening.cut_bo_h  = or2;
            p.screening.Cmin = 1.0; p.screening.Cmax = 3.0;
            p.screening.precompute();
        }
        pot.set_pair_params(Z1, Z2, p);
    };

    // Homospecies values (for mixing reference)
    constexpr Scalar A_CC = 1.3936e3, A_NN = 1.1e4, A_BB = 2.7702e2;
    constexpr Scalar B_CC = 3.4674e2, B_NN = 2.1945e2, B_BB = 1.8349e2;
    constexpr Scalar l_CC = 3.4879, l_NN = 5.7708, l_BB = 1.9922;
    constexpr Scalar m_CC = 2.2119, m_NN = 2.5115, m_BB = 1.5856;
    constexpr Scalar r1_CC = 1.80, r1_NN = 2.0, r1_BB = 1.8;
    constexpr Scalar r2_CC = 2.10, r2_NN = 2.3, r2_BB = 2.1;

    // xi (chi) values — provided directly, not mixed
    constexpr Scalar xi_CC = 1.0, xi_CN = 0.9685, xi_CB = 1.0025;
    constexpr Scalar xi_NN = 1.0, xi_NB = 1.1593, xi_BB = 1.0;

    if constexpr (!Scr) {
        // C-C
        set_pair(6, 6, A_CC, B_CC, l_CC, m_CC, xi_CC, r1_CC, r2_CC);
        // C-N (geometric A,B,r1,r2; arithmetic lambda,mu)
        set_pair(6, 7, std::sqrt(A_CC*A_NN), std::sqrt(B_CC*B_NN),
                 0.5*(l_CC+l_NN), 0.5*(m_CC+m_NN), xi_CN,
                 std::sqrt(r1_CC*r1_NN), std::sqrt(r2_CC*r2_NN));
        // C-B
        set_pair(6, 5, std::sqrt(A_CC*A_BB), std::sqrt(B_CC*B_BB),
                 0.5*(l_CC+l_BB), 0.5*(m_CC+m_BB), xi_CB,
                 std::sqrt(r1_CC*r1_BB), std::sqrt(r2_CC*r2_BB));
        // N-N
        set_pair(7, 7, A_NN, B_NN, l_NN, m_NN, xi_NN, r1_NN, r2_NN);
        // N-B
        set_pair(7, 5, std::sqrt(A_NN*A_BB), std::sqrt(B_NN*B_BB),
                 0.5*(l_NN+l_BB), 0.5*(m_NN+m_BB), xi_NB,
                 std::sqrt(r1_NN*r1_BB), std::sqrt(r2_NN*r2_BB));
        // B-B
        set_pair(5, 5, A_BB, B_BB, l_BB, m_BB, xi_BB, r1_BB, r2_BB);
    } else {
        // Screened version: different inner cutoffs, outer cutoffs from or1/or2
        // r1,r2 (inner): CC=2.0/2.4, NN=2.0/2.4, BB=1.8/2.16; cross: geometric
        // or1,or2 (outer): CC=2.0/4.0, NN=3.0/6.0, BB=1.8/3.6; cross: geometric
        constexpr Scalar sr1_CC=2.00, sr2_CC=2.40, sor1_CC=2.00, sor2_CC=4.00;
        constexpr Scalar sr1_NN=2.00, sr2_NN=2.40, sor1_NN=3.00, sor2_NN=6.00;
        constexpr Scalar sr1_BB=1.80, sr2_BB=2.16, sor1_BB=1.80, sor2_BB=3.60;

        // C-C
        set_pair(6, 6, A_CC, B_CC, l_CC, m_CC, xi_CC,
                 sr1_CC, sr2_CC, sor1_CC, sor2_CC);
        // C-N
        set_pair(6, 7, std::sqrt(A_CC*A_NN), std::sqrt(B_CC*B_NN),
                 0.5*(l_CC+l_NN), 0.5*(m_CC+m_NN), xi_CN,
                 std::sqrt(sr1_CC*sr1_NN), std::sqrt(sr2_CC*sr2_NN),
                 std::sqrt(sor1_CC*sor1_NN), std::sqrt(sor2_CC*sor2_NN));
        // C-B
        set_pair(6, 5, std::sqrt(A_CC*A_BB), std::sqrt(B_CC*B_BB),
                 0.5*(l_CC+l_BB), 0.5*(m_CC+m_BB), xi_CB,
                 std::sqrt(sr1_CC*sr1_BB), std::sqrt(sr2_CC*sr2_BB),
                 std::sqrt(sor1_CC*sor1_BB), std::sqrt(sor2_CC*sor2_BB));
        // N-N
        set_pair(7, 7, A_NN, B_NN, l_NN, m_NN, xi_NN,
                 sr1_NN, sr2_NN, sor1_NN, sor2_NN);
        // N-B
        set_pair(7, 5, std::sqrt(A_NN*A_BB), std::sqrt(B_NN*B_BB),
                 0.5*(l_NN+l_BB), 0.5*(m_NN+m_BB), xi_NB,
                 std::sqrt(sr1_NN*sr1_BB), std::sqrt(sr2_NN*sr2_BB),
                 std::sqrt(sor1_NN*sor1_BB), std::sqrt(sor2_NN*sor2_BB));
        // B-B
        set_pair(5, 5, A_BB, B_BB, l_BB, m_BB, xi_BB,
                 sr1_BB, sr2_BB, sor1_BB, sor2_BB);
    }
}

template<bool Screening>
void Tersoff<Screening>::load_parameters(const std::string& name) {
    if (name == "Tersoff_PRB_39_5566_Si_C") {
        load_tersoff_prb_39_5566_si_c(*this);
    } else if (name == "Goumri_Said_ChemPhys_302_135_Al_N") {
        load_goumri_said_chemphys_302_135_al_n(*this);
    } else if (name == "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N") {
        load_matsunaga_fisher_matsubara_b_c_n(*this);
    } else if (name == "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr") {
        if constexpr (Screening) {
            load_matsunaga_fisher_matsubara_b_c_n(*this);
        } else {
            throw std::runtime_error(
                "Matsunaga_Fisher_Matsubara__Scr requires TersoffScr, not Tersoff");
        }
    } else {
        throw std::runtime_error("Unknown parameter set: " + name);
    }
}

// Type aliases
using TersoffPotential = Tersoff<false>;
using TersoffScreened = Tersoff<true>;

} // namespace atomistica
