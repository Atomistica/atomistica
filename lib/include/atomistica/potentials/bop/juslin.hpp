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
#include "brenner.hpp"  // reuse BrennerPairParams and BrennerElementParams

namespace atomistica {

/**
 * @brief Juslin potential implementation
 *
 * The Juslin potential uses the same functional form as Brenner for pair
 * potentials and angular terms, but the distance-dependent bond-order
 * function uses per-TRIPLET parameters (alpha, omega, m) indexed by
 * (central atom, j-neighbor, k-neighbor), rather than per-pair parameters.
 *
 * Distance function: h = omega[i][j][k] * exp((alpha[i][j][k] * dr)^m[i][j][k])
 * where dr = r_ij - r_ik.
 *
 * References:
 * - Juslin et al., J. Appl. Phys. 98, 123520 (2005) [W-C-H]
 * - Kuopanportti et al., Comp. Mat. Sci. 111, 525 (2016) [Fe-C-H]
 *
 * @tparam Screening Enable screening (default: false)
 */
template<bool Screening = false>
class Juslin : public BOPBase<Juslin<Screening>, Screening> {
public:
    using Base = BOPBase<Juslin<Screening>, Screening>;
    friend Base;

    Juslin() = default;

    void add_element(int Z, const BrennerElementParams& params = BrennerElementParams{}) {
        int idx = static_cast<int>(element_params_.size());
        element_map_[Z] = idx;
        element_params_.push_back(params);
        nel_ = static_cast<int>(element_params_.size());
        update_pair_count();
        triplet_alpha_.resize(nel_ * nel_ * nel_, 0.0);
        triplet_omega_.resize(nel_ * nel_ * nel_, 1.0);
        triplet_m_.resize(nel_ * nel_ * nel_, 1);
    }

    void set_pair_params(int Z1, int Z2, const BrennerPairParams& params) {
        int el1 = element_index(Z1);
        int el2 = element_index(Z2);
        if (el1 < 0 || el2 < 0) return;

        int ptype = pair_type(el1, el2);
        ensure_pair_storage(ptype);

        pair_params_[ptype] = params;
        if (params.r1 < params.r2) {
            pair_params_[ptype].precompute();
        }

        update_max_cutoff();
    }

    /**
     * @brief Set per-triplet distance function parameters
     * @param eli Central atom element index (0-based)
     * @param elj j-neighbor element index (0-based)
     * @param elk k-neighbor element index (0-based)
     */
    void set_triplet_params(int eli, int elj, int elk,
                            Scalar alpha, Scalar omega, int m) {
        int idx = triple_index(eli, elj, elk);
        triplet_alpha_[idx] = alpha;
        triplet_omega_[idx] = omega;
        triplet_m_[idx] = m;
    }

    void load_parameters(const std::string& name);

    // Required interface for CRTP base
    int element_index(int Z) const {
        auto it = element_map_.find(Z);
        return (it != element_map_.end()) ? it->second : -1;
    }

    int pair_type(int eli, int elj) const {
        return pair_index(eli, elj, num_elements());
    }

    int num_elements() const { return nel_; }

    Scalar cutoff_impl() const { return max_cutoff_; }

    Scalar pair_cutoff(int ptype) const { return pair_params_[ptype].r2; }

    Scalar screened_cutoff(int ptype) const {
        return pair_params_[ptype].screening.cut_out_h;
    }

    CutoffResult cutoff_function(int ptype, Scalar r) const {
        return pair_params_[ptype].cutoff(r);
    }

    const ScreeningParams& screening_params(int ptype) const {
        return pair_params_[ptype].screening;
    }

    // Same pair functions as Brenner
    std::pair<Scalar, Scalar> repulsive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.expR * (r - p.r0));
        Scalar VR = p.VR_f * exp_val;
        return {VR, -p.expR * VR};
    }

    std::pair<Scalar, Scalar> attractive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.expA * (r - p.r0));
        Scalar VA = -p.VA_f * exp_val;
        return {VA, -p.expA * VA};
    }

    // Angular function: same as Brenner (gamma*(1 + c²/d² - c²/(d²+(h+cos)²)))
    std::pair<Scalar, Scalar> angular_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) const
    {
        const auto& p = pair_params_[ptype_ik];
        Scalar h_cos = p.h + cos_theta;
        Scalar h_cos2 = h_cos * h_cos;
        Scalar denom = p.d_sq + h_cos2;
        Scalar denom_inv = 1.0 / denom;

        Scalar g = p.gamma * (1.0 + p.c_d - p.c_sq * denom_inv);
        Scalar dg = 2.0 * p.gamma * p.c_sq * h_cos * denom_inv * denom_inv;
        return {g, dg};
    }

    /**
     * @brief Distance function h = omega * exp((alpha * dr)^m)
     *
     * Uses per-triplet (eli, elj, elk) parameters instead of per-pair.
     * Returns {h, dh/dr_ik, dh/dr_ij}.
     */
    std::tuple<Scalar, Scalar, Scalar> distance_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik,
        Scalar r_ij, Scalar r_ik) const
    {
        int tidx = triple_index(eli, elj, elk);
        Scalar alpha = triplet_alpha_[tidx];
        Scalar omega = triplet_omega_[tidx];
        int m = triplet_m_[tidx];

        Scalar dr = r_ij - r_ik;

        if (alpha == 0.0) {
            return {omega, 0.0, 0.0};
        }

        Scalar h, dh_dr;
        if (m == 1) {
            h = omega * std::exp(alpha * dr);
            dh_dr = alpha * h;
        } else if (m == 3) {
            Scalar arg = alpha * dr;
            h = omega * std::exp(arg * arg * arg);
            dh_dr = 3.0 * alpha * arg * arg * h;
        } else {
            Scalar arg = alpha * dr;
            Scalar argm = std::pow(std::abs(arg), m) * (arg >= 0 ? 1.0 : -1.0);
            h = omega * std::exp(argm);
            dh_dr = m * std::pow(std::abs(arg), m - 1) * alpha * h;
        }

        return {h, -dh_dr, dh_dr};
    }

    // Bond order: same as Brenner  (1 + z^n)^(-1/(2n))
    std::pair<Scalar, Scalar> bond_order(int eli, int ptype, Scalar z) const {
        const auto& p = pair_params_[ptype];

        if (p.n == 1.0) {
            Scalar arg = 1.0 + z;
            Scalar b = std::pow(arg, p.bo_exp);
            Scalar db = p.bo_exp * std::pow(arg, p.bo_exp - 1.0);
            return {b, db};
        }

        if (z < 1e-10) return {1.0, 0.0};

        Scalar z_n = std::pow(z, p.n);
        Scalar arg = 1.0 + z_n;
        Scalar b = std::pow(arg, p.bo_exp);
        Scalar db = -0.5 * std::pow(z, p.n - 1.0) * std::pow(arg, p.bo_exp - 1.0);
        return {b, db};
    }

private:
    int triple_index(int eli, int elj, int elk) const {
        return eli * nel_ * nel_ + elj * nel_ + elk;
    }

    void update_pair_count() {
        int np = atomistica::num_pairs(nel_);
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
            Scalar cut = Screening ? p.screening.cut_out_h : p.r2;
            if (cut > max_cutoff_) max_cutoff_ = cut;
        }
    }

    std::map<int, int> element_map_;
    std::vector<BrennerElementParams> element_params_;
    std::vector<BrennerPairParams> pair_params_;
    std::vector<Scalar> triplet_alpha_;
    std::vector<Scalar> triplet_omega_;
    std::vector<int> triplet_m_;
    int nel_ = 0;
    Scalar max_cutoff_ = 0.0;
};

// ============================================================================
// Built-in parameter sets
// ============================================================================

// Helper to make an inactive pair (r1=r2=0, D0=0)
inline BrennerPairParams make_inactive_pair() {
    BrennerPairParams p;
    p.D0 = 0.0;
    p.r0 = 0.0;
    p.S = 2.0;  // Must be > 1 to avoid error in precompute
    p.beta = 0.0;
    p.gamma = 0.0;
    p.c = 0.0;
    p.d = 1.0;
    p.h = 0.0;
    p.n = 1.0;
    p.mu = 0.0;
    p.m = 1;
    p.r1 = 0.0;
    p.r2 = 0.0;
    // Don't call precompute() since r1=r2=0 → no cutoff; VR_f=0 anyway
    p.VR_f = 0.0;
    p.VA_f = 0.0;
    p.expR = 0.0;
    p.expA = 0.0;
    p.bo_exp = -0.5;
    p.bo_fac = -0.25;
    p.bo_exp1 = -1.5;
    p.c_sq = 0.0;
    p.d_sq = 1.0;
    p.c_d = 0.0;
    // Leave cutoff uninitialized — pair_cutoff() returns r2=0 so bonds are always skipped
    return p;
}

/**
 * @brief Juslin W-C-H parameters (J. Appl. Phys. 98, 123520, 2005)
 *
 * Elements: W (Z=74, index 0), C (Z=6, index 1), H (Z=1, index 2)
 * Active pairs: W-W, C-C, C-H, H-H
 * Inactive pairs: W-C, W-H (D0=0, r1=r2=0)
 */
template<bool Scr>
inline void load_juslin_jap_98_123520_wch(Juslin<Scr>& pot) {
    // Add elements: W=0, C=1, H=2
    pot.add_element(74);  // W
    pot.add_element(6);   // C
    pot.add_element(1);   // H

    // nel=3, pair_index(eli, elj, 3):
    // W-W=0, W-C=1, W-H=2, C-C=3, C-H=4, H-H=5

    // W-W (pair 0)
    BrennerPairParams ww;
    ww.D0 = 5.41861; ww.r0 = 2.34095; ww.S = 1.92708; ww.beta = 1.38528;
    ww.gamma = 0.00188227; ww.c = 2.14969; ww.d = 0.17126; ww.h = -0.27780;
    ww.n = 1.0; ww.mu = 0.0; ww.m = 1;
    if constexpr (Scr) {
        ww.r1 = 3.20; ww.r2 = 3.80;
        ww.screening.cut_in_l = 3.20; ww.screening.cut_in_h = 3.80;
        ww.screening.cut_out_l = 3.20; ww.screening.cut_out_h = 3.80;
        ww.screening.cut_bo_l  = 3.20; ww.screening.cut_bo_h  = 3.80;
        ww.screening.Cmin = 1.0; ww.screening.Cmax = 3.0;
        ww.screening.precompute();
    } else {
        ww.r1 = 3.20; ww.r2 = 3.80;
    }
    pot.set_pair_params(74, 74, ww);

    // W-C (pair 1): INACTIVE
    BrennerPairParams wc = make_inactive_pair();
    if constexpr (Scr) {
        wc.screening.cut_in_l = 0.0; wc.screening.cut_in_h = 0.0;
        wc.screening.cut_out_l = 0.0; wc.screening.cut_out_h = 0.0;
        wc.screening.cut_bo_l = 0.0; wc.screening.cut_bo_h = 0.0;
        wc.screening.Cmin = 1.0; wc.screening.Cmax = 3.0;
        wc.screening.precompute();
    }
    pot.set_pair_params(74, 6, wc);

    // W-H (pair 2): INACTIVE
    BrennerPairParams wh = make_inactive_pair();
    if constexpr (Scr) {
        wh.screening = wc.screening;
    }
    pot.set_pair_params(74, 1, wh);

    // C-C (pair 3)
    BrennerPairParams cc;
    cc.D0 = 6.0; cc.r0 = 1.39; cc.S = 1.22; cc.beta = 2.1;
    cc.gamma = 0.00020813; cc.c = 330.0; cc.d = 3.5; cc.h = 1.0;
    cc.n = 1.0; cc.mu = 0.0; cc.m = 1;
    if constexpr (Scr) {
        cc.r1 = 1.70; cc.r2 = 2.00;
        cc.screening.cut_in_l = 1.70; cc.screening.cut_in_h = 2.00;
        cc.screening.cut_out_l = 1.70; cc.screening.cut_out_h = 2.00;
        cc.screening.cut_bo_l  = 1.70; cc.screening.cut_bo_h  = 2.00;
        cc.screening.Cmin = 1.0; cc.screening.Cmax = 3.0;
        cc.screening.precompute();
    } else {
        cc.r1 = 1.70; cc.r2 = 2.00;
    }
    pot.set_pair_params(6, 6, cc);

    // C-H (pair 4)
    BrennerPairParams ch;
    ch.D0 = 3.642; ch.r0 = 1.1199; ch.S = 1.69077; ch.beta = 1.9583;
    ch.gamma = 12.33; ch.c = 0.0; ch.d = 1.0; ch.h = 1.0;
    ch.n = 1.0; ch.mu = 0.0; ch.m = 1;
    if constexpr (Scr) {
        ch.r1 = 1.30; ch.r2 = 1.80;
        ch.screening.cut_in_l = 1.30; ch.screening.cut_in_h = 1.80;
        ch.screening.cut_out_l = 1.30; ch.screening.cut_out_h = 1.80;
        ch.screening.cut_bo_l  = 1.30; ch.screening.cut_bo_h  = 1.80;
        ch.screening.Cmin = 1.0; ch.screening.Cmax = 3.0;
        ch.screening.precompute();
    } else {
        ch.r1 = 1.30; ch.r2 = 1.80;
    }
    pot.set_pair_params(6, 1, ch);

    // H-H (pair 5)
    BrennerPairParams hh;
    hh.D0 = 4.7509; hh.r0 = 0.74144; hh.S = 2.3432; hh.beta = 1.9436;
    hh.gamma = 12.33; hh.c = 0.0; hh.d = 1.0; hh.h = 1.0;
    hh.n = 1.0; hh.mu = 0.0; hh.m = 1;
    if constexpr (Scr) {
        hh.r1 = 1.10; hh.r2 = 1.70;
        hh.screening.cut_in_l = 1.10; hh.screening.cut_in_h = 1.70;
        hh.screening.cut_out_l = 1.10; hh.screening.cut_out_h = 1.70;
        hh.screening.cut_bo_l  = 1.10; hh.screening.cut_bo_h  = 1.70;
        hh.screening.Cmin = 1.0; hh.screening.Cmax = 3.0;
        hh.screening.precompute();
    } else {
        hh.r1 = 1.10; hh.r2 = 1.70;
    }
    pot.set_pair_params(1, 1, hh);

    // --- Triplet parameters (alpha, omega, m) ---
    // Ordering: eli=0=W, eli=1=C, eli=2=H
    // triple_index(eli, elj, elk) = eli*9 + elj*3 + elk

    // W central atom (eli=0): alpha non-zero only for k=W (elk=0)
    pot.set_triplet_params(0, 0, 0, 0.45876, 1.0, 1);  // W central, W-W
    pot.set_triplet_params(0, 0, 1, 0.0,     1.0, 1);  // W central, W-C
    pot.set_triplet_params(0, 0, 2, 0.0,     1.0, 1);  // W central, W-H
    pot.set_triplet_params(0, 1, 0, 0.45876, 1.0, 1);  // W central, C-W
    pot.set_triplet_params(0, 1, 1, 0.0,     1.0, 1);  // W central, C-C
    pot.set_triplet_params(0, 1, 2, 0.0,     1.0, 1);  // W central, C-H
    pot.set_triplet_params(0, 2, 0, 0.45876, 1.0, 1);  // W central, H-W
    pot.set_triplet_params(0, 2, 1, 0.0,     1.0, 1);  // W central, H-C
    pot.set_triplet_params(0, 2, 2, 0.0,     1.0, 1);  // W central, H-H

    // C central atom (eli=1): alpha non-zero for C-H and H-* combos
    pot.set_triplet_params(1, 0, 0, 0.0, 1.0, 1);       // C central, W-W
    pot.set_triplet_params(1, 0, 1, 0.0, 1.0, 1);       // C central, W-C
    pot.set_triplet_params(1, 0, 2, 0.0, 1.0, 1);       // C central, W-H
    pot.set_triplet_params(1, 1, 0, 0.0, 1.0, 1);       // C central, C-W
    pot.set_triplet_params(1, 1, 1, 0.0, 1.0, 1);       // C central, C-C
    pot.set_triplet_params(1, 1, 2, 4.0, 1.0, 1);       // C central, C-H
    pot.set_triplet_params(1, 2, 0, 0.0, 1.0, 1);       // C central, H-W
    pot.set_triplet_params(1, 2, 1, 4.0, 2.94586, 1);   // C central, H-C
    pot.set_triplet_params(1, 2, 2, 4.0, 4.54415, 1);   // C central, H-H

    // H central atom (eli=2): alpha non-zero for C-* and H-C combos
    pot.set_triplet_params(2, 0, 0, 0.0, 1.0, 1);       // H central, W-W
    pot.set_triplet_params(2, 0, 1, 0.0, 1.0, 1);       // H central, W-C
    pot.set_triplet_params(2, 0, 2, 0.0, 1.0, 1);       // H central, W-H
    pot.set_triplet_params(2, 1, 0, 0.0, 1.0, 1);       // H central, C-W
    pot.set_triplet_params(2, 1, 1, 4.0, 0.33946, 1);   // H central, C-C
    pot.set_triplet_params(2, 1, 2, 4.0, 0.22006, 1);   // H central, C-H
    pot.set_triplet_params(2, 2, 0, 0.0, 1.0, 1);       // H central, H-W
    pot.set_triplet_params(2, 2, 1, 4.0, 1.0, 1);       // H central, H-C
    pot.set_triplet_params(2, 2, 2, 4.0, 1.0, 1);       // H central, H-H
}

template<bool Screening>
void Juslin<Screening>::load_parameters(const std::string& name) {
    if (name == "Juslin_JAP_98_123520_WCH") {
        load_juslin_jap_98_123520_wch(*this);
    } else {
        throw std::runtime_error("Juslin: unknown parameter set: " + name);
    }
}

using JuslinPotential = Juslin<false>;
using JuslinScreened  = Juslin<true>;

} // namespace atomistica
