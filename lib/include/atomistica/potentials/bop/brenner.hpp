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
 * @brief Brenner pair parameters
 *
 * In the Brenner potential, angular parameters (gamma, c, d, h) are
 * stored per-pair rather than per-element (unlike Tersoff).
 */
struct BrennerPairParams : public BOPPairParams {
    // Brenner-specific parameters
    Scalar D0 = 0.0;      ///< Dimer binding energy
    Scalar r0 = 0.0;      ///< Dimer equilibrium distance
    Scalar S = 2.0;       ///< Slope of Pauling plot (must be > 1)
    Scalar beta = 0.0;    ///< Dimer stiffness

    // Angular parameters (per-pair in Brenner, unlike Tersoff)
    Scalar gamma = 0.0;   ///< Scaling factor for bond-order angular term
    Scalar c = 0.0;       ///< Angular parameter
    Scalar d = 1.0;       ///< Angular parameter (must be != 0)
    Scalar h = 0.0;       ///< Angular parameter (cosine offset)

    // Distance-dependent bond order parameters
    Scalar n = 1.0;       ///< Bond-order exponent: (1 + z^n)^(-1/(2n))
    int m = 1;            ///< Distance function exponent: exp((2*mu*dr)^m)

    // Precomputed values
    Scalar c_sq = 0.0;
    Scalar d_sq = 1.0;
    Scalar c_d = 0.0;     ///< c²/d²

    Scalar VR_f = 0.0;    ///< D0/(S-1)
    Scalar VA_f = 0.0;    ///< S*D0/(S-1)
    Scalar expR = 0.0;    ///< beta*sqrt(2*S)
    Scalar expA = 0.0;    ///< beta*sqrt(2/S)

    Scalar bo_exp = 0.0;  ///< -0.5/n
    Scalar bo_fac = 0.0;  ///< 0.5 * bo_exp * n = -0.25
    Scalar bo_exp1 = 0.0; ///< bo_exp - 1

    // Screening parameters (only used when Screening=true)
    ScreeningParams screening;

    void precompute() {
        if (S <= 1.0) {
            throw std::runtime_error("Brenner: S must be > 1");
        }
        if (d == 0.0) {
            throw std::runtime_error("Brenner: d must be != 0");
        }

        c_sq = c * c;
        d_sq = d * d;
        c_d = c_sq / d_sq;

        VR_f = D0 / (S - 1.0);
        VA_f = S * D0 / (S - 1.0);
        expR = beta * std::sqrt(2.0 * S);
        expA = beta * std::sqrt(2.0 / S);

        bo_exp = -0.5 / n;
        bo_fac = 0.5 * bo_exp * n;  // = -0.25
        bo_exp1 = bo_exp - 1.0;

        init_cutoff();
    }
};

/**
 * @brief Brenner element parameters
 *
 * In the Brenner potential, most parameters are stored per-pair.
 * Element parameters are minimal.
 */
struct BrennerElementParams : public BOPElementParams {
    // Element-specific parameters can be added here if needed
    // For now, we inherit beta and n from base (though Brenner stores these per-pair)
};

/**
 * @brief Brenner potential implementation
 *
 * The Brenner potential is a bond-order potential similar to Tersoff but
 * with different functional forms:
 *
 * Energy: E = 0.5 * sum_{ij} fc(r_ij) * [V_R(r_ij) + b_ij * V_A(r_ij)]
 *
 * V_R(r) = D0/(S-1) * exp(-beta*sqrt(2*S)*(r-r0))
 * V_A(r) = -S*D0/(S-1) * exp(-beta*sqrt(2/S)*(r-r0))
 *
 * b_ij = (1 + zeta_ij^n)^(-1/(2n))
 * zeta_ij = sum_{k != j} fc(r_ik) * g(cos_theta_jik) * h(dr_ij - dr_ik)
 *
 * g(cos) = gamma * (1 + c²/d² - c²/(d² + (h + cos)²))
 * h(dr) = exp((2*mu*dr)^m)  [or 1 if mu=0]
 *
 * Key differences from Tersoff:
 * - Angular parameters (gamma, c, d, h) are per-pair, not per-element
 * - Different exponential forms for V_R and V_A
 * - The "h" angular parameter adds to cos_theta (not subtracts)
 *
 * References:
 * - Brenner, Phys. Rev. B 42, 9458 (1990) [original]
 * - Erhart & Albe, Phys. Rev. B 71, 035211 (2005) [Si-C]
 * - Albe, Nordlund & Averback, Phys. Rev. B 65, 195124 (2002) [Pt-C]
 *
 * @tparam Screening Enable screening (default: false)
 */
template<bool Screening = false>
class Brenner : public BOPBase<Brenner<Screening>, Screening> {
public:
    using Base = BOPBase<Brenner<Screening>, Screening>;
    friend Base;

    Brenner() = default;

    /**
     * @brief Add an element to the potential
     */
    void add_element(int Z, const BrennerElementParams& params = BrennerElementParams{}) {
        int idx = static_cast<int>(element_params_.size());
        element_map_[Z] = idx;
        element_params_.push_back(params);
        update_pair_count();
    }

    /**
     * @brief Set pair parameters
     */
    void set_pair_params(int Z1, int Z2, const BrennerPairParams& params) {
        int el1 = element_index(Z1);
        int el2 = element_index(Z2);
        if (el1 < 0 || el2 < 0) return;

        int ptype = pair_type(el1, el2);
        ensure_pair_storage(ptype);

        pair_params_[ptype] = params;
        pair_params_[ptype].precompute();

        update_max_cutoff();
    }

    /**
     * @brief Load parameters from built-in database
     *
     * @param name Parameter set name (e.g., "Erhart_PRB_71_035211_SiC")
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
     * @brief Repulsive potential V_R(r) = D0/(S-1) * exp(-beta*sqrt(2*S)*(r-r0))
     */
    std::pair<Scalar, Scalar> repulsive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.expR * (r - p.r0));
        Scalar VR = p.VR_f * exp_val;
        Scalar dVR = -p.expR * VR;
        return {VR, dVR};
    }

    /**
     * @brief Attractive potential V_A(r) = -S*D0/(S-1) * exp(-beta*sqrt(2/S)*(r-r0))
     */
    std::pair<Scalar, Scalar> attractive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.expA * (r - p.r0));
        Scalar VA = -p.VA_f * exp_val;
        Scalar dVA = -p.expA * VA;  // = p.VA_f * p.expA * exp_val
        return {VA, dVA};
    }

    /**
     * @brief Angular function g(cos_theta)
     *
     * g(cos) = gamma * (1 + c²/d² - c²/(d² + (h + cos)²))
     * dg/dcos = 2 * gamma * c² * (h + cos) / (d² + (h + cos)²)²
     *
     * Note: Brenner uses (h + cos) while Tersoff uses (h - cos)
     * Also: Angular parameters come from the i-k pair in Brenner
     */
    std::pair<Scalar, Scalar> angular_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) const
    {
        // Use parameters from i-k pair (not central atom like Tersoff)
        const auto& p = pair_params_[ptype_ik];

        Scalar h_cos = p.h + cos_theta;  // Note: + not - (unlike Tersoff)
        Scalar h_cos2 = h_cos * h_cos;
        Scalar denom = p.d_sq + h_cos2;
        Scalar denom_inv = 1.0 / denom;

        Scalar g = p.gamma * (1.0 + p.c_d - p.c_sq * denom_inv);
        // dg/dcos = 2*gamma*c²*(h+cos)/(d²+(h+cos)²)²
        Scalar dg = 2.0 * p.gamma * p.c_sq * h_cos * denom_inv * denom_inv;

        return {g, dg};
    }

    /**
     * @brief Distance-dependent function h(dr_ij, dr_ik)
     *
     * h = exp((2*mu*dr)^m) where dr = r_ij - r_ik
     *
     * Special cases:
     * - mu = 0: h = 1, dh = 0
     * - m = 1: h = exp(2*mu*dr), dh = 2*mu*h
     * - m = 3: h = exp((2*mu*dr)³), dh = 6*mu*(2*mu*dr)²*h
     */
    std::tuple<Scalar, Scalar, Scalar> distance_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik,
        Scalar r_ij, Scalar r_ik) const
    {
        // Use mu and m from i-k pair
        const auto& p = pair_params_[ptype_ik];

        if (p.mu == 0.0) {
            return {1.0, 0.0, 0.0};
        }

        Scalar dr = r_ij - r_ik;
        Scalar h, dh_dr;

        if (p.m == 1) {
            Scalar arg = 2.0 * p.mu * dr;
            h = std::exp(arg);
            dh_dr = 2.0 * p.mu * h;
        } else if (p.m == 3) {
            Scalar arg = 2.0 * p.mu * dr;
            Scalar arg3 = arg * arg * arg;
            h = std::exp(arg3);
            dh_dr = 6.0 * p.mu * arg * arg * h;
        } else {
            Scalar arg = 2.0 * p.mu * dr;
            Scalar argm = std::pow(arg, p.m);
            h = std::exp(argm);
            dh_dr = 2.0 * p.mu * p.m * std::pow(arg, p.m - 1) * h;
        }

        // dh/dr_ik = -dh_dr, dh/dr_ij = +dh_dr
        return {h, -dh_dr, dh_dr};
    }

    /**
     * @brief Bond order function b(z)
     *
     * b(z) = (1 + z^n)^(-1/(2n))
     * db/dz = -1/2 * z^(n-1) * (1 + z^n)^(-1/(2n) - 1)
     *
     * For n=1: db/dz = -0.5 * (1+z)^(-1.5)
     *
     * Note: No beta prefactor like in Tersoff; that's absorbed into gamma
     */
    std::pair<Scalar, Scalar> bond_order(int eli, int ptype, Scalar z) const {
        const auto& p = pair_params_[ptype];

        if (p.n == 1.0) {
            // Simplified case: b = (1 + z)^(-0.5)
            // db/dz = -0.5 * (1 + z)^(-1.5)
            Scalar arg = 1.0 + z;
            Scalar b = std::pow(arg, p.bo_exp);  // bo_exp = -0.5
            Scalar db = p.bo_exp * std::pow(arg, p.bo_exp - 1.0);  // -0.5 * arg^(-1.5)
            return {b, db};
        }

        if (z < 1e-10) {
            // Avoid numerical issues at z=0
            return {1.0, 0.0};
        }

        // General case: b = (1 + z^n)^(-1/(2n))
        // db/dz = -1/(2n) * n * z^(n-1) * (1 + z^n)^(-1/(2n) - 1)
        //       = -0.5 * z^(n-1) * (1 + z^n)^(bo_exp - 1)
        Scalar z_n = std::pow(z, p.n);
        Scalar arg = 1.0 + z_n;
        Scalar b = std::pow(arg, p.bo_exp);
        Scalar db = -0.5 * std::pow(z, p.n - 1.0) * std::pow(arg, p.bo_exp - 1.0);

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
    std::vector<BrennerElementParams> element_params_;
    std::vector<BrennerPairParams> pair_params_;
    Scalar max_cutoff_ = 0.0;
};

// ============================================================================
// Built-in parameter sets
// ============================================================================

/**
 * @brief Erhart-Albe Si-C parameters from PRB 71, 035211 (2005)
 */
template<bool Scr>
inline void load_erhart_prb_71_035211_sic(Brenner<Scr>& pot) {
    // Add elements: C (Z=6), Si (Z=14)
    pot.add_element(6);   // C -> index 0
    pot.add_element(14);  // Si -> index 1

    // C-C pair (index 0)
    BrennerPairParams c_c;
    c_c.D0 = 6.00;
    c_c.r0 = 1.4276;
    c_c.S = 2.167;
    c_c.beta = 2.0099;
    c_c.gamma = 0.11233;
    c_c.c = 181.910;
    c_c.d = 6.28433;
    c_c.h = 0.5556;

    if constexpr (Scr) {
        c_c.mu = 1.0 / 1.4276;
        c_c.n = 1.0;
        c_c.m = 3;
        c_c.r1 = 2.00;
        c_c.r2 = 2.00 * 1.2;
        c_c.screening.cut_in_l = c_c.r1;
        c_c.screening.cut_in_h = c_c.r2;
        c_c.screening.cut_out_l = 2.00;
        c_c.screening.cut_out_h = 2.00 * 2.0;
        c_c.screening.cut_bo_l = 2.00;
        c_c.screening.cut_bo_h = 2.00 * 2.0;
        c_c.screening.Cmin = 1.0;
        c_c.screening.Cmax = 3.0;
        c_c.screening.precompute();
    } else {
        c_c.mu = 0.0;
        c_c.n = 1.0;
        c_c.m = 1;
        c_c.r1 = 1.85;
        c_c.r2 = 2.15;
    }
    pot.set_pair_params(6, 6, c_c);

    // C-Si pair (index 1)
    BrennerPairParams c_si;
    c_si.D0 = 4.36;
    c_si.r0 = 1.79;
    c_si.S = 1.847;
    c_si.beta = 1.6991;
    c_si.gamma = 0.011877;
    c_si.c = 273987.0;
    c_si.d = 180.314;
    c_si.h = 0.68;

    if constexpr (Scr) {
        c_si.mu = 1.0 / 1.79;
        c_si.n = 1.0;
        c_si.m = 3;
        c_si.r1 = std::sqrt(2.00 * 2.50);
        c_si.r2 = std::sqrt(2.00 * 2.50) * 1.2;
        c_si.screening.cut_in_l = c_si.r1;
        c_si.screening.cut_in_h = c_si.r2;
        c_si.screening.cut_out_l = std::sqrt(2.00 * 3.00);
        c_si.screening.cut_out_h = std::sqrt(2.00 * 3.00) * 2.0;
        c_si.screening.cut_bo_l = std::sqrt(2.00 * 3.00);
        c_si.screening.cut_bo_h = std::sqrt(2.00 * 3.00) * 2.0;
        c_si.screening.Cmin = 1.0;
        c_si.screening.Cmax = 3.0;
        c_si.screening.precompute();
    } else {
        c_si.mu = 0.0;
        c_si.n = 1.0;
        c_si.m = 1;
        c_si.r1 = 2.20;
        c_si.r2 = 2.60;
    }
    pot.set_pair_params(6, 14, c_si);

    // Si-Si pair (index 2)
    BrennerPairParams si_si;
    si_si.D0 = 3.24;
    si_si.r0 = 2.232;
    si_si.S = 1.842;
    si_si.beta = 1.4761;
    si_si.gamma = 0.114354;
    si_si.c = 2.00494;
    si_si.d = 0.81472;
    si_si.h = 0.259;

    if constexpr (Scr) {
        si_si.mu = 1.0 / 1.842;
        si_si.n = 1.0;
        si_si.m = 3;
        si_si.r1 = 2.50;
        si_si.r2 = 2.50 * 1.2;
        si_si.screening.cut_in_l = si_si.r1;
        si_si.screening.cut_in_h = si_si.r2;
        si_si.screening.cut_out_l = 3.00;
        si_si.screening.cut_out_h = 3.00 * 2.0;
        si_si.screening.cut_bo_l = 3.00;
        si_si.screening.cut_bo_h = 3.00 * 2.0;
        si_si.screening.Cmin = 1.0;
        si_si.screening.Cmax = 3.0;
        si_si.screening.precompute();
    } else {
        si_si.mu = 0.0;
        si_si.n = 1.0;
        si_si.m = 1;
        si_si.r1 = 2.68;
        si_si.r2 = 2.96;
    }
    pot.set_pair_params(14, 14, si_si);
}

/**
 * @brief Albe Pt-C parameters from PRB 65, 195124 (2002)
 */
template<bool Scr>
inline void load_albe_prb_65_195124_ptc(Brenner<Scr>& pot) {
    // Add elements: Pt (Z=78), C (Z=6)
    pot.add_element(78);  // Pt -> index 0
    pot.add_element(6);   // C -> index 1

    // Pt-Pt pair (index 0)
    BrennerPairParams pt_pt;
    pt_pt.D0 = 3.683;
    pt_pt.r0 = 2.384;
    pt_pt.S = 2.24297;
    pt_pt.beta = 1.64249;
    pt_pt.gamma = 0.0008542;
    pt_pt.c = 34.0;
    pt_pt.d = 1.1;
    pt_pt.h = 1.0;
    pt_pt.mu = 1.335;
    pt_pt.n = 1.0;
    pt_pt.m = 1;

    if constexpr (Scr) {
        pt_pt.r1 = 2.9;
        pt_pt.r2 = 3.3;
        pt_pt.screening.cut_in_l = pt_pt.r1;
        pt_pt.screening.cut_in_h = pt_pt.r2;
        pt_pt.screening.cut_out_l = 2.9;
        pt_pt.screening.cut_out_h = 3.3;
        pt_pt.screening.cut_bo_l = 2.9;
        pt_pt.screening.cut_bo_h = 3.3;
        pt_pt.screening.Cmin = 1.0;
        pt_pt.screening.Cmax = 3.0;
        pt_pt.screening.precompute();
    } else {
        pt_pt.r1 = 2.9;
        pt_pt.r2 = 3.3;
    }
    pot.set_pair_params(78, 78, pt_pt);

    // Pt-C pair (index 1)
    BrennerPairParams pt_c;
    pt_c.D0 = 5.3;
    pt_c.r0 = 1.84;
    pt_c.S = 1.1965;
    pt_c.beta = 1.836;
    pt_c.gamma = 0.0097;
    pt_c.c = 1.23;
    pt_c.d = 0.36;
    pt_c.h = 1.0;
    pt_c.mu = 0.0;
    pt_c.n = 1.0;
    pt_c.m = 1;

    if constexpr (Scr) {
        pt_c.r1 = 2.5;
        pt_c.r2 = 2.8;
        pt_c.screening.cut_in_l = pt_c.r1;
        pt_c.screening.cut_in_h = pt_c.r2;
        pt_c.screening.cut_out_l = 2.5;
        pt_c.screening.cut_out_h = 2.8;
        pt_c.screening.cut_bo_l = 2.5;
        pt_c.screening.cut_bo_h = 2.8;
        pt_c.screening.Cmin = 1.0;
        pt_c.screening.Cmax = 3.0;
        pt_c.screening.precompute();
    } else {
        pt_c.r1 = 2.5;
        pt_c.r2 = 2.8;
    }
    pot.set_pair_params(78, 6, pt_c);

    // C-C pair (index 2)
    BrennerPairParams c_c;
    c_c.D0 = 6.0;
    c_c.r0 = 1.39;
    c_c.S = 1.22;
    c_c.beta = 2.1;
    c_c.gamma = 0.00020813;
    c_c.c = 330.0;
    c_c.d = 3.5;
    c_c.h = 1.0;
    c_c.mu = 0.0;
    c_c.n = 1.0;
    c_c.m = 1;

    if constexpr (Scr) {
        c_c.r1 = 1.7;
        c_c.r2 = 2.0;
        c_c.screening.cut_in_l = c_c.r1;
        c_c.screening.cut_in_h = c_c.r2;
        c_c.screening.cut_out_l = 1.7;
        c_c.screening.cut_out_h = 2.0;
        c_c.screening.cut_bo_l = 1.7;
        c_c.screening.cut_bo_h = 2.0;
        c_c.screening.Cmin = 1.0;
        c_c.screening.Cmax = 3.0;
        c_c.screening.precompute();
    } else {
        c_c.r1 = 1.7;
        c_c.r2 = 2.0;
    }
    pot.set_pair_params(6, 6, c_c);
}

/**
 * @brief Henriksson Fe-C parameters from PRB 79, 144107 (2009)
 */
template<bool Scr>
inline void load_henriksson_prb_79_144107_fec(Brenner<Scr>& pot) {
    // Add elements: Fe (Z=26), C (Z=6)
    pot.add_element(26);  // Fe -> index 0
    pot.add_element(6);   // C -> index 1

    // Fe-Fe pair (index 0)
    BrennerPairParams fe_fe;
    fe_fe.D0 = 1.5;
    fe_fe.r0 = 2.29;
    fe_fe.S = 2.0693109;
    fe_fe.beta = 1.4;
    fe_fe.gamma = 0.0115751;
    fe_fe.c = 1.2898716;
    fe_fe.d = 0.3413219;
    fe_fe.h = -0.26;
    fe_fe.mu = 0.0;
    fe_fe.n = 1.0;
    fe_fe.m = 1;

    if constexpr (Scr) {
        fe_fe.r1 = 2.95;
        fe_fe.r2 = 3.35;
        fe_fe.screening.cut_in_l = fe_fe.r1;
        fe_fe.screening.cut_in_h = fe_fe.r2;
        fe_fe.screening.cut_out_l = 100.0;  // Effectively disabled
        fe_fe.screening.cut_out_h = 3.35;
        fe_fe.screening.cut_bo_l = 100.0;
        fe_fe.screening.cut_bo_h = 3.35;
        fe_fe.screening.Cmin = 1.0;
        fe_fe.screening.Cmax = 3.0;
        fe_fe.screening.precompute();
    } else {
        fe_fe.r1 = 2.95;
        fe_fe.r2 = 3.35;
    }
    pot.set_pair_params(26, 26, fe_fe);

    // Fe-C pair (index 1)
    BrennerPairParams fe_c;
    fe_c.D0 = 4.82645134;
    fe_c.r0 = 1.47736510;
    fe_c.S = 1.43134755;
    fe_c.beta = 1.63208170;
    fe_c.gamma = 0.00205862;
    fe_c.c = 8.95583221;
    fe_c.d = 0.72062047;
    fe_c.h = 0.87099874;
    fe_c.mu = 0.0;
    fe_c.n = 1.0;
    fe_c.m = 1;

    if constexpr (Scr) {
        fe_c.r1 = 2.3;
        fe_c.r2 = 2.7;
        fe_c.screening.cut_in_l = fe_c.r1;
        fe_c.screening.cut_in_h = fe_c.r2;
        fe_c.screening.cut_out_l = 100.0;
        fe_c.screening.cut_out_h = 2.7;
        fe_c.screening.cut_bo_l = 100.0;
        fe_c.screening.cut_bo_h = 2.7;
        fe_c.screening.Cmin = 1.0;
        fe_c.screening.Cmax = 3.0;
        fe_c.screening.precompute();
    } else {
        fe_c.r1 = 2.3;
        fe_c.r2 = 2.7;
    }
    pot.set_pair_params(26, 6, fe_c);

    // C-C pair (index 2)
    BrennerPairParams c_c;
    c_c.D0 = 6.0;
    c_c.r0 = 1.39;
    c_c.S = 1.22;
    c_c.beta = 2.1;
    c_c.gamma = 0.00020813;
    c_c.c = 330.0;
    c_c.d = 3.5;
    c_c.h = 1.0;

    if constexpr (Scr) {
        c_c.mu = 1.0 / 1.315;
        c_c.n = 1.0;
        c_c.m = 3;
        c_c.r1 = 2.00;
        c_c.r2 = 1.2 * 2.00;
        c_c.screening.cut_in_l = c_c.r1;
        c_c.screening.cut_in_h = c_c.r2;
        c_c.screening.cut_out_l = 2.00;
        c_c.screening.cut_out_h = 2.0 * 2.00;
        c_c.screening.cut_bo_l = 1.20;
        c_c.screening.cut_bo_h = 2.0 * 2.00;
        c_c.screening.Cmin = 1.0;
        c_c.screening.Cmax = 3.0;
        c_c.screening.precompute();
    } else {
        c_c.mu = 0.0;
        c_c.n = 1.0;
        c_c.m = 1;
        c_c.r1 = 1.70;
        c_c.r2 = 2.00;
    }
    pot.set_pair_params(6, 6, c_c);
}

/**
 * @brief Kioseoglou Al-N parameters from Phys. Stat. Sol. (b) 245, 1118 (2008)
 */
template<bool Scr>
inline void load_kioseoglou_pssb_245_1118_aln(Brenner<Scr>& pot) {
    // Add elements: N (Z=7), Al (Z=13)
    pot.add_element(7);   // N -> index 0
    pot.add_element(13);  // Al -> index 1

    // N-N pair (index 0)
    BrennerPairParams n_n;
    n_n.D0 = 9.9100;
    n_n.r0 = 1.1100;
    n_n.S = 1.4922;
    n_n.beta = 2.05945;
    n_n.gamma = 0.76612;
    n_n.c = 0.178493;
    n_n.d = 0.20172;
    n_n.h = 0.045238;
    n_n.mu = 0.0;
    n_n.n = 1.0;
    n_n.m = 1;

    if constexpr (Scr) {
        n_n.r1 = 2.00;
        n_n.r2 = 2.40;
        n_n.screening.cut_in_l = n_n.r1;
        n_n.screening.cut_in_h = n_n.r2;
        n_n.screening.cut_out_l = 2.00;
        n_n.screening.cut_out_h = 2.40;
        n_n.screening.cut_bo_l = 2.00;
        n_n.screening.cut_bo_h = 2.40;
        n_n.screening.Cmin = 1.0;
        n_n.screening.Cmax = 3.0;
        n_n.screening.precompute();
    } else {
        n_n.r1 = 2.00;
        n_n.r2 = 2.40;
    }
    pot.set_pair_params(7, 7, n_n);

    // N-Al pair (index 1)
    BrennerPairParams n_al;
    n_al.D0 = 3.3407;
    n_al.r0 = 1.8616;
    n_al.S = 1.7269;
    n_al.beta = 1.7219;
    n_al.gamma = 1.1e-6;
    n_al.c = 100390.0;
    n_al.d = 16.2170;
    n_al.h = 0.5980;
    n_al.mu = 0.0;
    n_al.n = 0.7200;
    n_al.m = 1;

    if constexpr (Scr) {
        n_al.r1 = 2.19;
        n_al.r2 = 2.49;
        n_al.screening.cut_in_l = n_al.r1;
        n_al.screening.cut_in_h = n_al.r2;
        n_al.screening.cut_out_l = 2.19;
        n_al.screening.cut_out_h = 2.49;
        n_al.screening.cut_bo_l = 2.19;
        n_al.screening.cut_bo_h = 2.49;
        n_al.screening.Cmin = 1.0;
        n_al.screening.Cmax = 3.0;
        n_al.screening.precompute();
    } else {
        n_al.r1 = 2.19;
        n_al.r2 = 2.49;
    }
    pot.set_pair_params(7, 13, n_al);

    // Al-Al pair (index 2)
    BrennerPairParams al_al;
    al_al.D0 = 1.5000;
    al_al.r0 = 2.4660;
    al_al.S = 2.7876;
    al_al.beta = 1.0949;
    al_al.gamma = 0.3168;
    al_al.c = 0.0748;
    al_al.d = 19.5691;
    al_al.h = 0.6593;
    al_al.mu = 0.0;
    al_al.n = 6.0865;
    al_al.m = 1;

    if constexpr (Scr) {
        al_al.r1 = 2.60;
        al_al.r2 = 2.80;
        al_al.screening.cut_in_l = al_al.r1;
        al_al.screening.cut_in_h = al_al.r2;
        al_al.screening.cut_out_l = 3.40;
        al_al.screening.cut_out_h = 3.60;
        al_al.screening.cut_bo_l = 3.40;
        al_al.screening.cut_bo_h = 3.60;
        al_al.screening.Cmin = 1.0;
        al_al.screening.Cmax = 3.0;
        al_al.screening.precompute();
    } else {
        al_al.r1 = 3.40;
        al_al.r2 = 3.60;
    }
    pot.set_pair_params(13, 13, al_al);
}

// ============================================================================
// Brenner original C potentials (PRB 42, 9458, 1990)
// ============================================================================

/**
 * @brief Brenner original C potential, parameter set I
 *
 * D. Brenner, Phys. Rev. B 42, 9458 (1990) - potential I
 * n = 1/(2*0.80469) ≈ 0.62135
 */
template<bool Scr>
inline void load_brenner_prb_42_9458_c_i(Brenner<Scr>& pot) {
    pot.add_element(6);  // C -> index 0

    BrennerPairParams cc;
    cc.D0    = 6.325;
    cc.r0    = 1.315;
    cc.S     = 1.29;
    cc.beta  = 1.5;
    cc.gamma = 0.011304;
    cc.c     = 19.0;
    cc.d     = 2.5;
    cc.h     = 1.0;
    cc.mu    = 0.0;
    cc.n     = 1.0 / (2.0 * 0.80469);  // ≈ 0.62135
    cc.m     = 1;
    cc.r1    = 1.70;
    cc.r2    = 2.00;
    if constexpr (Scr) {
        cc.screening.cut_in_l  = 1.70; cc.screening.cut_in_h  = 2.00;
        cc.screening.cut_out_l = 1.70; cc.screening.cut_out_h = 4.00;
        cc.screening.cut_bo_l  = 1.70; cc.screening.cut_bo_h  = 4.00;
        cc.screening.Cmin = 1.0; cc.screening.Cmax = 3.0;
        cc.screening.precompute();
    }
    pot.set_pair_params(6, 6, cc);
}

/**
 * @brief Brenner original C potential, parameter set II
 *
 * D. Brenner, Phys. Rev. B 42, 9458 (1990) - potential II
 * n = 1/(2*0.5) = 1.0  (same C-C as Erhart SiC C-C component)
 */
template<bool Scr>
inline void load_brenner_prb_42_9458_c_ii(Brenner<Scr>& pot) {
    pot.add_element(6);  // C -> index 0

    BrennerPairParams cc;
    cc.D0    = 6.0;
    cc.r0    = 1.39;
    cc.S     = 1.22;
    cc.beta  = 2.1;
    cc.gamma = 0.00020813;
    cc.c     = 330.0;
    cc.d     = 3.5;
    cc.h     = 1.0;
    cc.mu    = 0.0;
    cc.n     = 1.0 / (2.0 * 0.5);  // = 1.0
    cc.m     = 1;
    cc.r1    = 1.70;
    cc.r2    = 2.00;
    if constexpr (Scr) {
        cc.screening.cut_in_l  = 1.70; cc.screening.cut_in_h  = 2.00;
        cc.screening.cut_out_l = 1.70; cc.screening.cut_out_h = 4.00;
        cc.screening.cut_bo_l  = 1.70; cc.screening.cut_bo_h  = 4.00;
        cc.screening.Cmin = 1.0; cc.screening.Cmax = 3.0;
        cc.screening.precompute();
    }
    pot.set_pair_params(6, 6, cc);
}

template<bool Screening>
void Brenner<Screening>::load_parameters(const std::string& name) {
    if (name == "Erhart_PRB_71_035211_SiC") {
        load_erhart_prb_71_035211_sic(*this);
    } else if (name == "Albe_PRB_65_195124_PtC") {
        load_albe_prb_65_195124_ptc(*this);
    } else if (name == "Henriksson_PRB_79_144107_FeC") {
        load_henriksson_prb_79_144107_fec(*this);
    } else if (name == "Kioseoglou_PSSb_245_1118_AlN") {
        load_kioseoglou_pssb_245_1118_aln(*this);
    } else if (name == "Brenner_PRB_42_9458_C_I") {
        load_brenner_prb_42_9458_c_i(*this);
    } else if (name == "Brenner_PRB_42_9458_C_II") {
        load_brenner_prb_42_9458_c_ii(*this);
    } else {
        throw std::runtime_error("Unknown parameter set: " + name);
    }
}

// Type aliases
using BrennerPotential = Brenner<false>;
using BrennerScreened = Brenner<true>;

} // namespace atomistica
