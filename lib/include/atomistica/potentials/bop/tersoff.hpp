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

template<bool Screening>
void Tersoff<Screening>::load_parameters(const std::string& name) {
    if (name == "Tersoff_PRB_39_5566_Si_C") {
        load_tersoff_prb_39_5566_si_c(*this);
    } else {
        throw std::runtime_error("Unknown parameter set: " + name);
    }
}

// Type aliases
using TersoffPotential = Tersoff<false>;
using TersoffScreened = Tersoff<true>;

} // namespace atomistica
