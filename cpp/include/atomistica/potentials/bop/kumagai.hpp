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
 * @brief Kumagai pair parameters
 *
 * Pair-specific parameters for the Kumagai potential.
 * A, B, lambda1, lambda2 follow Tersoff convention (exponential form).
 */
struct KumagaiPairParams : public BOPPairParams {
    // Note: A, B, lambda (=lambda1), mu (=lambda2), r1, r2 from base
    // Distance function parameters
    Scalar alpha = 0.0;  // Exponential prefactor for distance function
    int beta = 1;        // Exponent for distance function (integer)

    // Screening parameters (only used when Screening=true)
    ScreeningParams screening;

    /**
     * @brief Precompute derived quantities
     */
    void precompute() {
        init_cutoff();
    }
};

/**
 * @brief Kumagai element parameters
 *
 * Element-specific parameters for the Kumagai potential.
 * Angular function has the form:
 *   g(cos) = c1 + (h - cos) * g0 * (1 + g1)
 * where:
 *   g0 = c2 * (h - cos) / (c3 + (h - cos)^2)
 *   g1 = c4 * exp(-c5 * (h - cos)^2)
 *
 * Bond order has the form:
 *   b = (1 + z^eta)^(-delta)
 */
struct KumagaiElementParams : public BOPElementParams {
    // Angular function parameters
    Scalar c1 = 0.0;
    Scalar c2 = 0.0;
    Scalar c3 = 1.0;
    Scalar c4 = 0.0;
    Scalar c5 = 0.0;
    Scalar h = 0.0;

    // Bond order parameters
    Scalar eta = 1.0;
    Scalar delta = 0.5;

    // Precomputed (for bond order)
    Scalar neg_delta = -0.5;  // -delta

    void precompute() {
        neg_delta = -delta;
        // We don't use beta and n from base class in Kumagai
        // Kumagai uses eta and delta instead
    }
};

/**
 * @brief Kumagai potential implementation
 *
 * The Kumagai potential is a bond-order potential similar to Tersoff but
 * with different angular and distance function forms:
 *
 * Kumagai, Izumi, Hara, Sakai, Comp. Mater. Sci. 39, 457 (2007)
 *
 * Energy: E = 0.5 * sum_{ij} fc(r_ij) * [V_R(r_ij) + b_ij * V_A(r_ij)]
 *
 * V_R(r) = A * exp(-lambda1 * r)
 * V_A(r) = -B * exp(-lambda2 * r)
 *
 * b_ij = (1 + zeta_ij^eta)^(-delta)
 *
 * zeta_ij = sum_{k != j} fc(r_ik) * g(cos_theta) * h(r_ij - r_ik)
 *
 * Angular function:
 *   g(cos) = c1 + (h - cos) * g0 * (1 + g1)
 *   g0 = c2 * (h - cos) / (c3 + (h - cos)^2)
 *   g1 = c4 * exp(-c5 * (h - cos)^2)
 *
 * Distance function:
 *   h(dr) = exp(alpha * dr^beta)
 *
 * @tparam Screening Enable screening (default: false)
 */
template<bool Screening = false>
class Kumagai : public BOPBase<Kumagai<Screening>, Screening> {
public:
    using Base = BOPBase<Kumagai<Screening>, Screening>;
    friend Base;

    Kumagai() = default;

    /**
     * @brief Add element with given atomic number and parameters
     */
    void add_element(int Z, const KumagaiElementParams& params = KumagaiElementParams{}) {
        int idx = static_cast<int>(element_params_.size());
        element_map_[Z] = idx;
        element_params_.push_back(params);
        element_params_.back().precompute();
        update_pair_count();
    }

    /**
     * @brief Set pair parameters for element pair
     */
    void set_pair_params(int Z1, int Z2, const KumagaiPairParams& params) {
        auto it1 = element_map_.find(Z1);
        auto it2 = element_map_.find(Z2);
        if (it1 == element_map_.end() || it2 == element_map_.end()) {
            throw std::runtime_error("Element not found in potential");
        }
        int ptype = pair_type(it1->second, it2->second);
        pair_params_[ptype] = params;
        pair_params_[ptype].precompute();
        update_cutoff();
    }

    /**
     * @brief Load built-in parameter set by name
     */
    void load_parameters(const std::string& name);

    /**
     * @brief Get maximum cutoff radius
     */
    Scalar cutoff() const { return max_cutoff_; }

    /**
     * @brief Get number of elements defined
     */
    int num_elements() const {
        return static_cast<int>(element_params_.size());
    }

    /**
     * @brief Get internal element index for atomic number Z
     * @return Element index, or -1 if not found
     */
    int element_index(int Z) const {
        auto it = element_map_.find(Z);
        return it != element_map_.end() ? it->second : -1;
    }

    /**
     * @brief Get pair type index for element pair
     */
    int pair_type(int eli, int elj) const {
        return pair_index(eli, elj, num_elements());
    }

    // =========================================================================
    // Required BOP interface (called by BOPKernel via CRTP)
    // =========================================================================

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
     * @brief Repulsive pair function V_R(r)
     *
     * V_R(r) = A * exp(-lambda1 * r)
     * dV_R/dr = -lambda1 * V_R
     */
    std::pair<Scalar, Scalar> repulsive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.lambda * r);
        Scalar V = p.A * exp_val;
        Scalar dV = -p.lambda * V;
        return {V, dV};
    }

    /**
     * @brief Attractive pair function V_A(r)
     *
     * V_A(r) = -B * exp(-lambda2 * r)
     * dV_A/dr = lambda2 * B * exp(-lambda2 * r) = -lambda2 * V_A
     */
    std::pair<Scalar, Scalar> attractive(int ptype, Scalar r) const {
        const auto& p = pair_params_[ptype];
        Scalar exp_val = std::exp(-p.mu * r);
        Scalar V = -p.B * exp_val;
        Scalar dV = -p.mu * V;  // = p.mu * p.B * exp_val
        return {V, dV};
    }

    /**
     * @brief Angular function g(cos_theta)
     *
     * Kumagai angular function:
     *   g(cos) = c1 + (h - cos) * g0 * (1 + g1)
     * where:
     *   g0 = c2 * (h - cos) / (c3 + (h - cos)^2)
     *   g1 = c4 * exp(-c5 * (h - cos)^2)
     *
     * Let x = h - cos, then:
     *   g = c1 + x * c2 * x / (c3 + x^2) * (1 + c4 * exp(-c5 * x^2))
     *
     * dg/d(cos) = dg/dx * dx/d(cos) = -dg/dx
     *
     * Let f = x^2 / (c3 + x^2), so g0 = c2 * f / x * x = c2 * f * x / x = c2 * x * f/x
     * Actually: g0 = c2 * x / (c3 + x^2)
     *
     * g = c1 + x * g0 * (1 + g1)
     *   = c1 + c2 * x^2 / (c3 + x^2) * (1 + c4 * exp(-c5 * x^2))
     *
     * Note: Uses element i parameters (the central atom)
     */
    std::pair<Scalar, Scalar> angular_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik, Scalar cos_theta) const
    {
        const auto& p = element_params_[eli];

        Scalar h_cos = p.h - cos_theta;
        Scalar h_cos_sq = h_cos * h_cos;

        Scalar denom = p.c3 + h_cos_sq;
        Scalar tmp = h_cos / denom;  // x / (c3 + x^2)
        Scalar g0 = p.c2 * tmp;      // c2 * x / (c3 + x^2)
        Scalar g1 = p.c4 * std::exp(-p.c5 * h_cos_sq);

        // g = c1 + h_cos * g0 * (1 + g1)
        Scalar val = g0 * (1.0 + g1);
        // dval/d(h_cos) before adding c1:
        // Let u = h_cos * g0 * (1 + g1)
        // du/d(h_cos) = g0 * (1 + g1) + h_cos * dg0/d(h_cos) * (1 + g1) + h_cos * g0 * dg1/d(h_cos)
        //
        // g0 = c2 * h_cos / (c3 + h_cos^2)
        // dg0/d(h_cos) = c2 * (c3 + h_cos^2 - h_cos * 2 * h_cos) / (c3 + h_cos^2)^2
        //              = c2 * (c3 - h_cos^2) / (c3 + h_cos^2)^2
        //
        // g1 = c4 * exp(-c5 * h_cos^2)
        // dg1/d(h_cos) = -2 * c5 * h_cos * g1
        //
        // du/d(h_cos) = val + h_cos * [dg0/d(h_cos) * (1 + g1) + g0 * dg1/d(h_cos)]
        //
        // But dg/d(cos) = -du/d(h_cos)
        //
        // From Fortran:
        // val         = go*(1.0_DP+ga1)
        // dval_dcosth = -2*(1.0_DP-h_cos*tmp)*val+2*c5*h_cos_sq*go*ga1
        // val         = c1+h_cos*val
        //
        // Let's follow Fortran exactly:
        // After "val = go*(1.0+ga1)" and before final assignment:
        // dval_dcosth = -2*(1 - h_cos*tmp)*val + 2*c5*h_cos_sq*go*ga1
        //
        // Note: tmp = h_cos/(c3+h_cos_sq), so h_cos*tmp = h_cos^2/(c3+h_cos_sq)
        Scalar dval_dcosth = -2.0 * (1.0 - h_cos * tmp) * val + 2.0 * p.c5 * h_cos_sq * g0 * g1;

        // Final value
        Scalar g = p.c1 + h_cos * val;

        return {g, dval_dcosth};
    }

    /**
     * @brief Distance function h(r_ij - r_ik)
     *
     * h(dr) = exp(alpha * dr^beta)
     * dh/d(dr) = beta * alpha * dr^(beta-1) * h
     *
     * Note: dr = r_ij - r_ik
     * dh/d(r_ik) = -dh/d(dr)
     * dh/d(r_ij) = +dh/d(dr)
     *
     * Returns (h, dh/d(r_ik), dh/d(r_ij))
     */
    std::tuple<Scalar, Scalar, Scalar> distance_function(
        int eli, int elj, int elk, int ptype_ij, int ptype_ik,
        Scalar r_ij, Scalar r_ik) const
    {
        const auto& p = pair_params_[ptype_ik];

        if (p.alpha == 0.0) {
            return {1.0, 0.0, 0.0};
        }

        Scalar dr = r_ij - r_ik;

        Scalar h, dh_dr;
        if (p.beta == 1) {
            h = std::exp(p.alpha * dr);
            dh_dr = p.alpha * h;
        } else if (p.beta == 3) {
            Scalar dr3 = dr * dr * dr;
            h = std::exp(p.alpha * dr3);
            dh_dr = 3.0 * p.alpha * dr * dr * h;
        } else {
            Scalar dr_beta = std::pow(dr, p.beta);
            h = std::exp(p.alpha * dr_beta);
            dh_dr = p.beta * p.alpha * std::pow(dr, p.beta - 1) * h;
        }

        // dh/dr_ik = -dh_dr, dh/dr_ij = +dh_dr
        return {h, -dh_dr, dh_dr};
    }

    /**
     * @brief Bond order function b(z)
     *
     * b(z) = (1 + z^eta)^(-delta)
     * db/dz = -delta * eta * z^(eta-1) * (1 + z^eta)^(-delta - 1)
     *
     * Note: Kumagai uses delta directly (not -1/(2n) like Tersoff)
     */
    std::pair<Scalar, Scalar> bond_order(int eli, int ptype, Scalar z) const {
        const auto& p = element_params_[eli];

        if (z < 1e-10) {
            return {1.0, 0.0};
        }

        Scalar z_eta = std::pow(z, p.eta);
        Scalar arg = 1.0 + z_eta;
        Scalar b = std::pow(arg, p.neg_delta);  // neg_delta = -delta

        // db/dz = neg_delta * eta * z^(eta-1) * (1 + z^eta)^(neg_delta - 1)
        Scalar db = p.neg_delta * p.eta * std::pow(z, p.eta - 1.0) * std::pow(arg, p.neg_delta - 1.0);

        return {b, db};
    }

private:
    void update_pair_count() {
        int n = num_elements();
        int np = num_pairs(n);
        if (static_cast<int>(pair_params_.size()) < np) {
            pair_params_.resize(np);
        }
    }

    void update_cutoff() {
        max_cutoff_ = 0.0;
        for (const auto& p : pair_params_) {
            if constexpr (Screening) {
                max_cutoff_ = std::max(max_cutoff_, p.screening.cut_out_h);
            } else {
                max_cutoff_ = std::max(max_cutoff_, p.r2);
            }
        }
    }

    std::map<int, int> element_map_;
    std::vector<KumagaiElementParams> element_params_;
    std::vector<KumagaiPairParams> pair_params_;
    Scalar max_cutoff_ = 0.0;
};

// =========================================================================
// Parameter loading implementations
// =========================================================================

/**
 * @brief Load Kumagai Si parameters
 *
 * Kumagai, Izumi, Hara, Sakai, Comp. Mater. Sci. 39, 457 (2007)
 */
template<bool Screening>
void load_kumagai_si(Kumagai<Screening>& pot) {
    // Element: Si (Z=14)
    KumagaiElementParams si;
    si.c1 = 0.20173476;
    si.c2 = 730418.72;
    si.c3 = 1000000.0;
    si.c4 = 1.0;
    si.c5 = 26.0;
    si.h = -0.365;
    si.eta = 1.0;
    si.delta = 0.53298909;
    si.precompute();

    pot.add_element(14, si);  // Si

    // Si-Si pair parameters
    KumagaiPairParams si_si;
    si_si.A = 3281.5905;
    si_si.B = 121.00047;
    si_si.lambda = 3.2300135;  // lambda1
    si_si.mu = 1.3457970;      // lambda2
    si_si.alpha = 2.3890327;
    si_si.beta = 1;

    if constexpr (Screening) {
        si_si.r1 = 2.50;
        si_si.r2 = 2.50 * 1.2;  // = 3.0
        si_si.screening.cut_in_l = si_si.r1;
        si_si.screening.cut_in_h = si_si.r2;
        si_si.screening.cut_out_l = 3.0;
        si_si.screening.cut_out_h = 3.0 * 2.0;  // = 6.0
        si_si.screening.cut_bo_l = 3.0;
        si_si.screening.cut_bo_h = 3.0 * 2.0;  // = 6.0
        si_si.screening.Cmin = 1.0;
        si_si.screening.Cmax = 3.0;
        si_si.screening.precompute();
    } else {
        si_si.r1 = 2.70;
        si_si.r2 = 3.30;
    }
    si_si.precompute();

    pot.set_pair_params(14, 14, si_si);
}

template<bool Screening>
void Kumagai<Screening>::load_parameters(const std::string& name) {
    element_map_.clear();
    element_params_.clear();
    pair_params_.clear();
    max_cutoff_ = 0.0;

    if (name == "Kumagai_CompMaterSci_39_457_Si") {
        load_kumagai_si(*this);
    } else {
        throw std::runtime_error("Unknown parameter set: " + name);
    }
}

// Type aliases
using KumagaiPotential = Kumagai<false>;
using KumagaiScreened = Kumagai<true>;

} // namespace atomistica
