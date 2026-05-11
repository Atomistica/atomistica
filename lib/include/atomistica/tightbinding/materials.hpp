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

#include "../config.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Cubic spline for SK integral interpolation
 */
class SKSpline {
public:
    SKSpline() = default;

    /**
     * @brief Initialize spline from tabulated data
     *
     * @param x Distance grid points
     * @param y Function values at each grid point [n_points x n_columns]
     * @param n_columns Number of columns (e.g., 10 for SK integrals)
     */
    void init(const std::vector<Scalar>& x, const std::vector<std::vector<Scalar>>& y,
              int n_columns);

    /**
     * @brief Initialize from uniform grid
     */
    void init_uniform(Scalar x0, Scalar dx, const std::vector<std::vector<Scalar>>& y,
                      int n_columns);

    /**
     * @brief Evaluate spline at distance r
     *
     * @param r Distance
     * @param values Output array for interpolated values
     */
    inline void eval(Scalar r, std::array<Scalar, NUM_SK_INTEGRALS>& values) const {
        if (r >= cutoff_ || n_ < 2) {
            values.fill(0.0);
            return;
        }

        // Binary search for interval
        int lo = 0, hi = n_ - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x_[mid] > r) hi = mid;
            else lo = mid;
        }

        Scalar dx = x_[hi] - x_[lo];
        Scalar a = (x_[hi] - r) / dx;
        Scalar b = (r - x_[lo]) / dx;

        for (int col = 0; col < n_cols_ && col < NUM_SK_INTEGRALS; ++col) {
            values[col] = a * y_[lo][col] + b * y_[hi][col]
                        + ((a*a*a - a) * d2y_[lo][col] + (b*b*b - b) * d2y_[hi][col])
                          * dx * dx / 6.0;
        }
    }

    /**
     * @brief Evaluate spline and derivative at distance r
     */
    inline void eval_deriv(Scalar r, std::array<Scalar, NUM_SK_INTEGRALS>& values,
                           std::array<Scalar, NUM_SK_INTEGRALS>& derivatives) const {
        if (r >= cutoff_ || n_ < 2) {
            values.fill(0.0);
            derivatives.fill(0.0);
            return;
        }

        // Binary search for interval
        int lo = 0, hi = n_ - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x_[mid] > r) hi = mid;
            else lo = mid;
        }

        Scalar dx = x_[hi] - x_[lo];
        Scalar a = (x_[hi] - r) / dx;
        Scalar b = (r - x_[lo]) / dx;

        for (int col = 0; col < n_cols_ && col < NUM_SK_INTEGRALS; ++col) {
            values[col] = a * y_[lo][col] + b * y_[hi][col]
                        + ((a*a*a - a) * d2y_[lo][col] + (b*b*b - b) * d2y_[hi][col])
                          * dx * dx / 6.0;

            // Derivative
            derivatives[col] = (y_[hi][col] - y_[lo][col]) / dx
                             - (3.0*a*a - 1.0) / 6.0 * dx * d2y_[lo][col]
                             + (3.0*b*b - 1.0) / 6.0 * dx * d2y_[hi][col];
        }
    }

    Scalar cutoff() const { return cutoff_; }
    bool is_valid() const { return n_ > 1; }

private:
    int n_ = 0;
    int n_cols_ = 0;
    Scalar cutoff_ = 0.0;
    std::vector<Scalar> x_;
    std::vector<std::vector<Scalar>> y_;
    std::vector<std::vector<Scalar>> d2y_;
};

/**
 * @brief Simple 1D spline for repulsive potential
 */
class RepulsiveSpline {
public:
    RepulsiveSpline() = default;

    void init(const std::vector<Scalar>& x, const std::vector<Scalar>& y);

    inline Scalar eval(Scalar r) const {
        if (r >= cutoff_ || n_ < 2) return 0.0;

        int lo = 0, hi = n_ - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x_[mid] > r) hi = mid;
            else lo = mid;
        }

        Scalar dx = x_[hi] - x_[lo];
        Scalar a = (x_[hi] - r) / dx;
        Scalar b = (r - x_[lo]) / dx;

        return a * y_[lo] + b * y_[hi]
             + ((a*a*a - a) * d2y_[lo] + (b*b*b - b) * d2y_[hi]) * dx * dx / 6.0;
    }

    inline Scalar eval_deriv(Scalar r, Scalar& derivative) const {
        if (r >= cutoff_ || n_ < 2) {
            derivative = 0.0;
            return 0.0;
        }

        int lo = 0, hi = n_ - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (x_[mid] > r) hi = mid;
            else lo = mid;
        }

        Scalar dx = x_[hi] - x_[lo];
        Scalar a = (x_[hi] - r) / dx;
        Scalar b = (r - x_[lo]) / dx;

        Scalar val = a * y_[lo] + b * y_[hi]
                   + ((a*a*a - a) * d2y_[lo] + (b*b*b - b) * d2y_[hi]) * dx * dx / 6.0;

        derivative = (y_[hi] - y_[lo]) / dx
                   - (3.0*a*a - 1.0) / 6.0 * dx * d2y_[lo]
                   + (3.0*b*b - 1.0) / 6.0 * dx * d2y_[hi];

        return val;
    }

    Scalar cutoff() const { return cutoff_; }
    bool is_valid() const { return n_ > 1; }

private:
    int n_ = 0;
    Scalar cutoff_ = 0.0;
    std::vector<Scalar> x_;
    std::vector<Scalar> y_;
    std::vector<Scalar> d2y_;
};

/**
 * @brief Materials database for tight-binding calculations
 *
 * Stores element parameters and pair interactions from SKF files
 */
class MaterialsDatabase {
public:
    MaterialsDatabase() = default;

    /**
     * @brief Load SKF files from a directory
     *
     * @param path Directory containing SKF files (e.g., "mio-1-1/")
     */
    void load_skf_directory(const std::string& path);

    /**
     * @brief Add element parameters manually
     */
    void add_element(const TBElementParams& elem) {
        elements_[elem.atomic_number] = elem;
    }

    /**
     * @brief Get element parameters by atomic number
     */
    const TBElementParams& get_element(int Z) const {
        auto it = elements_.find(Z);
        if (it == elements_.end()) {
            throw std::runtime_error("Element Z=" + std::to_string(Z) + " not in database");
        }
        return it->second;
    }

    /**
     * @brief Check if element exists
     */
    bool has_element(int Z) const {
        return elements_.find(Z) != elements_.end();
    }

    /**
     * @brief Load pair parameters from SKF file
     *
     * @param Z1 Atomic number of first element
     * @param Z2 Atomic number of second element
     */
    void load_pair(int Z1, int Z2);

    /**
     * @brief Get H spline for pair
     */
    const SKSpline& get_H_spline(int Z1, int Z2) const {
        if (Z1 > Z2) std::swap(Z1, Z2);
        auto key = std::make_pair(Z1, Z2);
        auto it = H_splines_.find(key);
        if (it == H_splines_.end()) {
            throw std::runtime_error("H spline not loaded for pair");
        }
        return it->second;
    }

    /**
     * @brief Get S spline for pair
     */
    const SKSpline& get_S_spline(int Z1, int Z2) const {
        if (Z1 > Z2) std::swap(Z1, Z2);
        auto key = std::make_pair(Z1, Z2);
        auto it = S_splines_.find(key);
        if (it == S_splines_.end()) {
            throw std::runtime_error("S spline not loaded for pair");
        }
        return it->second;
    }

    /**
     * @brief Get repulsive spline for pair
     */
    const RepulsiveSpline& get_rep_spline(int Z1, int Z2) const {
        if (Z1 > Z2) std::swap(Z1, Z2);
        auto key = std::make_pair(Z1, Z2);
        auto it = rep_splines_.find(key);
        if (it == rep_splines_.end()) {
            throw std::runtime_error("Repulsive spline not loaded for pair");
        }
        return it->second;
    }

    /**
     * @brief Get cutoff for pair
     */
    Scalar get_cutoff(int Z1, int Z2) const {
        if (Z1 > Z2) std::swap(Z1, Z2);
        auto key = std::make_pair(Z1, Z2);
        auto it = cutoffs_.find(key);
        if (it == cutoffs_.end()) {
            return 0.0;
        }
        return it->second;
    }

    /**
     * @brief Get maximum cutoff across all loaded pairs
     */
    Scalar get_max_cutoff() const {
        Scalar max_cut = 0.0;
        for (const auto& kv : cutoffs_) {
            max_cut = std::max(max_cut, kv.second);
        }
        return max_cut;
    }

private:
    std::string folder_;
    std::map<int, TBElementParams> elements_;
    std::map<std::pair<int,int>, SKSpline> H_splines_;
    std::map<std::pair<int,int>, SKSpline> S_splines_;
    std::map<std::pair<int,int>, RepulsiveSpline> rep_splines_;
    std::map<std::pair<int,int>, Scalar> cutoffs_;
    std::map<std::pair<int,int>, bool> pair_loaded_;

    /**
     * @brief Get element symbol from atomic number
     */
    std::string get_element_symbol(int Z) const;

    /**
     * @brief Load SKF file in DFTB format
     */
    void load_skf_file(const std::string& filename, int Z1, int Z2);
};

/**
 * @brief Predefined parameter sets
 */
namespace parameters {

/**
 * @brief Set up carbon parameters for mio-1-1 DFTB
 */
inline TBElementParams carbon_mio() {
    TBElementParams c;
    c.symbol = "C";
    c.atomic_number = 6;
    c.num_orbitals = 4;  // sp
    c.l_max = 1;

    c.l = {0, 1, 1, 1, -1, -1, -1, -1, -1};

    // On-site energies from mio-1-1
    c.onsite[0] = -0.50489172;  // s
    c.onsite[1] = -0.19435511;  // p
    c.onsite[2] = -0.19435511;  // p
    c.onsite[3] = -0.19435511;  // p

    c.hubbard_U = 0.3647;
    c.valence_electrons = 4.0;

    return c;
}

/**
 * @brief Set up hydrogen parameters for mio-1-1 DFTB
 */
inline TBElementParams hydrogen_mio() {
    TBElementParams h;
    h.symbol = "H";
    h.atomic_number = 1;
    h.num_orbitals = 1;  // s only
    h.l_max = 0;

    h.l = {0, -1, -1, -1, -1, -1, -1, -1, -1};
    h.onsite[0] = -0.23855330;  // s

    h.hubbard_U = 0.4195;
    h.valence_electrons = 1.0;

    return h;
}

/**
 * @brief Set up oxygen parameters for mio-1-1 DFTB
 */
inline TBElementParams oxygen_mio() {
    TBElementParams o;
    o.symbol = "O";
    o.atomic_number = 8;
    o.num_orbitals = 4;  // sp
    o.l_max = 1;

    o.l = {0, 1, 1, 1, -1, -1, -1, -1, -1};
    o.onsite[0] = -0.87841725;  // s
    o.onsite[1] = -0.33231360;  // p
    o.onsite[2] = -0.33231360;  // p
    o.onsite[3] = -0.33231360;  // p

    o.hubbard_U = 0.4954;
    o.valence_electrons = 6.0;

    return o;
}

/**
 * @brief Set up nitrogen parameters for mio-1-1 DFTB
 */
inline TBElementParams nitrogen_mio() {
    TBElementParams n;
    n.symbol = "N";
    n.atomic_number = 7;
    n.num_orbitals = 4;  // sp
    n.l_max = 1;

    n.l = {0, 1, 1, 1, -1, -1, -1, -1, -1};
    n.onsite[0] = -0.63622628;  // s
    n.onsite[1] = -0.26004355;  // p
    n.onsite[2] = -0.26004355;  // p
    n.onsite[3] = -0.26004355;  // p

    n.hubbard_U = 0.4310;
    n.valence_electrons = 5.0;

    return n;
}

} // namespace parameters

} // namespace tb
} // namespace atomistica
