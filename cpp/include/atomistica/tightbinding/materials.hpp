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
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../config.hpp"
#include "../math/spline.hpp"
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
              int n_columns) {
        n_ = x.size();
        n_cols_ = n_columns;
        x_ = x;
        y_ = y;

        // Compute second derivatives for each column using natural spline conditions
        d2y_.resize(n_, std::vector<Scalar>(n_cols_, 0.0));

        std::vector<Scalar> u(n_);
        for (int col = 0; col < n_cols_; ++col) {
            // Natural spline: second derivative is zero at boundaries
            d2y_[0][col] = 0.0;
            u[0] = 0.0;

            // Forward pass
            for (int i = 1; i < n_ - 1; ++i) {
                Scalar sig = (x_[i] - x_[i-1]) / (x_[i+1] - x_[i-1]);
                Scalar p = sig * d2y_[i-1][col] + 2.0;
                d2y_[i][col] = (sig - 1.0) / p;
                u[i] = (y_[i+1][col] - y_[i][col]) / (x_[i+1] - x_[i])
                     - (y_[i][col] - y_[i-1][col]) / (x_[i] - x_[i-1]);
                u[i] = (6.0 * u[i] / (x_[i+1] - x_[i-1]) - sig * u[i-1]) / p;
            }

            // Backward pass
            d2y_[n_-1][col] = 0.0;
            for (int i = n_ - 2; i >= 0; --i) {
                d2y_[i][col] = d2y_[i][col] * d2y_[i+1][col] + u[i];
            }
        }

        cutoff_ = x_.back();
    }

    /**
     * @brief Initialize from uniform grid
     */
    void init_uniform(Scalar x0, Scalar dx, const std::vector<std::vector<Scalar>>& y,
                      int n_columns) {
        std::vector<Scalar> x(y.size());
        for (size_t i = 0; i < y.size(); ++i) {
            x[i] = x0 + i * dx;
        }
        init(x, y, n_columns);
    }

    /**
     * @brief Evaluate spline at distance r
     *
     * @param r Distance
     * @param values Output array for interpolated values
     */
    void eval(Scalar r, std::array<Scalar, NUM_SK_INTEGRALS>& values) const {
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
    void eval_deriv(Scalar r, std::array<Scalar, NUM_SK_INTEGRALS>& values,
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

    void init(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
        n_ = x.size();
        x_ = x;
        y_ = y;

        // Compute second derivatives
        d2y_.resize(n_, 0.0);
        std::vector<Scalar> u(n_);

        d2y_[0] = 0.0;
        u[0] = 0.0;

        for (int i = 1; i < n_ - 1; ++i) {
            Scalar sig = (x_[i] - x_[i-1]) / (x_[i+1] - x_[i-1]);
            Scalar p = sig * d2y_[i-1] + 2.0;
            d2y_[i] = (sig - 1.0) / p;
            u[i] = (y_[i+1] - y_[i]) / (x_[i+1] - x_[i])
                 - (y_[i] - y_[i-1]) / (x_[i] - x_[i-1]);
            u[i] = (6.0 * u[i] / (x_[i+1] - x_[i-1]) - sig * u[i-1]) / p;
        }

        d2y_[n_-1] = 0.0;
        for (int i = n_ - 2; i >= 0; --i) {
            d2y_[i] = d2y_[i] * d2y_[i+1] + u[i];
        }

        cutoff_ = x_.back();
    }

    Scalar eval(Scalar r) const {
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

    Scalar eval_deriv(Scalar r, Scalar& derivative) const {
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
    void load_skf_directory(const std::string& path) {
        folder_ = path;
        // Ensure trailing slash
        if (!folder_.empty() && folder_.back() != '/') {
            folder_ += '/';
        }
    }

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
    void load_pair(int Z1, int Z2) {
        if (Z1 > Z2) std::swap(Z1, Z2);
        auto key = std::make_pair(Z1, Z2);

        if (pair_loaded_.find(key) != pair_loaded_.end()) {
            return;  // Already loaded
        }

        std::string sym1 = get_element_symbol(Z1);
        std::string sym2 = get_element_symbol(Z2);

        std::string filename = folder_ + sym1 + "-" + sym2 + ".skf";
        load_skf_file(filename, Z1, Z2);

        pair_loaded_[key] = true;
    }

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
    std::string get_element_symbol(int Z) const {
        static const char* symbols[] = {
            "X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
            "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
            "In", "Sn", "Sb", "Te", "I", "Xe"
        };
        if (Z > 0 && Z < 55) return symbols[Z];
        return "X";
    }

    /**
     * @brief Load SKF file in DFTB format
     */
    void load_skf_file(const std::string& filename, int Z1, int Z2) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Cannot open SKF file: " + filename);
        }

        // Read first line: dx, n
        Scalar dx;
        int n;
        file >> dx >> n;

        // Skip rest of first line (may have additional values)
        std::string line;
        std::getline(file, line);

        // For diagonal elements, read element parameters
        if (Z1 == Z2) {
            std::getline(file, line);
            std::istringstream iss(line);

            TBElementParams elem;
            elem.atomic_number = Z1;
            elem.symbol = get_element_symbol(Z1);

            // Read onsite energies (d, p, s order in file)
            std::array<Scalar, 3> e_self;
            iss >> e_self[0] >> e_self[1] >> e_self[2];

            // Skip espin
            Scalar espin;
            iss >> espin;

            // Read Hubbard U (d, p, s order)
            std::array<Scalar, 3> u;
            iss >> u[0] >> u[1] >> u[2];

            // Read valence electrons (d, p, s order)
            std::array<Scalar, 3> q;
            iss >> q[0] >> q[1] >> q[2];

            // Determine orbital configuration based on non-zero entries
            elem.num_orbitals = 0;
            elem.l_max = -1;

            // s orbital (index 0 in our arrays)
            if (std::abs(e_self[2]) > 1e-10 || std::abs(q[2]) > 0.1) {
                elem.l[0] = 0;
                elem.onsite[0] = e_self[2];  // s orbital
                elem.num_orbitals = 1;
                elem.l_max = 0;
            }

            // p orbitals (indices 1,2,3)
            if (std::abs(e_self[1]) > 1e-10 || std::abs(q[1]) > 0.1) {
                for (int i = 1; i <= 3; ++i) {
                    elem.l[i] = 1;
                    elem.onsite[i] = e_self[1];  // p orbital
                }
                elem.num_orbitals = 4;
                elem.l_max = 1;
            }

            // d orbitals (indices 4-8)
            if (std::abs(e_self[0]) > 1e-10 || std::abs(q[0]) > 0.1) {
                for (int i = 4; i <= 8; ++i) {
                    elem.l[i] = 2;
                    elem.onsite[i] = e_self[0];  // d orbital
                }
                elem.num_orbitals = 9;
                elem.l_max = 2;
            }

            // Set Hubbard U (use average or s-orbital value)
            elem.hubbard_U = u[2];  // s orbital U

            // Set valence electrons
            elem.valence_electrons = q[0] + q[1] + q[2];

            elements_[Z1] = elem;
        }

        // Read H and S tables
        std::vector<std::vector<Scalar>> H_data(n, std::vector<Scalar>(NUM_SK_INTEGRALS, 0.0));
        std::vector<std::vector<Scalar>> S_data(n, std::vector<Scalar>(NUM_SK_INTEGRALS, 0.0));
        std::vector<Scalar> r_grid(n);

        for (int i = 0; i < n; ++i) {
            r_grid[i] = (i + 1) * dx;  // SKF uses 1-indexed grid

            std::getline(file, line);
            if (line.empty()) {
                std::getline(file, line);
            }
            std::istringstream iss(line);

            // Read H integrals (10 values)
            for (int j = 0; j < NUM_SK_INTEGRALS; ++j) {
                iss >> H_data[i][j];
            }
            // Read S integrals (10 values)
            for (int j = 0; j < NUM_SK_INTEGRALS; ++j) {
                iss >> S_data[i][j];
            }
        }

        // Create splines
        auto key = std::make_pair(std::min(Z1, Z2), std::max(Z1, Z2));

        H_splines_[key].init(r_grid, H_data, NUM_SK_INTEGRALS);
        S_splines_[key].init(r_grid, S_data, NUM_SK_INTEGRALS);
        cutoffs_[key] = r_grid.back();

        // Read repulsive potential (after "Spline" keyword)
        while (std::getline(file, line)) {
            if (line.find("Spline") != std::string::npos) {
                break;
            }
        }

        if (file.good()) {
            // Read repulsive spline data
            int n_rep;
            Scalar cutoff_rep;
            file >> n_rep >> cutoff_rep;

            // Read tail coefficients
            Scalar c1, c2, c3;
            file >> c1 >> c2 >> c3;

            // Read spline segments and tabulate
            const Scalar REP_DX = 0.005;
            int n_tab = static_cast<int>(cutoff_rep / REP_DX) + 1;
            std::vector<Scalar> r_rep(n_tab);
            std::vector<Scalar> v_rep(n_tab);

            // Read all segments
            std::vector<std::array<Scalar, 8>> segments;  // x1, x2, c0, c1, c2, c3, c4, c5
            for (int seg = 0; seg < n_rep; ++seg) {
                std::array<Scalar, 8> s = {0};
                file >> s[0] >> s[1];  // x1, x2

                // Read coefficients (4 or 6)
                int n_coeff = (seg == n_rep - 1) ? 6 : 4;
                for (int c = 0; c < n_coeff; ++c) {
                    file >> s[2 + c];
                }
                segments.push_back(s);
            }

            // Tabulate repulsive potential
            for (int i = 0; i < n_tab; ++i) {
                r_rep[i] = i * REP_DX;
                Scalar r = r_rep[i];

                if (r < segments[0][0]) {
                    // Exponential tail
                    v_rep[i] = c3 + std::exp(c2 - c1 * r);
                } else {
                    // Find segment
                    v_rep[i] = 0.0;
                    for (const auto& s : segments) {
                        if (r >= s[0] && r < s[1]) {
                            Scalar dr = r - s[0];
                            v_rep[i] = s[2] + dr * (s[3] + dr * (s[4] + dr * (s[5] + dr * (s[6] + dr * s[7]))));
                            break;
                        }
                    }
                }
            }

            rep_splines_[key].init(r_rep, v_rep);
        }
    }
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
