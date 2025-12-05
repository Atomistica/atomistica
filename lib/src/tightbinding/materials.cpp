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

#include <atomistica/tightbinding/materials.hpp>

#include <fstream>
#include <sstream>

namespace atomistica {
namespace tb {

// Helper: expand Fortran repeat notation (e.g., "5*0.0" -> "0.0 0.0 0.0 0.0 0.0")
static std::string expand_fortran_notation(const std::string& input) {
    std::string result;
    std::istringstream iss(input);
    std::string token;

    while (iss >> token) {
        // Check for repeat notation: "N*value"
        size_t star_pos = token.find('*');
        if (star_pos != std::string::npos) {
            std::string count_str = token.substr(0, star_pos);
            std::string value_str = token.substr(star_pos + 1);

            // Remove trailing commas from value
            if (!value_str.empty() && value_str.back() == ',') {
                value_str.pop_back();
            }

            int count = std::stoi(count_str);
            for (int i = 0; i < count; ++i) {
                if (!result.empty()) result += " ";
                result += value_str;
            }
        } else {
            // Remove trailing commas
            if (!token.empty() && token.back() == ',') {
                token.pop_back();
            }
            if (!result.empty()) result += " ";
            result += token;
        }
    }
    return result;
}

// Helper: parse SKF line - replaces commas with spaces and expands Fortran notation
static std::string parse_skf_line(std::string line) {
    // Replace commas with spaces
    for (auto& c : line) {
        if (c == ',') c = ' ';
    }
    // Expand Fortran repeat notation
    return expand_fortran_notation(line);
}

void SKSpline::init(const std::vector<Scalar>& x, const std::vector<std::vector<Scalar>>& y,
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

void SKSpline::init_uniform(Scalar x0, Scalar dx, const std::vector<std::vector<Scalar>>& y,
                            int n_columns) {
    std::vector<Scalar> x(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
        x[i] = x0 + i * dx;
    }
    init(x, y, n_columns);
}

void RepulsiveSpline::init(const std::vector<Scalar>& x, const std::vector<Scalar>& y) {
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

void MaterialsDatabase::load_skf_directory(const std::string& path) {
    folder_ = path;
    // Ensure trailing slash
    if (!folder_.empty() && folder_.back() != '/') {
        folder_ += '/';
    }
}

void MaterialsDatabase::load_pair(int Z1, int Z2) {
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

std::string MaterialsDatabase::get_element_symbol(int Z) const {
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

void MaterialsDatabase::load_skf_file(const std::string& filename, int Z1, int Z2) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open SKF file: " + filename);
    }

    // Read first line: dx, n (comma or space separated)
    std::string line;
    std::getline(file, line);
    std::istringstream first_line(parse_skf_line(line));

    Scalar dx;
    int n;
    first_line >> dx >> n;

    // For diagonal elements, read element parameters
    if (Z1 == Z2) {
        std::getline(file, line);
        std::istringstream iss(parse_skf_line(line));

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
        std::istringstream iss(parse_skf_line(line));

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

} // namespace tb
} // namespace atomistica
