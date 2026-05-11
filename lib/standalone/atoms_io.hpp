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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "../include/atomistica/core/atomic_system.hpp"
#include "../include/atomistica/config.hpp"

namespace atomistica {

// Unit conversion constant: 1 amu*(Å/fs)^2 = AMU_AFSQ_PER_EV eV
constexpr double AMU_AFSQ_PER_EV = 103.636;

// Standard atomic masses by atomic number (amu)
inline double standard_atomic_mass(int Z) {
    static const std::unordered_map<int, double> masses = {
        {1, 1.008}, {2, 4.003}, {3, 6.941}, {4, 9.012}, {5, 10.811},
        {6, 12.011}, {7, 14.007}, {8, 15.999}, {9, 18.998}, {10, 20.180},
        {11, 22.990}, {12, 24.305}, {13, 26.982}, {14, 28.086}, {15, 30.974},
        {16, 32.060}, {17, 35.453}, {18, 39.948}, {19, 39.098}, {20, 40.078},
        {26, 55.845}, {28, 58.693}, {29, 63.546}, {47, 107.868},
        {74, 183.840}, {78, 195.084}, {79, 196.967},
    };
    auto it = masses.find(Z);
    return (it != masses.end()) ? it->second : static_cast<double>(Z);
}

// Element symbol to atomic number
inline int symbol_to_Z(const std::string& sym) {
    static const std::unordered_map<std::string, int> sym2Z = {
        {"H",1},{"He",2},{"Li",3},{"Be",4},{"B",5},{"C",6},{"N",7},
        {"O",8},{"F",9},{"Ne",10},{"Na",11},{"Mg",12},{"Al",13},
        {"Si",14},{"P",15},{"S",16},{"Cl",17},{"Ar",18},{"K",19},
        {"Ca",20},{"Fe",26},{"Ni",28},{"Cu",29},{"Ag",47},
        {"W",74},{"Pt",78},{"Au",79},
    };
    auto it = sym2Z.find(sym);
    if (it != sym2Z.end()) return it->second;
    throw std::runtime_error("Unknown element symbol: " + sym);
}

// Atomic number to element symbol
inline std::string Z_to_symbol(int Z) {
    static const std::unordered_map<int, std::string> Z2sym = {
        {1,"H"},{2,"He"},{3,"Li"},{4,"Be"},{5,"B"},{6,"C"},{7,"N"},
        {8,"O"},{9,"F"},{10,"Ne"},{11,"Na"},{12,"Mg"},{13,"Al"},
        {14,"Si"},{15,"P"},{16,"S"},{17,"Cl"},{18,"Ar"},{19,"K"},
        {20,"Ca"},{26,"Fe"},{28,"Ni"},{29,"Cu"},{47,"Ag"},
        {74,"W"},{78,"Pt"},{79,"Au"},
    };
    auto it = Z2sym.find(Z);
    if (it != Z2sym.end()) return it->second;
    return "X";
}

// Replace Fortran D-notation exponents with E for C++ parsing
inline std::string fortran_to_cpp_float(const std::string& s) {
    std::string r = s;
    for (auto& c : r) {
        if (c == 'D' || c == 'd') c = 'E';
    }
    return r;
}

inline double parse_float(const std::string& s) {
    return std::stod(fortran_to_cpp_float(s));
}

// Strip whitespace from both ends
inline std::string strip(const std::string& s) {
    size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

// Case-insensitive string comparison
inline bool iequal(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (std::tolower(a[i]) != std::tolower(b[i])) return false;
    return true;
}

// Extract the section label from a "<--- ..." line
inline std::string section_label(const std::string& line) {
    std::string s = strip(line);
    // Remove leading "<---" and optional spaces/dashes
    size_t pos = 0;
    while (pos < s.size() && (s[pos] == '<' || s[pos] == '-' || s[pos] == ' '))
        ++pos;
    return strip(s.substr(pos));
}

struct AtomsData {
    AtomicSystem system;
    bool has_velocities = false;

    enum class UnitMode { eV_A, eV_A_fs } unit_mode = UnitMode::eV_A;
};

inline AtomsData read_atoms_dat(const std::string& filename) {
    std::ifstream f(filename);
    if (!f) throw std::runtime_error("Cannot open atoms file: " + filename);

    AtomsData data;
    int nat = 0;
    bool in_positions = false;
    bool in_velocities = false;
    bool in_forces = false;
    bool in_cell = false;
    bool positions_done = false;
    int pos_count = 0;
    int vel_count = 0;
    int force_count = 0;
    int cell_count = 0;

    // temporary storage
    std::vector<int> Zs;
    std::vector<int> groups;
    std::vector<std::array<double,3>> positions;
    std::vector<std::array<double,3>> velocities;
    Mat3 cell = Mat3::Identity();
    cell *= 100.0; // default large cell

    // Detect unit mode from masses: if first atom mass >> standard amu, it's eV_A_fs
    bool unit_detected = false;

    std::string line;
    while (std::getline(f, line)) {
        std::string s = strip(line);
        if (s.empty()) continue;

        // Detect section headers
        if (s.find("<---") != std::string::npos ||
            s.find("<--") != std::string::npos) {
            std::string label = section_label(s);
            in_positions = false;
            in_velocities = false;
            in_forces = false;
            in_cell = false;

            if (iequal(label, "Total number of atoms") ||
                iequal(label, "Number of atoms")) {
                // Next non-empty line is nat
                std::string nl;
                while (std::getline(f, nl)) {
                    nl = strip(nl);
                    if (!nl.empty() && nl.find("<---") == std::string::npos) {
                        nat = std::stoi(nl);
                        Zs.resize(nat);
                        groups.resize(nat, 1);
                        positions.resize(nat);
                        velocities.resize(nat, {0.0, 0.0, 0.0});
                        break;
                    }
                }
            } else if (iequal(label, "*** The following line is ignored ***") ||
                       iequal(label, "Number of occupied orbitals")) {
                // Skip next line
                std::string nl;
                std::getline(f, nl);
            } else if (iequal(label, "Element, atomic mass, coordinates, group, dissipation, temperature, (next)") ||
                       iequal(label, "Atom positions")) {
                in_positions = true;
                pos_count = 0;
            } else if (iequal(label, "Velocities")) {
                in_velocities = true;
                vel_count = 0;
                data.has_velocities = true;
            } else if (iequal(label, "Forces")) {
                in_forces = true;
                force_count = 0;
            } else if (iequal(label, "cell")) {
                in_cell = true;
                cell_count = 0;
            }
            // Other sections (shear_dx, continuous_coordinates, etc.) are skipped
            continue;
        }

        if (in_positions && pos_count < nat) {
            std::istringstream ss(s);
            std::string sym;
            double mass, x, y, z;
            ss >> sym >> mass >> x >> y >> z;
            if (!ss) throw std::runtime_error("Error parsing atom position line: " + s);

            int Z = symbol_to_Z(sym);
            Zs[pos_count] = Z;
            positions[pos_count] = {x, y, z};
            // Optional group field (6th value after coordinates; default 1)
            int grp = 1;
            ss >> grp;  // silently ignored if not present
            groups[pos_count] = grp;

            // Detect unit mode from first atom's mass
            if (!unit_detected) {
                double std_mass = standard_atomic_mass(Z);
                if (std_mass > 0.0 && mass / std_mass > 50.0) {
                    data.unit_mode = AtomsData::UnitMode::eV_A_fs;
                } else {
                    data.unit_mode = AtomsData::UnitMode::eV_A;
                }
                unit_detected = true;
            }
            ++pos_count;
            if (pos_count == nat) {
                in_positions = false;
                positions_done = true;
            }
        } else if (in_velocities && vel_count < nat) {
            std::istringstream ss(fortran_to_cpp_float(s));
            double vx, vy, vz;
            if (!(ss >> vx >> vy >> vz))
                throw std::runtime_error("Error parsing velocity line: " + s);
            velocities[vel_count] = {vx, vy, vz};
            ++vel_count;
            if (vel_count == nat) in_velocities = false;
        } else if (in_forces && force_count < nat) {
            ++force_count;
            if (force_count == nat) in_forces = false;
        } else if (in_cell && cell_count < 3) {
            std::istringstream ss(fortran_to_cpp_float(s));
            double a, b, c;
            if (!(ss >> a >> b >> c))
                throw std::runtime_error("Error parsing cell line: " + s);
            // Each row of the file is a column of the cell matrix (= lattice vector)
            cell.col(cell_count) = Vec3(a, b, c);
            ++cell_count;
            if (cell_count == 3) in_cell = false;
        }
    }

    if (nat == 0) throw std::runtime_error("No atoms found in: " + filename);

    // Build AtomicSystem
    data.system.resize(nat);
    data.system.set_cell(cell);
    data.system.pbc() = {true, true, true};

    double mass_scale = (data.unit_mode == AtomsData::UnitMode::eV_A_fs)
                        ? 1.0 / AMU_AFSQ_PER_EV : 1.0;
    double vel_scale = (data.unit_mode == AtomsData::UnitMode::eV_A_fs)
                       ? 1.0 / std::sqrt(AMU_AFSQ_PER_EV) : 1.0;

    // Store group flags in per-atom property
    if (!groups.empty()) {
        data.system.properties().add<ArrayXi>("group", static_cast<size_t>(nat));
        for (int i = 0; i < nat; ++i)
            data.system.properties().get<ArrayXi>("group")[i] = groups[i];
    }

    for (int i = 0; i < nat; ++i) {
        data.system.atomic_numbers()[i] = Zs[i];
        double m = standard_atomic_mass(Zs[i]) * mass_scale;
        // eV_A_fs: mass_scale divides out the 103.636 to get amu
        // eV_A: mass_scale = 1, already in amu
        // Actually for eV_A_fs: atoms.dat mass is amu*103.636 but we use standard mass directly
        // Just use standard amu regardless and apply mass_scale=1 always
        // (We ignore the file mass, use lookup)
        (void)mass_scale; // not used — always use standard amu
        m = standard_atomic_mass(Zs[i]);
        data.system.set_mass(i, m);

        Vec3 r(positions[i][0], positions[i][1], positions[i][2]);
        data.system.set_position(i, r);

        Vec3 v(velocities[i][0] * vel_scale,
               velocities[i][1] * vel_scale,
               velocities[i][2] * vel_scale);
        data.system.set_velocity(i, v);
    }

    return data;
}

inline void write_atoms_dat(const std::string& filename,
                             const AtomicSystem& system,
                             AtomsData::UnitMode mode = AtomsData::UnitMode::eV_A) {
    std::ofstream f(filename);
    if (!f) throw std::runtime_error("Cannot open for writing: " + filename);

    int nat = static_cast<int>(system.num_atoms());
    double mass_scale = (mode == AtomsData::UnitMode::eV_A_fs) ? AMU_AFSQ_PER_EV : 1.0;
    double vel_scale  = (mode == AtomsData::UnitMode::eV_A_fs) ? std::sqrt(AMU_AFSQ_PER_EV) : 1.0;

    f << "<--- Total number of atoms\n";
    f << " " << nat << "\n";
    f << "<--- *** The following line is ignored ***\n";
    f << " \n";
    f << "<--- Element, atomic mass, coordinates, group, dissipation, temperature, (next)\n";

    const ArrayXi* grp_arr = system.properties().has("group")
        ? &system.properties().get<ArrayXi>("group") : nullptr;

    for (int i = 0; i < nat; ++i) {
        int Z = system.atomic_number(i);
        double m = standard_atomic_mass(Z) * mass_scale;
        Vec3 r = system.position(i);
        int grp = grp_arr ? (*grp_arr)[i] : 1;
        char buf[256];
        std::snprintf(buf, sizeof(buf),
            " %-4s%20.10E%20.10E%20.10E%20.10E%5d%20.10E%20.10E\n",
            Z_to_symbol(Z).c_str(), m, r[0], r[1], r[2], grp, 0.0, 0.0);
        f << buf;
    }

    f << "<--- Velocities\n";
    for (int i = 0; i < nat; ++i) {
        Vec3 v = system.velocity(i) * vel_scale;
        char buf[128];
        std::snprintf(buf, sizeof(buf),
            " %20.10E%20.10E%20.10E\n", v[0], v[1], v[2]);
        f << buf;
    }

    f << "<--- Forces\n";
    for (int i = 0; i < nat; ++i) {
        Vec3 fv = system.forces().col(i).matrix();
        char buf[128];
        std::snprintf(buf, sizeof(buf),
            " %20.10E%20.10E%20.10E\n", fv[0], fv[1], fv[2]);
        f << buf;
    }

    f << " <--- cell\n";
    for (int j = 0; j < 3; ++j) {
        Vec3 col = system.cell().col(j);
        char buf[128];
        std::snprintf(buf, sizeof(buf),
            "     %20.10E%20.10E%20.10E\n", col[0], col[1], col[2]);
        f << buf;
    }
}

} // namespace atomistica
