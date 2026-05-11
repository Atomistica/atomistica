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
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../../config.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../../math/spline.hpp"
#include "../potential_base.hpp"

namespace atomistica {

/**
 * @brief Element information from EAM file
 */
struct EAMElementInfo {
    std::string symbol;     // Element symbol
    int atomic_number = 0;  // Atomic number Z
    Scalar mass = 0.0;      // Atomic mass
    Scalar lattice_constant = 0.0;  // Lattice constant a0
    std::string lattice_type;       // Lattice type (FCC, BCC, etc.)
};

/**
 * @brief Single-element EAM potential (funcfl format)
 *
 * Implements the Embedded Atom Method:
 *   E = sum_i [ F_i(rho_i) + 0.5 * sum_{j!=i} phi(r_ij) ]
 *
 * where:
 *   F(rho) = embedding energy
 *   rho_i = sum_{j!=i} rho(r_ij) = electron density at atom i
 *   phi(r) = pair repulsive potential = Z(r)^2 / r
 *
 * File format (funcfl):
 *   Line 1: Comment
 *   Line 2: Z, mass, lattice_constant, lattice_type
 *   Line 3: nF, dF, nr, dr, cutoff
 *   Lines 4+: F(rho) values (nF total)
 *   Next: Z(r) values (nr total)
 *   Next: rho(r) values (nr total)
 *
 * Reference: Foiles, Baskes, Daw, PRB 33, 7983 (1986)
 */
class TabulatedEAM : public PotentialBase<TabulatedEAM> {
public:
    TabulatedEAM() = default;

    /**
     * @brief Load potential from file
     */
    void load(const std::string& filename);

    /**
     * @brief Check if potential is loaded
     */
    bool is_valid() const { return embedding_.is_valid(); }

    /**
     * @brief Get element information
     */
    const EAMElementInfo& element_info() const { return element_info_; }

    /**
     * @brief Evaluate embedding function F(rho)
     */
    SplineResult embedding(Scalar rho) const { return embedding_.eval(rho); }

    /**
     * @brief Evaluate effective charge Z(r)
     */
    SplineResult effective_charge(Scalar r) const { return Z_.eval(r); }

    /**
     * @brief Evaluate electron density rho(r)
     */
    SplineResult density(Scalar r) const { return rho_.eval(r); }

    /**
     * @brief Evaluate pair potential phi(r) = Z(r)^2 / r
     *
     * Also returns derivative: d(phi)/dr = (2*Z*dZ - phi) / r
     */
    SplineResult pair_potential(Scalar r) const;

    // CRTP implementation
    Scalar cutoff_impl() const { return cutoff_; }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    EAMElementInfo element_info_;
    Scalar cutoff_ = 0.0;

    CubicSpline embedding_;  // F(rho)
    CubicSpline Z_;          // Z(r) - effective charge
    CubicSpline rho_;        // rho(r) - electron density
};

/**
 * @brief Multi-element EAM potential (setfl/alloy format)
 *
 * Supports multiple element types with element-specific functions:
 *   E = sum_i [ F_i(rho_i) + 0.5 * sum_{j!=i} phi_{ij}(r_ij) ]
 *
 * where:
 *   F_i(rho) = embedding energy for element i
 *   rho_i = sum_{j!=i} rho_j(r_ij) = electron density at atom i
 *   phi_{ij}(r) = pair potential between elements i and j
 *
 * File format (setfl):
 *   Lines 1-3: Comments
 *   Line 4: nel, element_1, element_2, ...
 *   Line 5: nF, dF, nr, dr, cutoff
 *   For each element i:
 *     Line: Z, mass, a0, lattice
 *     F_i(rho) values (nF total)
 *     rho_i(r) values (nr total)
 *   For each pair (i,j) with i >= j:
 *     phi_{ij}(r) values (nr total)
 */
class TabulatedAlloyEAM : public PotentialBase<TabulatedAlloyEAM> {
public:
    TabulatedAlloyEAM() = default;

    /**
     * @brief Load potential from file
     */
    void load(const std::string& filename);

    /**
     * @brief Check if potential is loaded
     */
    bool is_valid() const { return num_elements_ > 0; }

    /**
     * @brief Get number of elements
     */
    int num_elements() const { return num_elements_; }

    /**
     * @brief Get element information
     */
    const EAMElementInfo& element_info(int elem_idx) const {
        return element_info_.at(elem_idx);
    }

    /**
     * @brief Get element index from symbol
     * @return Element index or -1 if not found
     */
    int element_index(const std::string& symbol) const;

    /**
     * @brief Get element index from atomic number
     * @return Element index or -1 if not found
     */
    int element_index_by_Z(int Z) const;

    /**
     * @brief Get list of element symbols
     */
    std::vector<std::string> element_symbols() const;

    /**
     * @brief Evaluate embedding function F_i(rho) for element i
     */
    SplineResult embedding(int elem_idx, Scalar rho) const {
        return embedding_.at(elem_idx).eval(rho);
    }

    /**
     * @brief Evaluate electron density rho_i(r) for element i
     */
    SplineResult density(int elem_idx, Scalar r) const {
        return rho_.at(elem_idx).eval(r);
    }

    /**
     * @brief Evaluate pair potential phi_{ij}(r) between elements i and j
     */
    SplineResult pair_potential(int elem_i, int elem_j, Scalar r) const;

    // CRTP implementation
    Scalar cutoff_impl() const { return cutoff_; }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    int num_elements_ = 0;
    Scalar cutoff_ = 0.0;

    std::vector<EAMElementInfo> element_info_;
    std::map<std::string, int> symbol_to_index_;
    std::map<int, int> Z_to_index_;

    std::vector<CubicSpline> embedding_;  // F_i(rho) per element
    std::vector<CubicSpline> rho_;        // rho_i(r) per element
    std::vector<std::vector<CubicSpline>> phi_;  // phi_{ij}(r) pair potentials [i][j]
};

// ============================================================================
// Implementation
// ============================================================================

inline void TabulatedEAM::load(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open EAM file: " + filename);
    }

    std::string line;

    // Line 1: Comment
    std::getline(file, line);

    // Line 2: Z, mass, a0, lattice
    std::getline(file, line);
    {
        std::istringstream iss(line);
        iss >> element_info_.atomic_number >> element_info_.mass
            >> element_info_.lattice_constant >> element_info_.lattice_type;
    }

    // Line 3: nF, dF, nr, dr, cutoff
    int nF, nr;
    Scalar dF, dr;
    std::getline(file, line);
    {
        std::istringstream iss(line);
        iss >> nF >> dF >> nr >> dr >> cutoff_;
    }

    // Read F(rho) values
    std::vector<Scalar> F_values(nF);
    for (int i = 0; i < nF; ++i) {
        file >> F_values[i];
    }
    embedding_.init(0.0, static_cast<Scalar>(nF - 1) * dF, F_values);

    // Read Z(r) values and apply scaling
    // The funcfl format stores Z(r) in sqrt(Hartree*Bohr) units.
    // Pair potential: phi(r) = Z_i(r)*Z_j(r)/r in Hartree-Bohr units.
    // Converting to eV/Å: phi_eV = Z_file^2 * Hartree * Bohr / r_Å
    // So Z_scaled = Z_file * sqrt(Hartree * Bohr).
    // Hartree = 27.2114 eV, Bohr = 0.529177 Å
    const Scalar Z_scale = std::sqrt(27.2114 * 0.529177);
    std::vector<Scalar> Z_values(nr);
    for (int i = 0; i < nr; ++i) {
        file >> Z_values[i];
        Z_values[i] *= Z_scale;
    }
    Z_.init(0.0, static_cast<Scalar>(nr - 1) * dr, Z_values);

    // Read rho(r) values
    std::vector<Scalar> rho_values(nr);
    for (int i = 0; i < nr; ++i) {
        file >> rho_values[i];
    }
    rho_.init(0.0, static_cast<Scalar>(nr - 1) * dr, rho_values);
}

inline SplineResult TabulatedEAM::pair_potential(Scalar r) const {
    if (r < 1e-10) {
        return {0.0, 0.0};
    }
    SplineResult Z_result = Z_.eval(r);
    Scalar Z_val = Z_result.value;
    Scalar dZ = Z_result.derivative;

    Scalar phi = Z_val * Z_val / r;
    Scalar dphi = (2.0 * Z_val * dZ - phi) / r;

    return {phi, dphi};
}

inline PotentialResults TabulatedEAM::compute_impl(
    AtomicSystem& system,
    NeighborList& neighbors,
    bool compute_forces,
    bool compute_virial)
{
    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();
    const Mat3& cell = system.cell();
    const Scalar cutoff_sq = cutoff_ * cutoff_;

    // First pass: compute electron densities at each atom
    std::vector<Scalar> rho_sum(num_atoms, 0.0);

    for (std::size_t i = 0; i < num_atoms; ++i) {
        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq < cutoff_sq && r_sq > 1e-10) {
                Scalar r = std::sqrt(r_sq);
                rho_sum[i] += rho_.eval(r).value;
            }
        }
    }

    // Second pass: compute embedding energies and store derivatives
    std::vector<Scalar> dF_drho(num_atoms, 0.0);

    for (std::size_t i = 0; i < num_atoms; ++i) {
        SplineResult F_result = embedding_.eval(rho_sum[i]);
        results.energy += F_result.value;
        dF_drho[i] = F_result.derivative;
    }

    // Third pass: compute pair energies and forces
    for (std::size_t i = 0; i < num_atoms; ++i) {
        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq < cutoff_sq && r_sq > 1e-10) {
                Scalar r = std::sqrt(r_sq);

                // Pair potential phi(r) = Z(r)^2 / r
                SplineResult phi_result = pair_potential(r);
                Scalar phi = phi_result.value;
                Scalar dphi = phi_result.derivative;

                // Half energy (full neighbor list)
                results.energy += 0.5 * phi;

                if (compute_forces || compute_virial) {
                    // Electron density derivative
                    Scalar drho = rho_.eval(r).derivative;

                    // Total force on atom i from interaction with atom j
                    // The total potential energy is:
                    //   E = sum_i F(rho_i) + 0.5 * sum_i sum_{j!=i} phi(r_ij)
                    // Force on atom i:
                    //   F_i = -dE/dr_i = -dF/drho_i * drho_i/dr_i - 0.5 * sum_j dphi/dr_ij * dr_ij/dr_i
                    //                   - sum_j dF/drho_j * drho_j/dr_i
                    // where drho_i/dr_i = sum_j drho(r_ij)/dr * (-r_hat_ij)
                    //       drho_j/dr_i = drho(r_ij)/dr * r_hat_ij (contribution from r_ij to rho_j)
                    // Full neighbor list: visit each pair twice, so halve pair contributions
                    // For force contributions, we add to atom i when visiting (i,j), and this
                    // is half the total force from pair ij plus the density contributions
                    Scalar force_over_r = -(dphi + (dF_drho[i] + dF_drho[j]) * drho) / r;

                    Vec3 force = force_over_r * dr;

                    if (compute_forces) {
                        // Full neighbor list: add to atom i only
                        // When visiting (j,i), the opposite force is added to atom j
                        system.forces().col(i) -= force.array();
                    }

                    if (compute_virial) {
                        // Halve virial because pairs are counted twice
                        results.virial -= 0.5 * dr * force.transpose();
                    }
                }
            }
        }
    }

    return results;
}

// TabulatedAlloyEAM implementation

inline void TabulatedAlloyEAM::load(const std::string& filename) {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Cannot open EAM alloy file: " + filename);
    }

    std::string line;

    // Lines 1-3: Comments
    for (int i = 0; i < 3; ++i) {
        std::getline(file, line);
    }

    // Line 4: nel, element_1, element_2, ...
    std::getline(file, line);
    {
        std::istringstream iss(line);
        iss >> num_elements_;
        element_info_.resize(num_elements_);
        for (int i = 0; i < num_elements_; ++i) {
            iss >> element_info_[i].symbol;
            symbol_to_index_[element_info_[i].symbol] = i;
        }
    }

    // Line 5: nF, dF, nr, dr, cutoff
    int nF, nr;
    Scalar dF, dr;
    std::getline(file, line);
    {
        std::istringstream iss(line);
        iss >> nF >> dF >> nr >> dr >> cutoff_;
    }

    // Allocate arrays
    embedding_.resize(num_elements_);
    rho_.resize(num_elements_);
    phi_.resize(num_elements_, std::vector<CubicSpline>(num_elements_));

    // Read per-element data
    for (int elem = 0; elem < num_elements_; ++elem) {
        // Element info line: Z, mass, a0, lattice
        std::getline(file, line);
        {
            std::istringstream iss(line);
            iss >> element_info_[elem].atomic_number
                >> element_info_[elem].mass
                >> element_info_[elem].lattice_constant
                >> element_info_[elem].lattice_type;
        }
        Z_to_index_[element_info_[elem].atomic_number] = elem;

        // Read F(rho) values
        std::vector<Scalar> F_values(nF);
        for (int i = 0; i < nF; ++i) {
            file >> F_values[i];
        }
        embedding_[elem].init(0.0, static_cast<Scalar>(nF - 1) * dF, F_values);

        // Read rho(r) values
        std::vector<Scalar> rho_values(nr);
        for (int i = 0; i < nr; ++i) {
            file >> rho_values[i];
        }
        rho_[elem].init(0.0, static_cast<Scalar>(nr - 1) * dr, rho_values);
    }

    // Read pair potentials phi_{ij}(r) for i >= j
    // Note: In setfl format, phi values are stored as r*phi(r)
    for (int i = 0; i < num_elements_; ++i) {
        for (int j = 0; j <= i; ++j) {
            std::vector<Scalar> phi_values(nr);
            for (int k = 0; k < nr; ++k) {
                file >> phi_values[k];
                // Convert from r*phi(r) to phi(r)
                Scalar r = static_cast<Scalar>(k) * dr;
                if (r > 1e-10) {
                    phi_values[k] /= r;
                } else {
                    phi_values[k] = 0.0;
                }
            }
            phi_[i][j].init(0.0, static_cast<Scalar>(nr - 1) * dr, phi_values);

            // Symmetric: phi_[j][i] = phi_[i][j]
            if (i != j) {
                phi_[j][i] = phi_[i][j];
            }
        }
    }
}

inline int TabulatedAlloyEAM::element_index(const std::string& symbol) const {
    auto it = symbol_to_index_.find(symbol);
    if (it != symbol_to_index_.end()) {
        return it->second;
    }
    return -1;
}

inline int TabulatedAlloyEAM::element_index_by_Z(int Z) const {
    auto it = Z_to_index_.find(Z);
    if (it != Z_to_index_.end()) {
        return it->second;
    }
    return -1;
}

inline std::vector<std::string> TabulatedAlloyEAM::element_symbols() const {
    std::vector<std::string> symbols;
    symbols.reserve(num_elements_);
    for (const auto& info : element_info_) {
        symbols.push_back(info.symbol);
    }
    return symbols;
}

inline SplineResult TabulatedAlloyEAM::pair_potential(int elem_i, int elem_j, Scalar r) const {
    return phi_.at(elem_i).at(elem_j).eval(r);
}

inline PotentialResults TabulatedAlloyEAM::compute_impl(
    AtomicSystem& system,
    NeighborList& neighbors,
    bool compute_forces,
    bool compute_virial)
{
    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();
    const Mat3& cell = system.cell();
    const Scalar cutoff_sq = cutoff_ * cutoff_;

    // Map atomic numbers to element indices
    std::vector<int> type_map(num_atoms, -1);
    for (std::size_t i = 0; i < num_atoms; ++i) {
        int Z = system.atomic_numbers()(i);
        type_map[i] = element_index_by_Z(Z);
    }

    // First pass: compute electron densities at each atom
    std::vector<Scalar> rho_sum(num_atoms, 0.0);

    for (std::size_t i = 0; i < num_atoms; ++i) {
        int type_i = type_map[i];
        if (type_i < 0) continue;

        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;
            int type_j = type_map[j];
            if (type_j < 0) continue;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq < cutoff_sq && r_sq > 1e-10) {
                Scalar r = std::sqrt(r_sq);
                // Density from atom j contributes to atom i
                rho_sum[i] += rho_[type_j].eval(r).value;
            }
        }
    }

    // Second pass: compute embedding energies and store derivatives
    std::vector<Scalar> dF_drho(num_atoms, 0.0);

    for (std::size_t i = 0; i < num_atoms; ++i) {
        int type_i = type_map[i];
        if (type_i < 0) continue;

        SplineResult F_result = embedding_[type_i].eval(rho_sum[i]);
        results.energy += F_result.value;
        dF_drho[i] = F_result.derivative;
    }

    // Third pass: compute pair energies and forces
    for (std::size_t i = 0; i < num_atoms; ++i) {
        int type_i = type_map[i];
        if (type_i < 0) continue;

        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;
            int type_j = type_map[j];
            if (type_j < 0) continue;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq < cutoff_sq && r_sq > 1e-10) {
                Scalar r = std::sqrt(r_sq);

                // Pair potential for this element pair
                SplineResult phi_result = phi_[type_i][type_j].eval(r);
                Scalar phi = phi_result.value;
                Scalar dphi = phi_result.derivative;

                // Half energy (full neighbor list)
                results.energy += 0.5 * phi;

                if (compute_forces || compute_virial) {
                    // Electron density derivatives
                    Scalar drho_j = rho_[type_j].eval(r).derivative;
                    Scalar drho_i = rho_[type_i].eval(r).derivative;

                    // Total force on atom i from interaction with atom j
                    // Same logic as TabulatedEAM above
                    Scalar force_over_r = -(dphi + dF_drho[i] * drho_j + dF_drho[j] * drho_i) / r;

                    Vec3 force = force_over_r * dr;

                    if (compute_forces) {
                        // Full neighbor list: add to atom i only
                        system.forces().col(i) -= force.array();
                    }

                    if (compute_virial) {
                        // Halve virial because pairs are counted twice
                        results.virial -= 0.5 * dr * force.transpose();
                    }
                }
            }
        }
    }

    return results;
}

} // namespace atomistica
