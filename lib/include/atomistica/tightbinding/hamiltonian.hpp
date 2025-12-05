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

#include <map>
#include <vector>

#include "../config.hpp"
#include "../core/atomic_system.hpp"
#include "../core/neighbor_list.hpp"
#include "materials.hpp"
#include "slater_koster.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Compute distance vector between atoms i and j using neighbor info
 *
 * @param system Atomic system
 * @param i Central atom index
 * @param neighbor Neighbor information (index and cell_shift)
 * @return Distance vector r_j - r_i accounting for periodic images
 */
inline Vec3 neighbor_distance_vector(const AtomicSystem& system, std::size_t i,
                                     const Neighbor& neighbor) {
    Vec3 r_i = system.position(i);
    Vec3 r_j = system.position(neighbor.index);
    Vec3 shift = system.cell() * Vec3(neighbor.cell_shift[0],
                                       neighbor.cell_shift[1],
                                       neighbor.cell_shift[2]);
    return r_j + shift - r_i;
}

/**
 * @brief Tight-binding Hamiltonian builder
 *
 * Constructs H and S matrices from atomic positions using
 * Slater-Koster transformations and tabulated integrals.
 */
class TBHamiltonian {
public:
    TBHamiltonian() = default;

    /**
     * @brief Set materials database
     */
    void set_materials(MaterialsDatabase* db) { materials_ = db; }

    /**
     * @brief Initialize Hamiltonian for a given atomic system
     *
     * Sets up orbital indices and allocates matrices
     *
     * @param system Atomic system
     * @param elements Vector of element parameters for each atom type
     */
    void init(const AtomicSystem& system, const std::vector<TBElementParams>& elements);

    /**
     * @brief Build H and S matrices
     *
     * @param system Atomic system
     * @param neighbors Neighbor list
     */
    void build_matrices(const AtomicSystem& system, const NeighborList& neighbors);

    /**
     * @brief Add SCC (self-consistent charge) correction to Hamiltonian
     *
     * Modifies H based on Mulliken charges: H_scc = H + shift * S
     * where shift depends on charge difference from neutral
     *
     * @param gamma Gamma matrix (Coulomb interaction between atoms)
     */
    void add_scc_correction(const MatX& gamma);

    /**
     * @brief Compute repulsive energy
     */
    Scalar compute_repulsive_energy(const AtomicSystem& system,
                                    const NeighborList& neighbors);

    /**
     * @brief Compute repulsive forces
     */
    void compute_repulsive_forces(const AtomicSystem& system,
                                  const NeighborList& neighbors,
                                  MatX3& forces);

    /**
     * @brief Get Hamiltonian data structure
     */
    DenseHamiltonian& hamiltonian() { return ham_; }
    const DenseHamiltonian& hamiltonian() const { return ham_; }

    /**
     * @brief Get H matrix
     */
    MatX& H() { return ham_.H; }
    const MatX& H() const { return ham_.H; }

    /**
     * @brief Get S matrix
     */
    MatX& S() { return ham_.S; }
    const MatX& S() const { return ham_.S; }

private:
    MaterialsDatabase* materials_ = nullptr;
    DenseHamiltonian ham_;
    std::vector<TBElementParams> element_params_;
    std::map<int, int> z_to_elem_;
};

/**
 * @brief Compute gamma matrix for SCC-DFTB
 *
 * gamma_ij = short-range Coulomb interaction between atoms i and j
 * Uses Hubbard U parameters and distance-dependent function.
 *
 * @param system Atomic system
 * @param elements Element parameters
 * @param use_periodic Include periodic images
 * @return Gamma matrix
 */
MatX compute_gamma_matrix(const AtomicSystem& system,
                          const std::vector<TBElementParams>& elements,
                          const std::vector<int>& elem_index,
                          bool use_periodic = false);

} // namespace tb
} // namespace atomistica
