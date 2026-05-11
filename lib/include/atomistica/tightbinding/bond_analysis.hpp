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

#include <vector>
#include <cmath>
#include <stdexcept>

#include <Eigen/Eigenvalues>

#include "../config.hpp"
#include "../core/neighbor_list.hpp"
#include "types.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Result of bond analysis for a single bond
 */
struct BondProperties {
    int atom_i;                    // First atom index
    int atom_j;                    // Second atom index
    std::array<int, 3> cell_shift; // Periodic image shift

    Scalar overlap_population;     // Mulliken overlap population: Tr(rho_ij * S_ji)
    Scalar loewdin_bond_order;     // Loewdin bond order: Tr(P_ij + P_ji) / 2
    Scalar covalent_energy;        // Covalent bond energy (Bornsen et al.)
};

/**
 * @brief Compute matrix square root using eigendecomposition
 *
 * Computes sqrt(S) where S is symmetric positive definite.
 * Uses S = V * D * V^T, so sqrt(S) = V * sqrt(D) * V^T
 *
 * @param S Symmetric positive definite matrix
 * @return Matrix square root
 */
inline MatX matrix_sqrt(const MatX& S) {
    Eigen::SelfAdjointEigenSolver<MatX> solver(S);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed in matrix_sqrt");
    }

    VecX eigenvalues = solver.eigenvalues();
    MatX eigenvectors = solver.eigenvectors();

    // Check for positive definiteness
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues[i] < -1e-10) {
            throw std::runtime_error("Matrix is not positive definite in matrix_sqrt");
        }
        eigenvalues[i] = std::sqrt(std::max(eigenvalues[i], 0.0));
    }

    // Reconstruct: sqrt(S) = V * sqrt(D) * V^T
    return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

/**
 * @brief Compute matrix inverse square root using eigendecomposition
 *
 * Computes S^(-1/2) where S is symmetric positive definite.
 *
 * @param S Symmetric positive definite matrix
 * @return Matrix inverse square root
 */
inline MatX matrix_inv_sqrt(const MatX& S) {
    Eigen::SelfAdjointEigenSolver<MatX> solver(S);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed in matrix_inv_sqrt");
    }

    VecX eigenvalues = solver.eigenvalues();
    MatX eigenvectors = solver.eigenvectors();

    // Compute inverse square root of eigenvalues
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues[i] < 1e-10) {
            throw std::runtime_error("Matrix is singular in matrix_inv_sqrt");
        }
        eigenvalues[i] = 1.0 / std::sqrt(eigenvalues[i]);
    }

    // Reconstruct: S^(-1/2) = V * D^(-1/2) * V^T
    return eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose();
}

/**
 * @brief Bond analysis for tight-binding calculations
 *
 * Provides various bond analysis methods:
 * - Mulliken overlap population
 * - Loewdin bond order
 * - Covalent bond energy (Bornsen et al., J.Phys.: Cond. Mat. 11, L287 (1999))
 */
class BondAnalyzer {
public:
    BondAnalyzer() = default;

    /**
     * @brief Perform full bond analysis
     *
     * Computes overlap population, Loewdin bond order, and covalent energy
     * for all bonds in the neighbor list.
     *
     * @param ham Hamiltonian with computed density matrix
     * @param neighbors Neighbor list
     * @param compute_loewdin If true, compute Loewdin bond order (requires matrix sqrt)
     * @return Vector of bond properties for each bond in neighbor list
     */
    std::vector<BondProperties> analyze(const DenseHamiltonian& ham,
                                        const NeighborList& neighbors,
                                        bool compute_loewdin = true);

    /**
     * @brief Compute only Mulliken overlap populations
     *
     * Faster than full analysis when only overlap population is needed.
     */
    std::vector<Scalar> compute_overlap_populations(const DenseHamiltonian& ham,
                                                    const NeighborList& neighbors);

    /**
     * @brief Compute Loewdin charges
     *
     * Loewdin charges are based on the orthogonalized density matrix:
     * q_i^Loewdin = sum_a P_aa where P = S^(1/2) * rho * S^(1/2)
     *
     * @param ham Hamiltonian with computed density matrix
     * @return Loewdin charges for each atom
     */
    VecX compute_loewdin_charges(const DenseHamiltonian& ham);

    /**
     * @brief Compute orbital-resolved density of states contribution
     *
     * Returns the diagonal of rho * S which gives the local orbital occupations.
     *
     * @param ham Hamiltonian with computed density matrix
     * @return Diagonal of rho * S (orbital occupations)
     */
    VecX compute_orbital_occupations(const DenseHamiltonian& ham);

    /**
     * @brief Compute total bond order between two atoms
     *
     * Sum of Loewdin bond orders over all bonds between atoms i and j
     * (including periodic images).
     *
     * @param bonds Vector of bond properties from analyze()
     * @param atom_i First atom index
     * @param atom_j Second atom index
     * @return Total Loewdin bond order
     */
    static Scalar total_bond_order(const std::vector<BondProperties>& bonds,
                                   int atom_i, int atom_j);

    /**
     * @brief Compute total covalent energy for an atom
     *
     * Sum of half the covalent bond energies for all bonds involving atom i.
     *
     * @param bonds Vector of bond properties from analyze()
     * @param atom_i Atom index
     * @return Total covalent energy contribution for atom i
     */
    static Scalar atom_covalent_energy(const std::vector<BondProperties>& bonds,
                                       int atom_i);
};

} // namespace tb
} // namespace atomistica
