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
#include <string>
#include <vector>
#include <cstdint>

#include "../config.hpp"

namespace atomistica {
namespace tb {

/**
 * @brief Orbital angular momentum types
 */
enum class OrbitalType : int {
    S = 0,    // l = 0
    Px = 1,   // l = 1, m = +1
    Py = 2,   // l = 1, m = -1
    Pz = 3,   // l = 1, m = 0
    Dxy = 4,  // l = 2
    Dyz = 5,
    Dzx = 6,
    Dx2y2 = 7,
    D3z2r2 = 8
};

/**
 * @brief Maximum number of orbitals per atom (s + 3p + 5d)
 */
constexpr int MAX_ORBITALS = 9;

/**
 * @brief Slater-Koster integral types
 */
enum class SKIntegralType : int {
    dds = 0,  // d-d sigma
    ddp = 1,  // d-d pi
    ddd = 2,  // d-d delta
    pds = 3,  // p-d sigma
    pdp = 4,  // p-d pi
    pps = 5,  // p-p sigma
    ppp = 6,  // p-p pi
    sds = 7,  // s-d sigma
    sps = 8,  // s-p sigma
    sss = 9   // s-s sigma
};

/**
 * @brief Number of independent Slater-Koster integrals
 */
constexpr int NUM_SK_INTEGRALS = 10;

/**
 * @brief Angular momentum for each orbital type
 */
inline constexpr std::array<int, MAX_ORBITALS> ORBITAL_L = {0, 1, 1, 1, 2, 2, 2, 2, 2};

/**
 * @brief Element parameters for tight-binding
 */
struct TBElementParams {
    std::string symbol;           // Element symbol (e.g., "C")
    int atomic_number = 0;        // Z
    int num_orbitals = 0;         // Number of valence orbitals (1, 4, or 9)
    int l_max = 0;                // Maximum angular momentum (0=s, 1=sp, 2=spd)

    std::array<int, MAX_ORBITALS> l;           // Angular momentum per orbital
    std::array<Scalar, MAX_ORBITALS> onsite;   // On-site energies (eV)

    Scalar hubbard_U = 0.0;       // Hubbard U parameter (for SCC)
    Scalar valence_electrons = 0.0; // Number of valence electrons (q0)

    // Spin-polarization parameters (optional)
    bool has_spin = false;
    std::array<std::array<Scalar, 3>, 3> W_spin = {}; // W parameters for spin

    TBElementParams() {
        l.fill(-1);
        onsite.fill(0.0);
    }

    /**
     * @brief Check if element has s orbitals only
     */
    bool is_s_only() const { return num_orbitals == 1; }

    /**
     * @brief Check if element has sp orbitals
     */
    bool is_sp() const { return num_orbitals == 4; }

    /**
     * @brief Check if element has spd orbitals
     */
    bool is_spd() const { return num_orbitals == 9; }
};

/**
 * @brief Pair parameters for tight-binding (H, S, V_rep tables)
 */
struct TBPairParams {
    int Z1 = 0;
    int Z2 = 0;
    Scalar cutoff = 0.0;          // Cutoff for H/S tables
    Scalar cutoff_rep = 0.0;      // Cutoff for repulsive potential

    // Tabulated Slater-Koster integrals H(r) and S(r)
    // Each array contains spline data for interpolation
    std::vector<Scalar> r_grid;   // Distance grid points
    std::vector<std::array<Scalar, NUM_SK_INTEGRALS>> H_table;  // H integrals
    std::vector<std::array<Scalar, NUM_SK_INTEGRALS>> S_table;  // S integrals

    // Repulsive potential V_rep(r)
    std::vector<Scalar> r_rep_grid;
    std::vector<Scalar> V_rep;

    bool is_valid() const { return Z1 > 0 && Z2 > 0 && !r_grid.empty(); }
};

/**
 * @brief Hamiltonian storage for tight-binding calculations
 */
struct DenseHamiltonian {
    int num_atoms = 0;
    int num_orbitals = 0;         // Total orbitals across all atoms

    // Matrix storage (column-major for LAPACK compatibility)
    MatX H;                       // Hamiltonian [norb x norb]
    MatX S;                       // Overlap [norb x norb]
    MatX rho;                     // Density matrix [norb x norb]
    MatX e_matrix;                // H * rho (for force calculations)

    // Eigenvalue results
    VecX eigenvalues;             // [norb]
    MatX eigenvectors;            // [norb x norb]
    VecX occupation;              // Occupation numbers [norb]

    // Per-atom data
    std::vector<int> orbitals_per_atom;    // Number of orbitals on each atom
    std::vector<int> orbital_offset;        // First orbital index for each atom
    std::vector<int> element_index;         // Element type for each atom

    // Mulliken charges
    VecX charges;                 // Net charges (q0 - q)
    VecX neutral_charges;         // Reference neutral charges (q0)

    // Energies
    Scalar band_energy = 0.0;     // E_bs
    Scalar repulsive_energy = 0.0; // E_rep
    Scalar fermi_level = 0.0;     // Chemical potential mu

    void resize(int nat, int norb) {
        num_atoms = nat;
        num_orbitals = norb;

        H = MatX::Zero(norb, norb);
        S = MatX::Zero(norb, norb);
        rho = MatX::Zero(norb, norb);
        e_matrix = MatX::Zero(norb, norb);

        eigenvalues = VecX::Zero(norb);
        eigenvectors = MatX::Zero(norb, norb);
        occupation = VecX::Zero(norb);

        orbitals_per_atom.resize(nat, 0);
        orbital_offset.resize(nat, 0);
        element_index.resize(nat, -1);

        charges = VecX::Zero(nat);
        neutral_charges = VecX::Zero(nat);
    }

    void clear_matrices() {
        H.setZero();
        S.setZero();
        rho.setZero();
        e_matrix.setZero();
    }
};

/**
 * @brief SCC (Self-Consistent Charges) parameters
 */
struct SCCParams {
    int max_iterations = 200;
    Scalar convergence_threshold = 1e-4;  // Max charge change
    Scalar mixing_parameter = 0.2;        // Beta for simple mixing
    int anderson_memory = 3;              // History length for Anderson mixing

    // DFTB3 extensions
    bool enable_dftb3 = false;
    Scalar zeta = 0.0;                    // DFTB3 damping parameter
};

/**
 * @brief Solver parameters
 */
struct SolverParams {
    Scalar electronic_temperature = 0.01; // kT in eV (default ~116 K)
    bool use_divide_and_conquer = true;   // Use dsygvd vs dsygv
};

} // namespace tb
} // namespace atomistica
