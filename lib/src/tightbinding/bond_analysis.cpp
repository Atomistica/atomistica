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

#include <atomistica/tightbinding/bond_analysis.hpp>

namespace atomistica {
namespace tb {

std::vector<BondProperties> BondAnalyzer::analyze(const DenseHamiltonian& ham,
                                                   const NeighborList& neighbors,
                                                   bool compute_loewdin) {
    std::vector<BondProperties> results;

    int nat = ham.num_atoms;

    // Compute Loewdin orthogonalized density if requested
    MatX loewdin_rho;
    if (compute_loewdin) {
        MatX sqrt_S = matrix_sqrt(ham.S);
        loewdin_rho = sqrt_S * ham.rho * sqrt_S;
    }

    // Iterate over all bonds
    for (std::size_t i = 0; i < static_cast<std::size_t>(nat); ++i) {
        int offset_i = ham.orbital_offset[i];
        int norb_i = ham.orbitals_per_atom[i];

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;

            // Only process each bond once (i < j)
            // But keep all bonds for full neighbor information
            if (j < i) continue;

            int offset_j = ham.orbital_offset[j];
            int norb_j = ham.orbitals_per_atom[j];

            BondProperties bond;
            bond.atom_i = static_cast<int>(i);
            bond.atom_j = static_cast<int>(j);
            bond.cell_shift = it->cell_shift;

            // Compute overlap population: sum over orbital pairs
            // overlap_population = sum_{a in i, b in j} rho(a,b) * S(b,a)
            Scalar overlap_pop = 0.0;
            Scalar loewdin_bo = 0.0;
            Scalar e_cov = 0.0;

            for (int a = 0; a < norb_i; ++a) {
                int ia = offset_i + a;

                for (int b = 0; b < norb_j; ++b) {
                    int jb = offset_j + b;

                    // Mulliken overlap population
                    overlap_pop += ham.rho(ia, jb) * ham.S(jb, ia);

                    // Loewdin bond order
                    if (compute_loewdin) {
                        loewdin_bo += loewdin_rho(ia, jb) + loewdin_rho(jb, ia);
                    }

                    // Covalent bond energy (Bornsen et al.)
                    // E_cov = rho(a,b) * (H(b,a) - 0.5*S(b,a)*(H(a,a) + H(b,b)))
                    // Note: For SCC-DFTB, the electrostatic shift cancels out
                    Scalar H_ab = ham.H(jb, ia);
                    Scalar S_ab = ham.S(jb, ia);
                    Scalar H_aa = ham.H(ia, ia);
                    Scalar H_bb = ham.H(jb, jb);

                    e_cov += ham.rho(ia, jb) * (H_ab - 0.5 * S_ab * (H_aa + H_bb));
                }
            }

            bond.overlap_population = overlap_pop;
            bond.loewdin_bond_order = 0.5 * loewdin_bo;  // Factor of 0.5 from definition
            bond.covalent_energy = e_cov;

            results.push_back(bond);
        }
    }

    return results;
}

std::vector<Scalar> BondAnalyzer::compute_overlap_populations(const DenseHamiltonian& ham,
                                                               const NeighborList& neighbors) {
    std::vector<Scalar> results;

    int nat = ham.num_atoms;

    for (std::size_t i = 0; i < static_cast<std::size_t>(nat); ++i) {
        int offset_i = ham.orbital_offset[i];
        int norb_i = ham.orbitals_per_atom[i];

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            std::size_t j = it->index;
            if (j < i) continue;

            int offset_j = ham.orbital_offset[j];
            int norb_j = ham.orbitals_per_atom[j];

            Scalar overlap_pop = 0.0;

            for (int a = 0; a < norb_i; ++a) {
                int ia = offset_i + a;
                for (int b = 0; b < norb_j; ++b) {
                    int jb = offset_j + b;
                    overlap_pop += ham.rho(ia, jb) * ham.S(jb, ia);
                }
            }

            results.push_back(overlap_pop);
        }
    }

    return results;
}

VecX BondAnalyzer::compute_loewdin_charges(const DenseHamiltonian& ham) {
    int nat = ham.num_atoms;

    // Compute Loewdin density
    MatX sqrt_S = matrix_sqrt(ham.S);
    MatX loewdin_rho = sqrt_S * ham.rho * sqrt_S;

    VecX charges = VecX::Zero(nat);

    for (int i = 0; i < nat; ++i) {
        int offset = ham.orbital_offset[i];
        int norb = ham.orbitals_per_atom[i];

        Scalar q = 0.0;
        for (int a = 0; a < norb; ++a) {
            q += loewdin_rho(offset + a, offset + a);
        }

        // Net charge = neutral - actual
        charges[i] = ham.neutral_charges[i] - q;
    }

    return charges;
}

VecX BondAnalyzer::compute_orbital_occupations(const DenseHamiltonian& ham) {
    MatX rhoS = ham.rho * ham.S;
    return rhoS.diagonal();
}

Scalar BondAnalyzer::total_bond_order(const std::vector<BondProperties>& bonds,
                                       int atom_i, int atom_j) {
    Scalar total = 0.0;

    for (const auto& bond : bonds) {
        if ((bond.atom_i == atom_i && bond.atom_j == atom_j) ||
            (bond.atom_i == atom_j && bond.atom_j == atom_i)) {
            total += bond.loewdin_bond_order;
        }
    }

    return total;
}

Scalar BondAnalyzer::atom_covalent_energy(const std::vector<BondProperties>& bonds,
                                           int atom_i) {
    Scalar total = 0.0;

    for (const auto& bond : bonds) {
        if (bond.atom_i == atom_i || bond.atom_j == atom_i) {
            // Each bond is counted once, but shared between two atoms
            total += 0.5 * bond.covalent_energy;
        }
    }

    return total;
}

} // namespace tb
} // namespace atomistica
