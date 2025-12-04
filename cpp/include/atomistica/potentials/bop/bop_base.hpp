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
#include <functional>
#include <map>
#include <tuple>
#include <vector>

#include "../../config.hpp"
#include "../../math/cutoff_functions.hpp"
#include "../potential_base.hpp"
#include "screening.hpp"

namespace atomistica {

/**
 * @brief Maximum number of element types for BOP potentials
 */
constexpr int BOP_MAX_ELEMENTS = 10;

/**
 * @brief Compute pair index for symmetric pair (i,j) with i <= j
 */
inline int pair_index(int i, int j, int n_elements) {
    if (i > j) std::swap(i, j);
    return i * n_elements - (i * (i - 1)) / 2 + (j - i);
}

/**
 * @brief Number of unique pairs for n elements
 */
inline int num_pairs(int n_elements) {
    return n_elements * (n_elements + 1) / 2;
}

/**
 * @brief Base parameters common to all BOP potentials
 *
 * These are parameters that depend on pair type (element i, element j)
 */
struct BOPPairParams {
    // Pair potential parameters
    Scalar A = 0.0;       // Repulsive amplitude
    Scalar lambda = 0.0;  // Repulsive decay
    Scalar B = 0.0;       // Attractive amplitude
    Scalar mu = 0.0;      // Attractive decay

    // Cutoff parameters
    Scalar r1 = 0.0;      // Inner cutoff
    Scalar r2 = 0.0;      // Outer cutoff

    // Precomputed cutoff
    TrigOffCutoff cutoff;

    void init_cutoff() {
        cutoff.init(r1, r2);
    }
};

/**
 * @brief Angular parameters for BOP potentials
 *
 * These are parameters that depend on triplet type (i-j-k)
 */
struct BOPAngularParams {
    Scalar gamma = 1.0;   // Angular function amplitude
    Scalar c = 0.0;       // Angular function numerator
    Scalar d = 1.0;       // Angular function denominator
    Scalar h = 0.0;       // Angular function cos offset
    Scalar c2 = 0.0;      // Precomputed c*c
    Scalar d2 = 0.0;      // Precomputed d*d
    Scalar c2_d2 = 0.0;   // Precomputed c*c/d*d

    void precompute() {
        c2 = c * c;
        d2 = d * d;
        c2_d2 = c2 / d2;
    }
};

/**
 * @brief Element-specific bond-order parameters
 */
struct BOPElementParams {
    Scalar beta = 1.0;    // Bond-order parameter
    Scalar n = 1.0;       // Bond-order exponent
    Scalar xi = 1.0;      // Bond-order scaling
    Scalar omega = 1.0;   // Angular modulation factor

    // Precomputed values
    Scalar half_n = 0.5;
    Scalar minus_half_over_n = -0.5;

    void precompute() {
        half_n = 0.5 * n;
        if (std::abs(n) > 1e-10) {
            minus_half_over_n = -0.5 / n;
        } else {
            minus_half_over_n = 0.0;
        }
    }
};

/**
 * @brief Internal bond data computed during neighbor list traversal
 */
struct BondData {
    std::size_t j;              // Neighbor index
    int pair_type;              // Pair type index
    Scalar r;                   // Distance
    Vec3 dr;                    // Distance vector (rj - ri)
    Vec3 unit;                  // Unit vector
    Scalar fc;                  // Cutoff function value (pair/attractive)
    Scalar dfc;                 // Cutoff function derivative
    Scalar fc_bo;               // Bond-order cutoff (may differ with screening)
    Scalar dfc_bo;              // Bond-order cutoff derivative
    std::array<int, 3> shift;   // Periodic shift

    // Screening data (only used when Screening=true)
    Scalar S;                   // Screening factor
    Scalar dS_drij;             // dS/dr_ij
    std::vector<ScreeningNeighbor> screening_neighbors;  // Atoms contributing to screening
};

/**
 * @brief CRTP base class for Bond-Order Potentials
 *
 * Implements the common BOP algorithm:
 *   E = 0.5 * sum_{i,j} fc(r_ij) * [V_R(r_ij) + b_ij * V_A(r_ij)]
 *
 * where b_ij is the bond order computed from angular terms:
 *   b_ij = f(z_ij)
 *   z_ij = sum_k fc(r_ik) * g(cos(theta_jik)) * h(r_ik, r_ij)
 *
 * @tparam Derived The derived potential class (CRTP)
 * @tparam Screening Whether screening is enabled (compile-time toggle)
 */
template<typename Derived, bool Screening = false>
class BOPBase : public PotentialBase<BOPBase<Derived, Screening>> {
public:
    using Base = PotentialBase<BOPBase<Derived, Screening>>;
    friend Base;

    BOPBase() = default;

    /**
     * @brief Get maximum cutoff radius
     */
    Scalar cutoff() const {
        return derived().cutoff_impl();
    }

    /**
     * @brief Compute energy, forces, and virial
     */
    PotentialResults compute(AtomicSystem& system,
                            NeighborList& neighbors,
                            bool compute_forces = true,
                            bool compute_virial = true) {
        return compute_impl(system, neighbors, compute_forces, compute_virial);
    }

protected:
    /**
     * @brief Main BOP computation kernel
     */
    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial) {
        PotentialResults results;
        const std::size_t num_atoms = system.num_atoms();
        const Mat3& cell = system.cell();

        // Thread-local bond storage
        std::vector<BondData> bonds;
        bonds.reserve(50);  // Typical coordination

        for (std::size_t i = 0; i < num_atoms; ++i) {
            int Zi = system.atomic_numbers()(i);
            int eli = derived().element_index(Zi);
            if (eli < 0) continue;

            Vec3 ri = system.position(i).matrix();

            // Build bond list for atom i
            bonds.clear();
            auto [nb_begin, nb_end] = neighbors.neighbors(i);

            // First pass: collect all potential bonds with basic data
            for (auto it = nb_begin; it != nb_end; ++it) {
                const auto& neigh = *it;
                std::size_t j = neigh.index;
                int Zj = system.atomic_numbers()(j);
                int elj = derived().element_index(Zj);
                if (elj < 0) continue;

                int ptype = derived().pair_type(eli, elj);

                // Compute distance
                Vec3 rj = system.position(j).matrix();
                Vec3 dr = rj - ri;
                dr += cell.col(0) * neigh.cell_shift[0];
                dr += cell.col(1) * neigh.cell_shift[1];
                dr += cell.col(2) * neigh.cell_shift[2];

                Scalar r = dr.norm();
                Scalar cutoff_r = derived().pair_cutoff(ptype);

                if (r >= cutoff_r || r < 1e-10) continue;

                // Evaluate cutoff function
                auto [fc, dfc] = derived().cutoff_function(ptype, r);
                if (fc < 1e-15) continue;

                BondData bond;
                bond.j = j;
                bond.pair_type = ptype;
                bond.r = r;
                bond.dr = dr;
                bond.unit = dr / r;
                bond.fc = fc;
                bond.dfc = dfc;
                bond.fc_bo = fc;    // Default: same as pair cutoff
                bond.dfc_bo = dfc;
                bond.shift = neigh.cell_shift;
                bond.S = 1.0;
                bond.dS_drij = 0.0;

                bonds.push_back(bond);
            }

            // Second pass (screening only): compute screening factors
            if constexpr (Screening) {
                for (std::size_t b_ij = 0; b_ij < bonds.size(); ++b_ij) {
                    auto& bond_ij = bonds[b_ij];

                    // Get screening parameters for this pair type
                    const auto& scr_params = derived().screening_params(bond_ij.pair_type);

                    // Collect r_ik vectors for all potential screening atoms
                    std::vector<Vec3> rik_vectors;
                    rik_vectors.reserve(bonds.size());
                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik != b_ij) {
                            rik_vectors.push_back(bonds[b_ik].dr);
                        }
                    }

                    // Compute screening (simple version without neighbor tracking for now)
                    auto scr_result = compute_screening_simple(
                        scr_params, bond_ij.dr, bond_ij.r * bond_ij.r, rik_vectors);

                    bond_ij.S = scr_result.fully_screened ? 0.0 : std::exp(scr_result.S);
                    bond_ij.dS_drij = scr_result.dS_drij * bond_ij.S;

                    // Apply screening to cutoff functions
                    // Following Fortran: cutfcnbo = (1-fCin)*S*fCbo + fCin
                    // For simplicity, we use: fc_bo = S * fc
                    bond_ij.fc_bo = bond_ij.S * bond_ij.fc;
                    bond_ij.dfc_bo = bond_ij.S * bond_ij.dfc + bond_ij.dS_drij * bond_ij.fc / bond_ij.r;

                    // Also apply to pair cutoff for energy
                    bond_ij.fc = bond_ij.fc_bo;
                    bond_ij.dfc = bond_ij.dfc_bo;

                    // Skip fully screened bonds
                    if (bond_ij.S < 1e-10) {
                        bond_ij.fc = 0.0;
                        bond_ij.fc_bo = 0.0;
                    }
                }
            }

            // Compute pair energies and bond orders
            for (std::size_t b_ij = 0; b_ij < bonds.size(); ++b_ij) {
                const auto& bond_ij = bonds[b_ij];
                std::size_t j = bond_ij.j;
                int elj = derived().element_index(system.atomic_numbers()(j));

                // Pair potentials
                auto [VR, dVR] = derived().repulsive(bond_ij.pair_type, bond_ij.r);
                auto [VA, dVA] = derived().attractive(bond_ij.pair_type, bond_ij.r);

                // Compute bond order z_ij = sum_k fc_ik * g(cos_jik) * h(...)
                // Note: Use fc_bo (bond-order cutoff) which may be screened
                Scalar zij = 0.0;
                std::vector<Scalar> dz_dcos(bonds.size(), 0.0);
                std::vector<Scalar> dz_drik(bonds.size(), 0.0);
                std::vector<Scalar> dz_drij_via_h(bonds.size(), 0.0);

                for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                    if (b_ik == b_ij) continue;

                    const auto& bond_ik = bonds[b_ik];
                    int elk = derived().element_index(system.atomic_numbers()(bond_ik.j));

                    // Cosine of angle j-i-k
                    Scalar cos_jik = bond_ij.unit.dot(bond_ik.unit);

                    // Angular function g(cos_jik)
                    auto [g_val, dg] = derived().angular_function(
                        eli, elj, elk, bond_ij.pair_type, bond_ik.pair_type, cos_jik);

                    // Distance-dependent function h
                    auto [h_val, dh_drik, dh_drij] = derived().distance_function(
                        eli, elj, elk, bond_ij.pair_type, bond_ik.pair_type,
                        bond_ij.r, bond_ik.r);

                    // Contribution to z_ij (using bond-order cutoff fc_bo)
                    Scalar contrib = bond_ik.fc_bo * g_val * h_val;
                    zij += contrib;

                    // Store derivatives for force calculation (using dfc_bo)
                    dz_dcos[b_ik] = bond_ik.fc_bo * dg * h_val;
                    dz_drik[b_ik] = bond_ik.dfc_bo * g_val * h_val + bond_ik.fc_bo * g_val * dh_drik;
                    dz_drij_via_h[b_ik] = bond_ik.fc_bo * g_val * dh_drij;
                }

                // Bond order function b(z)
                auto [bij, dbij] = derived().bond_order(eli, bond_ij.pair_type, zij);

                // Total pair energy (factor 0.5 for half contribution)
                Scalar E_pair = 0.5 * bond_ij.fc * (VR + bij * VA);
                results.energy += E_pair;

                if (compute_forces || compute_virial) {
                    // Following Fortran BOP kernel structure:
                    // E = 0.5 * fc * (VR + b(z) * VA)
                    // F = -dE/dr
                    //
                    // Pair contribution (without bond-order derivative):
                    // dE/dr_ij = 0.5 * [dfc/dr * (VR + b*VA) + fc * (dVR/dr + b*dVA/dr)]
                    //
                    // Bond-order contribution:
                    // dE/dz = 0.5 * fc * db/dz * VA
                    // F_x = -dE/dz * dz/dr_x for x = i, j, k

                    // Prefactor for bond-order derivative term
                    // Note: dbij_dzij in Fortran = db/dz * VA * fc (scaled by 0.5 later)
                    Scalar dbij_dzij = 0.5 * bond_ij.fc * dbij * VA;

                    // Accumulate dz_ij/dr_i, dz_ij/dr_j, dz_ij/dr_k
                    // Following Fortran: dbidi, dbidj, dbidk
                    Vec3 dbidi = Vec3::Zero();  // dz/dr_i
                    Vec3 dbidj = Vec3::Zero();  // dz/dr_j
                    std::vector<Vec3> dbidk(bonds.size(), Vec3::Zero());  // dz/dr_k for each k

                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik == b_ij) continue;

                        const auto& bond_ik = bonds[b_ik];
                        Scalar cos_jik = bond_ij.unit.dot(bond_ik.unit);

                        // Angular coordinate derivatives for cos(theta_jik)
                        // cos = unit_ij · unit_ik
                        //
                        // Using the fact that cos = (r_ij · r_ik) / (r_ij * r_ik):
                        // d(cos)/d(r_i) = -(unit_ik - cos*unit_ij)/r_ij - (unit_ij - cos*unit_ik)/r_ik
                        // d(cos)/d(r_j) = (unit_ik - cos*unit_ij)/r_ij
                        // d(cos)/d(r_k) = (unit_ij - cos*unit_ik)/r_ik
                        //
                        // Note: These are the pure geometric derivatives assuming j and k are
                        // independent. The Fortran dcsdjk term accounts for when j and k have
                        // a fixed relationship, which is not the case in our full neighbor list.

                        Vec3 term_ij = (bond_ik.unit - cos_jik * bond_ij.unit) / bond_ij.r;
                        Vec3 term_ik = (bond_ij.unit - cos_jik * bond_ik.unit) / bond_ik.r;

                        Vec3 dcsdi = -term_ij - term_ik;
                        Vec3 dcsdj = term_ij;
                        Vec3 dcsdk = term_ik;

                        // dzfac = fc_ik * dg/dcos * h (Fortran line 1260)
                        Scalar dzfac = dz_dcos[b_ik];  // Already contains fc_ik * dg * h

                        // Angular contributions to dz/dr (Fortran lines 1262-1264):
                        // dgdi = dzfac * dcsdi, etc.
                        Vec3 dgdi = dzfac * dcsdi;
                        Vec3 dgdj = dzfac * dcsdj;
                        Vec3 dgdk = dzfac * dcsdk;

                        // Radial contributions from h and fc_ik
                        // dzdrij = g * fc_ik * dh/dr_ij (from dz_drij_via_h)
                        // dzdrik = g * (dfc_ik/dr_ik * h + fc_ik * dh/dr_ik) (from dz_drik)
                        Scalar dzdrij = dz_drij_via_h[b_ik];
                        Scalar dzdrik = dz_drik[b_ik];

                        // Accumulate (Fortran lines 1305, 1311-1312, 1319):
                        // dbidi = dbidi - dzdrij*unit_ij - dzdrik*unit_ik + dgdi
                        // dbidj = dbidj + dzdrij*unit_ij + dgdj
                        // dbidk = dzdrik*unit_ik + dgdk
                        dbidi += -dzdrij * bond_ij.unit - dzdrik * bond_ik.unit + dgdi;
                        dbidj += dzdrij * bond_ij.unit + dgdj;
                        dbidk[b_ik] = dzdrik * bond_ik.unit + dgdk;
                    }

                    // Pair radial force (without bond-order term)
                    // dffac = 0.5 * (dVR/dr * fc + b*dVA/dr * fc + VR * dfc/dr + b*VA * dfc/dr)
                    Scalar dffac = 0.5 * (dVR * bond_ij.fc + bij * dVA * bond_ij.fc +
                                          VR * bond_ij.dfc + bij * VA * bond_ij.dfc);

                    // df = dffac * unit_ij (Fortran line 1400)
                    // fi = fi + df (Fortran line 1401)
                    // fj = fj - df (Fortran line 1402)
                    Vec3 df_pair = dffac * bond_ij.unit;
                    Vec3 force_on_i = df_pair;
                    Vec3 force_on_j = -df_pair;

                    // Bond-order forces (Fortran lines 1417-1426):
                    // fi = fi - dbij_dzij * dbidi
                    // fj = fj - dbij_dzij * dbidj
                    force_on_i -= dbij_dzij * dbidi;
                    force_on_j -= dbij_dzij * dbidj;

                    // Forces on neighbors k (Fortran lines 1432-1467)
                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik == b_ij) continue;

                        const auto& bond_ik = bonds[b_ik];
                        std::size_t k = bond_ik.j;

                        // fk = fk - dbij_dzij * dbidk (Fortran line 1464)
                        Vec3 force_on_k = -dbij_dzij * dbidk[b_ik];

                        if (compute_forces) {
                            system.forces().col(k) += force_on_k.array();
                        }

                        if (compute_virial) {
                            // Virial: -outer_product(r_ik, force_on_k)
                            results.virial -= bond_ik.dr * force_on_k.transpose();
                        }
                    }

                    if (compute_forces) {
                        system.forces().col(i) += force_on_i.array();
                        system.forces().col(j) += force_on_j.array();
                    }

                    if (compute_virial) {
                        // Virial from pair force and bond-order on j
                        // Fortran: wij = wij + outer_product(rij, df) - dbij_dzij*wijb
                        // For simplicity, just add the r_ij * f_j contribution
                        results.virial -= bond_ij.dr * force_on_j.transpose();
                    }
                }
            }
        }

        return results;
    }

private:
    Derived& derived() { return static_cast<Derived&>(*this); }
    const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

} // namespace atomistica
