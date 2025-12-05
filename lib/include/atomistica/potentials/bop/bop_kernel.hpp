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
 */
struct BOPPairParams {
    Scalar A = 0.0;       // Repulsive amplitude
    Scalar lambda = 0.0;  // Repulsive decay
    Scalar B = 0.0;       // Attractive amplitude
    Scalar mu = 0.0;      // Attractive decay

    Scalar r1 = 0.0;      // Inner cutoff
    Scalar r2 = 0.0;      // Outer cutoff

    TrigOffCutoff cutoff;

    void init_cutoff() {
        cutoff.init(r1, r2);
    }
};

/**
 * @brief Element-specific bond-order parameters
 */
struct BOPElementParams {
    Scalar beta = 1.0;
    Scalar n = 1.0;
    Scalar xi = 1.0;
    Scalar omega = 1.0;

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
 * @brief Internal bond data for unscreened BOP
 */
struct BondData {
    std::size_t j;              // Neighbor index
    int pair_type;              // Pair type index
    Scalar r;                   // Distance
    Vec3 dr;                    // Distance vector (rj - ri)
    Vec3 unit;                  // Unit vector
    Scalar fc;                  // Cutoff function value
    Scalar dfc;                 // Cutoff function derivative
    std::array<int, 3> shift;   // Periodic shift
};

/**
 * @brief Internal bond data for screened BOP (includes screening info)
 */
struct ScreenedBondData : public BondData {
    Scalar S;                   // Screening factor (0 to 1)
    Scalar dS_drij;             // dS/dr_ij

    // Per-neighbor screening derivatives
    // screening_k[m] corresponds to bonds[m] that screens this bond
    std::vector<std::size_t> screening_k_idx;  // Index into bonds array
    std::vector<Scalar> dS_drik;               // dS/dr_ik for each screening atom
    std::vector<Scalar> dS_drjk;               // dS/dr_jk for each screening atom
};

/**
 * @brief CRTP base class for unscreened Bond-Order Potentials
 */
template<typename Derived>
class BOPKernel : public PotentialBase<BOPKernel<Derived>> {
public:
    using Base = PotentialBase<BOPKernel<Derived>>;
    friend Base;

    BOPKernel() = default;

    Scalar cutoff() const {
        return derived().cutoff_impl();
    }

    PotentialResults compute(AtomicSystem& system,
                            NeighborList& neighbors,
                            bool compute_forces = true,
                            bool compute_virial = true) {
        return compute_unscreened(system, neighbors, compute_forces, compute_virial);
    }

protected:
    PotentialResults compute_unscreened(AtomicSystem& system,
                                        NeighborList& neighbors,
                                        bool compute_forces,
                                        bool compute_virial) {
        PotentialResults results;
        const std::size_t num_atoms = system.num_atoms();
        const Mat3& cell = system.cell();

        std::vector<BondData> bonds;
        bonds.reserve(50);

        for (std::size_t i = 0; i < num_atoms; ++i) {
            int Zi = system.atomic_numbers()(i);
            int eli = derived().element_index(Zi);
            if (eli < 0) continue;

            Vec3 ri = system.position(i);

            // Build bond list for atom i
            bonds.clear();
            auto [nb_begin, nb_end] = neighbors.neighbors(i);

            for (auto it = nb_begin; it != nb_end; ++it) {
                const auto& neigh = *it;
                std::size_t j = neigh.index;
                int Zj = system.atomic_numbers()(j);
                int elj = derived().element_index(Zj);
                if (elj < 0) continue;

                int ptype = derived().pair_type(eli, elj);

                Vec3 rj = system.position(j);
                Vec3 dr = rj - ri;
                dr += cell.col(0) * neigh.cell_shift[0];
                dr += cell.col(1) * neigh.cell_shift[1];
                dr += cell.col(2) * neigh.cell_shift[2];

                Scalar r = dr.norm();
                Scalar cutoff_r = derived().pair_cutoff(ptype);

                if (r >= cutoff_r || r < 1e-10) continue;

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
                bond.shift = neigh.cell_shift;

                bonds.push_back(bond);
            }

            // Compute energies and forces for each bond
            for (std::size_t b_ij = 0; b_ij < bonds.size(); ++b_ij) {
                const auto& bond_ij = bonds[b_ij];
                std::size_t j = bond_ij.j;
                int elj = derived().element_index(system.atomic_numbers()(j));

                auto [VR, dVR] = derived().repulsive(bond_ij.pair_type, bond_ij.r);
                auto [VA, dVA] = derived().attractive(bond_ij.pair_type, bond_ij.r);

                // Compute bond order z_ij
                Scalar zij = 0.0;
                std::vector<Scalar> dz_dcos(bonds.size(), 0.0);
                std::vector<Scalar> dz_drik(bonds.size(), 0.0);
                std::vector<Scalar> dz_drij_via_h(bonds.size(), 0.0);

                for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                    if (b_ik == b_ij) continue;

                    const auto& bond_ik = bonds[b_ik];
                    int elk = derived().element_index(system.atomic_numbers()(bond_ik.j));

                    Scalar cos_jik = bond_ij.unit.dot(bond_ik.unit);

                    auto [g_val, dg] = derived().angular_function(
                        eli, elj, elk, bond_ij.pair_type, bond_ik.pair_type, cos_jik);

                    auto [h_val, dh_drik, dh_drij] = derived().distance_function(
                        eli, elj, elk, bond_ij.pair_type, bond_ik.pair_type,
                        bond_ij.r, bond_ik.r);

                    Scalar contrib = bond_ik.fc * g_val * h_val;
                    zij += contrib;

                    dz_dcos[b_ik] = bond_ik.fc * dg * h_val;
                    dz_drik[b_ik] = bond_ik.dfc * g_val * h_val + bond_ik.fc * g_val * dh_drik;
                    dz_drij_via_h[b_ik] = bond_ik.fc * g_val * dh_drij;
                }

                auto [bij, dbij] = derived().bond_order(eli, bond_ij.pair_type, zij);

                Scalar E_pair = 0.5 * bond_ij.fc * (VR + bij * VA);
                results.energy += E_pair;

                if (compute_forces || compute_virial) {
                    Scalar dbij_dzij = 0.5 * bond_ij.fc * dbij * VA;

                    Vec3 dbidi = Vec3::Zero();
                    Vec3 dbidj = Vec3::Zero();
                    std::vector<Vec3> dbidk(bonds.size(), Vec3::Zero());

                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik == b_ij) continue;

                        const auto& bond_ik = bonds[b_ik];
                        Scalar cos_jik = bond_ij.unit.dot(bond_ik.unit);

                        Vec3 term_ij = (bond_ik.unit - cos_jik * bond_ij.unit) / bond_ij.r;
                        Vec3 term_ik = (bond_ij.unit - cos_jik * bond_ik.unit) / bond_ik.r;

                        Vec3 dcsdi = -term_ij - term_ik;
                        Vec3 dcsdj = term_ij;
                        Vec3 dcsdk = term_ik;

                        Scalar dzfac = dz_dcos[b_ik];

                        Vec3 dgdi = dzfac * dcsdi;
                        Vec3 dgdj = dzfac * dcsdj;
                        Vec3 dgdk = dzfac * dcsdk;

                        Scalar dzdrij = dz_drij_via_h[b_ik];
                        Scalar dzdrik = dz_drik[b_ik];

                        dbidi += -dzdrij * bond_ij.unit - dzdrik * bond_ik.unit + dgdi;
                        dbidj += dzdrij * bond_ij.unit + dgdj;
                        dbidk[b_ik] = dzdrik * bond_ik.unit + dgdk;
                    }

                    Scalar dffac = 0.5 * (dVR * bond_ij.fc + bij * dVA * bond_ij.fc +
                                          VR * bond_ij.dfc + bij * VA * bond_ij.dfc);

                    Vec3 df_pair = dffac * bond_ij.unit;
                    Vec3 force_on_i = df_pair;
                    Vec3 force_on_j = -df_pair;

                    force_on_i -= dbij_dzij * dbidi;
                    force_on_j -= dbij_dzij * dbidj;

                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik == b_ij) continue;

                        const auto& bond_ik = bonds[b_ik];
                        std::size_t k = bond_ik.j;

                        Vec3 force_on_k = -dbij_dzij * dbidk[b_ik];

                        if (compute_forces) {
                            system.forces().col(k) += force_on_k.array();
                        }

                        if (compute_virial) {
                            results.virial -= bond_ik.dr * force_on_k.transpose();
                        }
                    }

                    if (compute_forces) {
                        system.forces().col(i) += force_on_i.array();
                        system.forces().col(j) += force_on_j.array();
                    }

                    if (compute_virial) {
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

/**
 * @brief CRTP base class for screened Bond-Order Potentials
 *
 * This class implements the full screening algorithm including
 * force derivatives from screening atom movements.
 */
template<typename Derived>
class ScreenedBOPKernel : public PotentialBase<ScreenedBOPKernel<Derived>> {
public:
    using Base = PotentialBase<ScreenedBOPKernel<Derived>>;
    friend Base;

    ScreenedBOPKernel() = default;

    Scalar cutoff() const {
        return derived().cutoff_impl();
    }

    PotentialResults compute(AtomicSystem& system,
                            NeighborList& neighbors,
                            bool compute_forces = true,
                            bool compute_virial = true) {
        return compute_screened(system, neighbors, compute_forces, compute_virial);
    }

protected:
    PotentialResults compute_screened(AtomicSystem& system,
                                      NeighborList& neighbors,
                                      bool compute_forces,
                                      bool compute_virial) {
        PotentialResults results;
        const std::size_t num_atoms = system.num_atoms();
        const Mat3& cell = system.cell();

        std::vector<ScreenedBondData> bonds;
        bonds.reserve(50);

        for (std::size_t i = 0; i < num_atoms; ++i) {
            int Zi = system.atomic_numbers()(i);
            int eli = derived().element_index(Zi);
            if (eli < 0) continue;

            Vec3 ri = system.position(i);

            // Build bond list for atom i
            bonds.clear();
            auto [nb_begin, nb_end] = neighbors.neighbors(i);

            // First pass: collect all bonds with basic data
            for (auto it = nb_begin; it != nb_end; ++it) {
                const auto& neigh = *it;
                std::size_t j = neigh.index;
                int Zj = system.atomic_numbers()(j);
                int elj = derived().element_index(Zj);
                if (elj < 0) continue;

                int ptype = derived().pair_type(eli, elj);

                Vec3 rj = system.position(j);
                Vec3 dr = rj - ri;
                dr += cell.col(0) * neigh.cell_shift[0];
                dr += cell.col(1) * neigh.cell_shift[1];
                dr += cell.col(2) * neigh.cell_shift[2];

                Scalar r = dr.norm();
                Scalar cutoff_r = derived().screened_cutoff(ptype);

                if (r >= cutoff_r || r < 1e-10) continue;

                auto [fc, dfc] = derived().cutoff_function(ptype, r);

                ScreenedBondData bond;
                bond.j = j;
                bond.pair_type = ptype;
                bond.r = r;
                bond.dr = dr;
                bond.unit = dr / r;
                bond.fc = fc;
                bond.dfc = dfc;
                bond.shift = neigh.cell_shift;
                bond.S = 1.0;
                bond.dS_drij = 0.0;

                bonds.push_back(bond);
            }

            // Second pass: compute screening for each bond
            for (std::size_t b_ij = 0; b_ij < bonds.size(); ++b_ij) {
                auto& bond_ij = bonds[b_ij];
                const auto& scr_params = derived().screening_params(bond_ij.pair_type);

                Scalar rij_sq = bond_ij.r * bond_ij.r;
                const Scalar C_dr_cut_rij_sq = scr_params.C_dr_cut * rij_sq;

                bond_ij.S = 0.0;  // Log-space accumulation
                bond_ij.dS_drij = 0.0;
                bond_ij.screening_k_idx.clear();
                bond_ij.dS_drik.clear();
                bond_ij.dS_drjk.clear();

                bool fully_screened = false;

                // Check all other bonds for screening contribution
                for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                    if (b_ik == b_ij) continue;

                    const auto& bond_ik = bonds[b_ik];
                    Scalar rik_sq = bond_ik.r * bond_ik.r;

                    // Skip atoms too far away to screen
                    if (rik_sq >= C_dr_cut_rij_sq) continue;

                    // r_jk = r_ik - r_ij
                    Vec3 rjk_vec = bond_ik.dr - bond_ij.dr;
                    Scalar rjk_sq = rjk_vec.squaredNorm();

                    // Geometric check: k must be between i and j
                    Scalar dot_ij_ik = bond_ij.dr.dot(bond_ik.dr);
                    Scalar dot_ij_jk = bond_ij.dr.dot(rjk_vec);

                    if (dot_ij_ik <= scr_params.dot_threshold ||
                        dot_ij_jk >= -scr_params.dot_threshold) {
                        continue;
                    }

                    // Compute screening parameter C
                    Scalar xik = rik_sq / rij_sq;
                    Scalar xjk = rjk_sq / rij_sq;

                    Scalar xik_m_xjk = xik - xjk;
                    Scalar xik_p_xjk = xik + xjk;

                    Scalar denom = 1.0 - xik_m_xjk * xik_m_xjk;
                    if (std::abs(denom) < 1e-15) continue;

                    Scalar fac = 1.0 / denom;
                    Scalar C = (2.0 * xik_p_xjk - xik_m_xjk * xik_m_xjk - 1.0) * fac;

                    if (C <= scr_params.Cmin) {
                        fully_screened = true;
                        break;
                    }

                    if (C < scr_params.Cmax) {
                        // Partial screening contribution
                        Scalar Cmax_C = scr_params.Cmax - C;
                        Scalar C_Cmin = C - scr_params.Cmin;
                        Scalar ratio = Cmax_C / C_Cmin;

                        bond_ij.S -= ratio * ratio;

                        // Derivatives of C
                        Scalar dCdxik = 4.0 * xik * fac * (1.0 + (C - 1.0) * xik_m_xjk);
                        Scalar dCdxjk = 4.0 * xjk * fac * (1.0 - (C - 1.0) * xik_m_xjk);

                        // dS/dC
                        Scalar dSdC = 2.0 * Cmax_C * scr_params.dC / (C_Cmin * C_Cmin * C_Cmin);

                        // Accumulate dS/d(rij^2)
                        Scalar dCdrij_sq = -(dCdxik * xik + dCdxjk * xjk) / rij_sq;
                        bond_ij.dS_drij += dSdC * dCdrij_sq;

                        // Store per-neighbor derivatives
                        bond_ij.screening_k_idx.push_back(b_ik);
                        // dS/d(rik) = dS/dC * dC/dxik * d(xik)/d(rik^2) * 2*rik
                        //           = dSdC * dCdxik / rij_sq * 2
                        bond_ij.dS_drik.push_back(dSdC * dCdxik * 2.0 / rij_sq);
                        // dS/d(rjk) via xjk
                        bond_ij.dS_drjk.push_back(dSdC * dCdxjk * 2.0 / rij_sq);
                    }
                }

                if (fully_screened || bond_ij.S < scr_params.screening_threshold) {
                    bond_ij.S = 0.0;
                    bond_ij.dS_drij = 0.0;
                    bond_ij.screening_k_idx.clear();
                    bond_ij.dS_drik.clear();
                    bond_ij.dS_drjk.clear();
                } else {
                    // Convert from log-space: S_actual = exp(S_log)
                    Scalar S_exp = std::exp(bond_ij.S);
                    // dS_actual/dr = S_actual * dS_log/dr
                    bond_ij.dS_drij *= S_exp;
                    for (auto& ds : bond_ij.dS_drik) ds *= S_exp;
                    for (auto& ds : bond_ij.dS_drjk) ds *= S_exp;
                    bond_ij.S = S_exp;

                    // Convert dS_drij from d/d(rij^2) to d/d(rij)
                    bond_ij.dS_drij *= 2.0 * bond_ij.r;
                }

                // Apply screening to cutoff
                bond_ij.fc *= bond_ij.S;
                bond_ij.dfc = bond_ij.S * bond_ij.dfc + bond_ij.dS_drij * bond_ij.fc / bond_ij.r;
            }

            // Compute energies and forces
            for (std::size_t b_ij = 0; b_ij < bonds.size(); ++b_ij) {
                const auto& bond_ij = bonds[b_ij];
                if (bond_ij.S < 1e-15) continue;  // Skip fully screened bonds

                std::size_t j = bond_ij.j;
                int elj = derived().element_index(system.atomic_numbers()(j));

                auto [VR, dVR] = derived().repulsive(bond_ij.pair_type, bond_ij.r);
                auto [VA, dVA] = derived().attractive(bond_ij.pair_type, bond_ij.r);

                // Compute bond order z_ij (using screened cutoff)
                Scalar zij = 0.0;
                std::vector<Scalar> dz_dcos(bonds.size(), 0.0);
                std::vector<Scalar> dz_drik(bonds.size(), 0.0);
                std::vector<Scalar> dz_drij_via_h(bonds.size(), 0.0);

                for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                    if (b_ik == b_ij) continue;

                    const auto& bond_ik = bonds[b_ik];
                    if (bond_ik.S < 1e-15) continue;

                    int elk = derived().element_index(system.atomic_numbers()(bond_ik.j));

                    Scalar cos_jik = bond_ij.unit.dot(bond_ik.unit);

                    auto [g_val, dg] = derived().angular_function(
                        eli, elj, elk, bond_ij.pair_type, bond_ik.pair_type, cos_jik);

                    auto [h_val, dh_drik, dh_drij] = derived().distance_function(
                        eli, elj, elk, bond_ij.pair_type, bond_ik.pair_type,
                        bond_ij.r, bond_ik.r);

                    // Use screened cutoff for bond order
                    Scalar fc_ik = bond_ik.fc;
                    Scalar dfc_ik = bond_ik.dfc;

                    Scalar contrib = fc_ik * g_val * h_val;
                    zij += contrib;

                    dz_dcos[b_ik] = fc_ik * dg * h_val;
                    dz_drik[b_ik] = dfc_ik * g_val * h_val + fc_ik * g_val * dh_drik;
                    dz_drij_via_h[b_ik] = fc_ik * g_val * dh_drij;
                }

                auto [bij, dbij] = derived().bond_order(eli, bond_ij.pair_type, zij);

                // Use screened cutoff for energy
                Scalar fc_ij = bond_ij.fc;
                Scalar E_pair = 0.5 * fc_ij * (VR + bij * VA);
                results.energy += E_pair;

                if (compute_forces || compute_virial) {
                    Scalar dfc_ij = bond_ij.dfc;
                    Scalar dbij_dzij = 0.5 * fc_ij * dbij * VA;

                    Vec3 dbidi = Vec3::Zero();
                    Vec3 dbidj = Vec3::Zero();
                    std::vector<Vec3> dbidk(bonds.size(), Vec3::Zero());

                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik == b_ij) continue;

                        const auto& bond_ik = bonds[b_ik];
                        if (bond_ik.S < 1e-15) continue;

                        Scalar cos_jik = bond_ij.unit.dot(bond_ik.unit);

                        Vec3 term_ij = (bond_ik.unit - cos_jik * bond_ij.unit) / bond_ij.r;
                        Vec3 term_ik = (bond_ij.unit - cos_jik * bond_ik.unit) / bond_ik.r;

                        Vec3 dcsdi = -term_ij - term_ik;
                        Vec3 dcsdj = term_ij;
                        Vec3 dcsdk = term_ik;

                        Scalar dzfac = dz_dcos[b_ik];

                        Vec3 dgdi = dzfac * dcsdi;
                        Vec3 dgdj = dzfac * dcsdj;
                        Vec3 dgdk = dzfac * dcsdk;

                        Scalar dzdrij = dz_drij_via_h[b_ik];
                        Scalar dzdrik = dz_drik[b_ik];

                        dbidi += -dzdrij * bond_ij.unit - dzdrik * bond_ik.unit + dgdi;
                        dbidj += dzdrij * bond_ij.unit + dgdj;
                        dbidk[b_ik] = dzdrik * bond_ik.unit + dgdk;
                    }

                    // Pair radial force (including screened cutoff derivative)
                    Scalar dffac = 0.5 * (dVR * fc_ij + bij * dVA * fc_ij +
                                          VR * dfc_ij + bij * VA * dfc_ij);

                    Vec3 df_pair = dffac * bond_ij.unit;
                    Vec3 force_on_i = df_pair;
                    Vec3 force_on_j = -df_pair;

                    force_on_i -= dbij_dzij * dbidi;
                    force_on_j -= dbij_dzij * dbidj;

                    // Forces on bond-order neighbors k
                    for (std::size_t b_ik = 0; b_ik < bonds.size(); ++b_ik) {
                        if (b_ik == b_ij) continue;

                        const auto& bond_ik = bonds[b_ik];
                        if (bond_ik.S < 1e-15) continue;

                        std::size_t k = bond_ik.j;

                        Vec3 force_on_k = -dbij_dzij * dbidk[b_ik];

                        if (compute_forces) {
                            system.forces().col(k) += force_on_k.array();
                        }

                        if (compute_virial) {
                            results.virial -= bond_ik.dr * force_on_k.transpose();
                        }
                    }

                    // Screening forces: dE/dS * dS/dr_k
                    // E_pair = 0.5 * S * fc_base * (VR + b*VA)
                    // dE/dS = 0.5 * fc_base * (VR + b*VA) = E_pair / S
                    if (bond_ij.S > 1e-10) {
                        Scalar dE_dS = E_pair / bond_ij.S;

                        // Forces on screening atoms
                        for (std::size_t s = 0; s < bond_ij.screening_k_idx.size(); ++s) {
                            std::size_t b_ik = bond_ij.screening_k_idx[s];
                            const auto& bond_ik = bonds[b_ik];
                            std::size_t k = bond_ik.j;

                            // r_jk = r_ik - r_ij
                            Vec3 rjk_vec = bond_ik.dr - bond_ij.dr;
                            Scalar rjk = rjk_vec.norm();
                            Vec3 rjk_unit = (rjk > 1e-10) ? Vec3(rjk_vec / rjk) : Vec3::Zero();

                            // dS/dr_ik already scaled by 2/rij_sq
                            // Actual derivative: dS/d(r_ik) = dS_drik * r_ik
                            Scalar dS_drik = bond_ij.dS_drik[s] * bond_ik.r;
                            Scalar dS_drjk = bond_ij.dS_drjk[s] * rjk;

                            // Force from r_ik contribution
                            // F_i += dE/dS * dS/dr_ik * (-unit_ik)
                            // F_k += dE/dS * dS/dr_ik * (+unit_ik)
                            Vec3 f_ik = dE_dS * dS_drik * bond_ik.unit;
                            force_on_i -= f_ik;
                            Vec3 force_k_from_ik = f_ik;

                            // Force from r_jk contribution
                            // r_jk = r_k - r_j, so dr_jk/dr_k = +1, dr_jk/dr_j = -1
                            // F_j += dE/dS * dS/dr_jk * (-unit_jk)
                            // F_k += dE/dS * dS/dr_jk * (+unit_jk)
                            Vec3 f_jk = dE_dS * dS_drjk * rjk_unit;
                            force_on_j -= f_jk;
                            Vec3 force_k_from_jk = f_jk;

                            Vec3 total_force_k = force_k_from_ik + force_k_from_jk;

                            if (compute_forces) {
                                system.forces().col(k) += total_force_k.array();
                            }

                            if (compute_virial) {
                                // Virial from screening forces
                                results.virial -= bond_ik.dr * force_k_from_ik.transpose();
                                // For jk contribution, use r_jk vector (from j to k)
                                results.virial -= rjk_vec * force_k_from_jk.transpose();
                            }
                        }
                    }

                    if (compute_forces) {
                        system.forces().col(i) += force_on_i.array();
                        system.forces().col(j) += force_on_j.array();
                    }

                    if (compute_virial) {
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
