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
#include <complex>
#include <stdexcept>
#include <vector>
#include <array>
#include <algorithm>

#include "../../config.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../potential_base.hpp"
#include "coulomb.hpp"

namespace atomistica {

/**
 * @brief Fast Multipole Method (FMM) for long-range Coulomb interactions
 *
 * Implements a balanced-tree fast multipole method based on spherical
 * harmonic expansions. This method computes electrostatic interactions
 * with O(N) complexity.
 *
 * The algorithm consists of:
 * 1. Sort particles into an octree
 * 2. M2M: Compute multipoles at leaves, propagate upward
 * 3. M2L: Convert distant multipoles to local expansions
 * 4. L2L: Propagate local expansions downward
 * 5. Near-field: Direct computation for nearby particles
 *
 * References:
 * - Greengard & Rokhlin, J. Comput. Phys. 73, 325 (1987)
 * - Lambert et al., J. Comput. Phys. 126, 274 (1996) [periodicity]
 *
 * @note This is a simplified serial implementation without MPI parallelization.
 */
class FMMCoulomb : public PotentialBase<FMMCoulomb> {
public:
    /**
     * @brief Construct FMM Coulomb solver
     * @param l_max Maximum angular momentum for multipole expansion (default 8)
     * @param n_level Number of tree levels (default 3)
     * @param leaf_size Maximum particles per leaf (default 200)
     * @param periodic_images Periodicity parameter k for summing 3^k images (default 1)
     */
    explicit FMMCoulomb(int l_max = 8, int n_level = 3,
                        int leaf_size = 200, int periodic_images = 1);

    ~FMMCoulomb();

    // Disable copy
    FMMCoulomb(const FMMCoulomb&) = delete;
    FMMCoulomb& operator=(const FMMCoulomb&) = delete;

    // Move operations
    FMMCoulomb(FMMCoulomb&& other) noexcept;
    FMMCoulomb& operator=(FMMCoulomb&& other) noexcept;

    // Parameters
    void set_l_max(int l_max);
    void set_n_level(int n_level);
    void set_leaf_size(int leaf_size);
    void set_periodic_images(int k);

    int l_max() const { return l_max_; }
    int n_level() const { return n_level_; }
    int leaf_size() const { return leaf_size_; }

    void set_charges(const std::vector<Scalar>& charges) {
        charges_ = charges;
    }

    const std::vector<Scalar>& charges() const { return charges_; }

    // CRTP implementation
    Scalar cutoff_impl() const {
        // FMM doesn't use a cutoff - it handles all interactions
        return 0.0;
    }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    // Tree level structure
    struct TreeLevel {
        std::array<int, 3> n;           // Number of cells in each direction
        int n_tot;                       // Total cells on this level
        std::vector<int> np;             // Number of particles per cell
        std::vector<std::vector<Scalar>> Ml0;    // Multipole moments (real, l=0..l_max)
        std::vector<std::vector<std::complex<Scalar>>> Mlm;  // Multipole moments (complex)
        std::vector<std::vector<Scalar>> Ll0;    // Local expansion (real)
        std::vector<std::vector<std::complex<Scalar>>> Llm;  // Local expansion (complex)
    };

    // Solid harmonics
    void compute_solid_harmonic_R(Scalar r, Scalar costh, Scalar phi,
                                   std::vector<Scalar>& Rl0,
                                   std::vector<std::complex<Scalar>>& Rlm) const;

    void compute_solid_harmonic_I(Scalar r, Scalar costh, Scalar phi,
                                   int max_l,
                                   std::vector<Scalar>& Il0,
                                   std::vector<std::complex<Scalar>>& Ilm) const;

    // Coordinate transforms
    void cartesian_to_spherical(const Vec3& r, Scalar& rr, Scalar& costh, Scalar& phi) const;

    // Multipole operations
    void multipole_to_multipole(const Vec3& dr, int l_max,
                                const std::vector<Scalar>& Ml0_child,
                                const std::vector<std::complex<Scalar>>& Mlm_child,
                                std::vector<Scalar>& Ml0,
                                std::vector<std::complex<Scalar>>& Mlm) const;

    void multipole_to_local(const Vec3& dr, int l_max,
                           const std::vector<Scalar>& Ml0,
                           const std::vector<std::complex<Scalar>>& Mlm,
                           std::vector<Scalar>& Ll0,
                           std::vector<std::complex<Scalar>>& Llm) const;

    void local_to_local(const Vec3& dr, int l_max_in,
                       const std::vector<Scalar>& Ll0_in,
                       const std::vector<std::complex<Scalar>>& Llm_in,
                       int l_max_out,
                       std::vector<Scalar>& Ll0_out,
                       std::vector<std::complex<Scalar>>& Llm_out) const;

    // Index functions
    int lm_index(int l, int m) const {
        // Maps (l,m) with m>0 to linear index
        // l=1: m=1 -> index 0
        // l=2: m=1 -> index 1, m=2 -> index 2
        // l=l: m=m -> index l*(l-1)/2 + m - 1
        return l * (l - 1) / 2 + m - 1;
    }

    int n_off_diag() const {
        return l_max_ * (l_max_ + 1) / 2;
    }

    // Cell indexing
    int cell_index(int x, int y, int z, const std::array<int, 3>& n) const {
        return z + n[2] * (y + n[1] * x);
    }

    std::array<int, 3> cell_coords(int idx, const std::array<int, 3>& n) const {
        int x = idx / (n[1] * n[2]);
        idx -= x * n[1] * n[2];
        int y = idx / n[2];
        int z = idx % n[2];
        return {x, y, z};
    }

    Vec3 cell_center(int x, int y, int z, int level) const;

    // Algorithm steps
    void initialize_tree(AtomicSystem& system);
    void sort_particles(AtomicSystem& system);
    void M2M_pass(AtomicSystem& system);
    void M2L_pass(AtomicSystem& system);
    void L2L_pass();
    Scalar near_field_and_local(AtomicSystem& system,
                                bool compute_forces,
                                bool compute_virial,
                                Mat3& virial);

    // Parameters
    int l_max_ = 8;
    int n_level_ = 3;
    int leaf_size_ = 200;
    int k_ = 1;  // periodic images

    std::vector<Scalar> charges_;

    // Tree structure
    Vec3 lower_, upper_;
    std::array<int, 3> n0_ = {2, 2, 2};  // Root decomposition
    std::vector<TreeLevel> tree_;

    // Particle -> cell mapping
    std::vector<int> np_;                          // Particles per leaf
    std::vector<std::vector<int>> p_in_leaf_;      // Particle indices in each leaf

    // Root multipole moments and local expansion
    std::vector<Scalar> root_Ml0_;
    std::vector<std::complex<Scalar>> root_Mlm_;
    std::vector<Scalar> root_Ll0_;
    std::vector<std::complex<Scalar>> root_Llm_;

    // Neighbor offsets (26 neighbors excluding self)
    static constexpr int NUM_NEIGHBORS = 26;
    std::array<std::array<int, 3>, NUM_NEIGHBORS> neighbor_offsets_;

    bool initialized_ = false;
};

// ============================================================================
// Implementation
// ============================================================================

inline FMMCoulomb::FMMCoulomb(int l_max, int n_level, int leaf_size, int periodic_images)
    : l_max_(l_max)
    , n_level_(n_level)
    , leaf_size_(leaf_size)
    , k_(periodic_images) {

    // Initialize neighbor offsets
    int idx = 0;
    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                if (dx != 0 || dy != 0 || dz != 0) {
                    neighbor_offsets_[idx++] = {dx, dy, dz};
                }
            }
        }
    }
}

inline FMMCoulomb::~FMMCoulomb() = default;
inline FMMCoulomb::FMMCoulomb(FMMCoulomb&& other) noexcept = default;
inline FMMCoulomb& FMMCoulomb::operator=(FMMCoulomb&& other) noexcept = default;

inline void FMMCoulomb::set_l_max(int l_max) {
    if (l_max < 1) throw std::invalid_argument("l_max must be >= 1");
    l_max_ = l_max;
    initialized_ = false;
}

inline void FMMCoulomb::set_n_level(int n_level) {
    if (n_level < 1) throw std::invalid_argument("n_level must be >= 1");
    n_level_ = n_level;
    initialized_ = false;
}

inline void FMMCoulomb::set_leaf_size(int leaf_size) {
    if (leaf_size < 1) throw std::invalid_argument("leaf_size must be >= 1");
    leaf_size_ = leaf_size;
    initialized_ = false;
}

inline void FMMCoulomb::set_periodic_images(int k) {
    if (k < 1) throw std::invalid_argument("k must be >= 1");
    k_ = k;
    initialized_ = false;
}

inline void FMMCoulomb::cartesian_to_spherical(const Vec3& r, Scalar& rr,
                                                Scalar& costh, Scalar& phi) const {
    rr = r.norm();
    if (rr < 1e-20) {
        costh = 1.0;
        phi = 0.0;
    } else {
        costh = r[2] / rr;
        phi = std::atan2(r[1], r[0]);
    }
}

inline void FMMCoulomb::compute_solid_harmonic_R(
    Scalar r, Scalar costh, Scalar phi,
    std::vector<Scalar>& Rl0,
    std::vector<std::complex<Scalar>>& Rlm) const {

    // R_l^m = r^l * Y_l^m (unnormalized solid harmonics)
    // Using recurrence relations for associated Legendre polynomials

    Rl0.assign(l_max_ + 1, 0.0);
    Rlm.assign(n_off_diag(), std::complex<Scalar>(0.0, 0.0));

    Scalar sinth = std::sqrt(std::max(0.0, 1.0 - costh * costh));

    // P_l^m for Legendre polynomials
    std::vector<std::vector<Scalar>> P(l_max_ + 1, std::vector<Scalar>(l_max_ + 1, 0.0));

    // P_0^0 = 1
    P[0][0] = 1.0;

    // P_l^l = -(2l-1) * sinth * P_{l-1}^{l-1}
    for (int l = 1; l <= l_max_; ++l) {
        P[l][l] = -(2 * l - 1) * sinth * P[l - 1][l - 1];
    }

    // P_l^{l-1} = (2l-1) * costh * P_{l-1}^{l-1}
    for (int l = 1; l <= l_max_; ++l) {
        P[l][l - 1] = (2 * l - 1) * costh * P[l - 1][l - 1];
    }

    // Recurrence: (l-m) P_l^m = (2l-1) costh P_{l-1}^m - (l+m-1) P_{l-2}^m
    for (int m = 0; m < l_max_ - 1; ++m) {
        for (int l = m + 2; l <= l_max_; ++l) {
            P[l][m] = ((2 * l - 1) * costh * P[l - 1][m] - (l + m - 1) * P[l - 2][m]) / (l - m);
        }
    }

    // Compute R_l^m
    Scalar r_power = 1.0;
    for (int l = 0; l <= l_max_; ++l) {
        // m = 0
        Rl0[l] = r_power * P[l][0];

        // m > 0
        for (int m = 1; m <= l; ++m) {
            Scalar phase = std::exp(std::complex<Scalar>(0.0, m * phi)).real();
            Scalar phase_i = std::exp(std::complex<Scalar>(0.0, m * phi)).imag();
            Rlm[lm_index(l, m)] = r_power * P[l][m] *
                std::complex<Scalar>(phase, phase_i);
        }

        r_power *= r;
    }
}

inline void FMMCoulomb::compute_solid_harmonic_I(
    Scalar r, Scalar costh, Scalar phi,
    int max_l,
    std::vector<Scalar>& Il0,
    std::vector<std::complex<Scalar>>& Ilm) const {

    // I_l^m = r^{-(l+1)} * Y_l^m (irregular solid harmonics)
    int n_off = max_l * (max_l + 1) / 2;
    Il0.assign(max_l + 1, 0.0);
    Ilm.assign(n_off, std::complex<Scalar>(0.0, 0.0));

    if (r < 1e-20) return;

    Scalar sinth = std::sqrt(std::max(0.0, 1.0 - costh * costh));

    // Legendre polynomials
    std::vector<std::vector<Scalar>> P(max_l + 1, std::vector<Scalar>(max_l + 1, 0.0));
    P[0][0] = 1.0;

    for (int l = 1; l <= max_l; ++l) {
        P[l][l] = -(2 * l - 1) * sinth * P[l - 1][l - 1];
    }

    for (int l = 1; l <= max_l; ++l) {
        if (l - 1 >= 0 && l <= max_l) {
            P[l][l - 1] = (2 * l - 1) * costh * P[l - 1][l - 1];
        }
    }

    for (int m = 0; m < max_l - 1; ++m) {
        for (int l = m + 2; l <= max_l; ++l) {
            P[l][m] = ((2 * l - 1) * costh * P[l - 1][m] - (l + m - 1) * P[l - 2][m]) / (l - m);
        }
    }

    // Compute I_l^m = r^{-(l+1)} * P_l^m * exp(i*m*phi)
    Scalar r_inv = 1.0 / r;
    Scalar r_inv_power = r_inv;  // Start with r^{-1}

    for (int l = 0; l <= max_l; ++l) {
        Il0[l] = r_inv_power * P[l][0];

        for (int m = 1; m <= l; ++m) {
            int idx = l * (l - 1) / 2 + m - 1;
            Scalar phase_r = std::cos(m * phi);
            Scalar phase_i = std::sin(m * phi);
            Ilm[idx] = r_inv_power * P[l][m] * std::complex<Scalar>(phase_r, phase_i);
        }

        r_inv_power *= r_inv;
    }
}

inline void FMMCoulomb::multipole_to_multipole(
    const Vec3& dr, int l_max,
    const std::vector<Scalar>& Ml0_child,
    const std::vector<std::complex<Scalar>>& Mlm_child,
    std::vector<Scalar>& Ml0,
    std::vector<std::complex<Scalar>>& Mlm) const {

    // M2M translation: shifts multipole expansion from child center to parent center

    Scalar rr, costh, phi;
    cartesian_to_spherical(dr, rr, costh, phi);

    std::vector<Scalar> Rl0;
    std::vector<std::complex<Scalar>> Rlm;
    compute_solid_harmonic_R(rr, costh, phi, Rl0, Rlm);

    // For each target (l,m) in parent
    for (int l = 0; l <= l_max; ++l) {
        // m = 0
        for (int lambda = 0; lambda <= l; ++lambda) {
            Ml0[l] += Ml0_child[l - lambda] * Rl0[lambda];

            for (int mu = 1; mu <= std::min(lambda, l - lambda); ++mu) {
                int sign = ((mu % 2) == 0) ? 1 : -1;
                Ml0[l] += 2 * sign * std::real(
                    Mlm_child[lm_index(l - lambda, mu)] *
                    Rlm[lm_index(lambda, mu)]);
            }
        }

        // m > 0
        for (int m = 1; m <= l; ++m) {
            int j = lm_index(l, m);

            for (int lambda = 0; lambda <= l; ++lambda) {
                for (int mu = std::max(-lambda, lambda - l + m);
                     mu <= std::min(lambda, -lambda + l + m); ++mu) {

                    std::complex<Scalar> cur_M, cur_R;

                    // Get M_{l-lambda}^{m-mu}
                    int lm = l - lambda;
                    int mm = m - mu;
                    if (mm < 0) {
                        int sign = ((-mm) % 2 == 0) ? 1 : -1;
                        cur_M = static_cast<Scalar>(sign) * std::conj(Mlm_child[lm_index(lm, -mm)]);
                    } else if (mm == 0) {
                        cur_M = Ml0_child[lm];
                    } else {
                        cur_M = Mlm_child[lm_index(lm, mm)];
                    }

                    // Get R_lambda^mu
                    if (mu < 0) {
                        int sign = ((-mu) % 2 == 0) ? 1 : -1;
                        cur_R = static_cast<Scalar>(sign) * std::conj(Rlm[lm_index(lambda, -mu)]);
                    } else if (mu == 0) {
                        cur_R = Rl0[lambda];
                    } else {
                        cur_R = Rlm[lm_index(lambda, mu)];
                    }

                    Mlm[j] += cur_M * std::conj(cur_R);
                }
            }
        }
    }
}

inline void FMMCoulomb::multipole_to_local(
    const Vec3& dr, int l_max,
    const std::vector<Scalar>& Ml0,
    const std::vector<std::complex<Scalar>>& Mlm,
    std::vector<Scalar>& Ll0,
    std::vector<std::complex<Scalar>>& Llm) const {

    // M2L translation: converts multipole expansion to local expansion

    Scalar rr, costh, phi;
    cartesian_to_spherical(dr, rr, costh, phi);

    std::vector<Scalar> Il0;
    std::vector<std::complex<Scalar>> Ilm;
    compute_solid_harmonic_I(rr, costh, phi, 2 * l_max, Il0, Ilm);

    for (int l = 0; l <= l_max; ++l) {
        // m = 0
        for (int lambda = 0; lambda <= l_max; ++lambda) {
            Ll0[l] += Ml0[lambda] * Il0[l + lambda];

            for (int mu = 1; mu <= lambda; ++mu) {
                int idx_m = lm_index(lambda, mu);
                int idx_i = (l + lambda) * (l + lambda - 1) / 2 + mu - 1;
                Ll0[l] += 2 * std::real(Mlm[idx_m] * Ilm[idx_i]);
            }
        }

        // m > 0
        for (int m = 1; m <= l; ++m) {
            int k = lm_index(l, m);

            for (int lambda = 0; lambda <= l_max; ++lambda) {
                for (int mu = -lambda; mu <= lambda; ++mu) {
                    std::complex<Scalar> cur_M, cur_I;

                    // Get I_{l+lambda}^{m+mu}
                    int li = l + lambda;
                    int mi = m + mu;
                    if (mi < 0) {
                        int sign = ((-mi) % 2 == 0) ? 1 : -1;
                        int idx = li * (li - 1) / 2 + (-mi) - 1;
                        cur_I = static_cast<Scalar>(sign) * std::conj(Ilm[idx]);
                    } else if (mi == 0) {
                        cur_I = Il0[li];
                    } else {
                        int idx = li * (li - 1) / 2 + mi - 1;
                        cur_I = Ilm[idx];
                    }

                    // Get M_lambda^mu
                    if (mu < 0) {
                        int sign = ((-mu) % 2 == 0) ? 1 : -1;
                        cur_M = static_cast<Scalar>(sign) * std::conj(Mlm[lm_index(lambda, -mu)]);
                    } else if (mu == 0) {
                        cur_M = Ml0[lambda];
                    } else {
                        cur_M = Mlm[lm_index(lambda, mu)];
                    }

                    Llm[k] += cur_M * cur_I;
                }
            }
        }
    }
}

inline void FMMCoulomb::local_to_local(
    const Vec3& dr, int l_max_in,
    const std::vector<Scalar>& Ll0_in,
    const std::vector<std::complex<Scalar>>& Llm_in,
    int l_max_out,
    std::vector<Scalar>& Ll0_out,
    std::vector<std::complex<Scalar>>& Llm_out) const {

    // L2L translation: shifts local expansion from parent center to child center

    Scalar rr, costh, phi;
    cartesian_to_spherical(dr, rr, costh, phi);

    std::vector<Scalar> Rl0;
    std::vector<std::complex<Scalar>> Rlm;
    compute_solid_harmonic_R(rr, costh, phi, Rl0, Rlm);

    // Conjugate Rlm
    for (auto& val : Rlm) {
        val = std::conj(val);
    }

    for (int l = 0; l <= l_max_out; ++l) {
        // m = 0
        for (int lambda = 0; lambda <= l_max_in - l; ++lambda) {
            Ll0_out[l] += Ll0_in[l + lambda] * Rl0[lambda];

            for (int mu = 1; mu <= lambda; ++mu) {
                int idx_l = lm_index(l + lambda, mu);
                int idx_r = lm_index(lambda, mu);
                Ll0_out[l] += 2 * std::real(Llm_in[idx_l] * Rlm[idx_r]);
            }
        }

        // m > 0
        for (int m = 1; m <= l; ++m) {
            int k = lm_index(l, m);

            for (int lambda = 0; lambda <= l_max_in - l; ++lambda) {
                for (int mu = -lambda; mu <= lambda; ++mu) {
                    std::complex<Scalar> cur_L, cur_R;

                    // Get L_{l+lambda}^{m+mu}
                    int ll = l + lambda;
                    int mm = m + mu;
                    if (mm < 0) {
                        int sign = ((-mm) % 2 == 0) ? 1 : -1;
                        cur_L = static_cast<Scalar>(sign) * std::conj(Llm_in[lm_index(ll, -mm)]);
                    } else if (mm == 0) {
                        cur_L = Ll0_in[ll];
                    } else {
                        cur_L = Llm_in[lm_index(ll, mm)];
                    }

                    // Get R_lambda^mu
                    if (mu < 0) {
                        int sign = ((-mu) % 2 == 0) ? 1 : -1;
                        cur_R = static_cast<Scalar>(sign) * std::conj(Rlm[lm_index(lambda, -mu)]);
                    } else if (mu == 0) {
                        cur_R = Rl0[lambda];
                    } else {
                        cur_R = Rlm[lm_index(lambda, mu)];
                    }

                    Llm_out[k] += cur_L * cur_R;
                }
            }
        }
    }
}

inline Vec3 FMMCoulomb::cell_center(int x, int y, int z, int level) const {
    const auto& n = tree_[level].n;
    Vec3 h = (upper_ - lower_);
    h[0] /= n[0];
    h[1] /= n[1];
    h[2] /= n[2];

    Vec3 center;
    center[0] = lower_[0] + (x + 0.5) * h[0];
    center[1] = lower_[1] + (y + 0.5) * h[1];
    center[2] = lower_[2] + (z + 0.5) * h[2];

    return center;
}

inline void FMMCoulomb::initialize_tree(AtomicSystem& system) {
    // Set domain bounds
    const Mat3& cell = system.cell();
    lower_ = Vec3::Zero();
    upper_[0] = cell(0, 0);
    upper_[1] = cell(1, 1);
    upper_[2] = cell(2, 2);

    // Initialize tree levels
    tree_.resize(n_level_);

    for (int level = 0; level < n_level_; ++level) {
        if (level == 0) {
            tree_[0].n = n0_;
        } else {
            tree_[level].n[0] = 2 * tree_[level - 1].n[0];
            tree_[level].n[1] = 2 * tree_[level - 1].n[1];
            tree_[level].n[2] = 2 * tree_[level - 1].n[2];
        }

        tree_[level].n_tot = tree_[level].n[0] * tree_[level].n[1] * tree_[level].n[2];

        // Allocate arrays
        tree_[level].np.resize(tree_[level].n_tot, 0);

        tree_[level].Ml0.resize(tree_[level].n_tot,
                                std::vector<Scalar>(l_max_ + 1, 0.0));
        tree_[level].Mlm.resize(tree_[level].n_tot,
                                std::vector<std::complex<Scalar>>(n_off_diag(),
                                    std::complex<Scalar>(0.0, 0.0)));

        tree_[level].Ll0.resize(tree_[level].n_tot,
                                std::vector<Scalar>(l_max_ + 1, 0.0));
        tree_[level].Llm.resize(tree_[level].n_tot,
                                std::vector<std::complex<Scalar>>(n_off_diag(),
                                    std::complex<Scalar>(0.0, 0.0)));
    }

    // Leaf-level particle lists
    int n_leaves = tree_[n_level_ - 1].n_tot;
    np_.resize(n_leaves, 0);
    p_in_leaf_.resize(n_leaves);

    // Root multipole/local
    root_Ml0_.resize(l_max_ + 1, 0.0);
    root_Mlm_.resize(n_off_diag(), std::complex<Scalar>(0.0, 0.0));
    root_Ll0_.resize(l_max_ + 1, 0.0);
    root_Llm_.resize(n_off_diag(), std::complex<Scalar>(0.0, 0.0));

    initialized_ = true;
}

inline void FMMCoulomb::sort_particles(AtomicSystem& system) {
    const std::size_t num_atoms = system.num_atoms();
    const auto& n = tree_[n_level_ - 1].n;

    // Clear particle lists
    for (auto& list : p_in_leaf_) {
        list.clear();
    }
    std::fill(np_.begin(), np_.end(), 0);

    Vec3 h = upper_ - lower_;
    h[0] /= n[0];
    h[1] /= n[1];
    h[2] /= n[2];

    for (std::size_t i = 0; i < num_atoms; ++i) {
        Vec3 pos = system.position(i).matrix();

        // Compute cell indices
        int ix = static_cast<int>((pos[0] - lower_[0]) / h[0]);
        int iy = static_cast<int>((pos[1] - lower_[1]) / h[1]);
        int iz = static_cast<int>((pos[2] - lower_[2]) / h[2]);

        // Clamp to valid range
        ix = std::max(0, std::min(ix, n[0] - 1));
        iy = std::max(0, std::min(iy, n[1] - 1));
        iz = std::max(0, std::min(iz, n[2] - 1));

        int idx = cell_index(ix, iy, iz, n);
        p_in_leaf_[idx].push_back(static_cast<int>(i));
        np_[idx]++;
    }
}

inline void FMMCoulomb::M2M_pass(AtomicSystem& system) {
    // Compute multipoles at leaf level
    int level = n_level_ - 1;
    const auto& n = tree_[level].n;

    // Clear multipoles
    for (auto& Ml0 : tree_[level].Ml0) {
        std::fill(Ml0.begin(), Ml0.end(), 0.0);
    }
    for (auto& Mlm : tree_[level].Mlm) {
        std::fill(Mlm.begin(), Mlm.end(), std::complex<Scalar>(0.0, 0.0));
    }

    for (int x = 0; x < n[0]; ++x) {
        for (int y = 0; y < n[1]; ++y) {
            for (int z = 0; z < n[2]; ++z) {
                int j = cell_index(x, y, z, n);
                Vec3 r0 = cell_center(x, y, z, level);

                for (int ii = 0; ii < np_[j]; ++ii) {
                    int i = p_in_leaf_[j][ii];
                    Vec3 pos = system.position(i).matrix();
                    Vec3 dr = pos - r0;
                    Scalar q = charges_[i];

                    if (dr.squaredNorm() < 1e-20) {
                        tree_[level].Ml0[j][0] += q;
                    } else {
                        Scalar rr, costh, phi;
                        cartesian_to_spherical(dr, rr, costh, phi);

                        std::vector<Scalar> Rl0;
                        std::vector<std::complex<Scalar>> Rlm;
                        compute_solid_harmonic_R(rr, costh, phi, Rl0, Rlm);

                        for (int l = 0; l <= l_max_; ++l) {
                            tree_[level].Ml0[j][l] += q * Rl0[l];
                        }
                        for (int idx = 0; idx < n_off_diag(); ++idx) {
                            tree_[level].Mlm[j][idx] += q * std::conj(Rlm[idx]);
                        }
                    }
                }

                tree_[level].np[j] = np_[j];
            }
        }
    }

    // Propagate upward
    for (int level = n_level_ - 2; level >= 0; --level) {
        const auto& n = tree_[level].n;
        const auto& n_child = tree_[level + 1].n;

        // Clear multipoles
        for (auto& Ml0 : tree_[level].Ml0) {
            std::fill(Ml0.begin(), Ml0.end(), 0.0);
        }
        for (auto& Mlm : tree_[level].Mlm) {
            std::fill(Mlm.begin(), Mlm.end(), std::complex<Scalar>(0.0, 0.0));
        }
        for (auto& np : tree_[level].np) {
            np = 0;
        }

        for (int x = 0; x < n[0]; ++x) {
            for (int y = 0; y < n[1]; ++y) {
                for (int z = 0; z < n[2]; ++z) {
                    int j = cell_index(x, y, z, n);
                    Vec3 r0 = cell_center(x, y, z, level);

                    int np = 0;

                    // Sum contributions from 8 children
                    for (int dx = 0; dx < 2; ++dx) {
                        for (int dy = 0; dy < 2; ++dy) {
                            for (int dz = 0; dz < 2; ++dz) {
                                int xc = 2 * x + dx;
                                int yc = 2 * y + dy;
                                int zc = 2 * z + dz;
                                int jc = cell_index(xc, yc, zc, n_child);

                                if (tree_[level + 1].np[jc] > 0) {
                                    Vec3 rc = cell_center(xc, yc, zc, level + 1);
                                    Vec3 dr = rc - r0;

                                    if (dr.squaredNorm() < 1e-20) {
                                        for (int l = 0; l <= l_max_; ++l) {
                                            tree_[level].Ml0[j][l] +=
                                                tree_[level + 1].Ml0[jc][l];
                                        }
                                        for (int idx = 0; idx < n_off_diag(); ++idx) {
                                            tree_[level].Mlm[j][idx] +=
                                                tree_[level + 1].Mlm[jc][idx];
                                        }
                                    } else {
                                        multipole_to_multipole(
                                            dr, l_max_,
                                            tree_[level + 1].Ml0[jc],
                                            tree_[level + 1].Mlm[jc],
                                            tree_[level].Ml0[j],
                                            tree_[level].Mlm[j]);
                                    }

                                    np += tree_[level + 1].np[jc];
                                }
                            }
                        }
                    }

                    tree_[level].np[j] = np;
                }
            }
        }
    }

    // Compute root multipole
    std::fill(root_Ml0_.begin(), root_Ml0_.end(), 0.0);
    std::fill(root_Mlm_.begin(), root_Mlm_.end(), std::complex<Scalar>(0.0, 0.0));

    Vec3 r0 = (lower_ + upper_) / 2.0;
    const auto& n0 = tree_[0].n;

    for (int x = 0; x < n0[0]; ++x) {
        for (int y = 0; y < n0[1]; ++y) {
            for (int z = 0; z < n0[2]; ++z) {
                int j = cell_index(x, y, z, n0);
                Vec3 r = cell_center(x, y, z, 0);
                Vec3 dr = r - r0;

                if (dr.squaredNorm() < 1e-20) {
                    for (int l = 0; l <= l_max_; ++l) {
                        root_Ml0_[l] += tree_[0].Ml0[j][l];
                    }
                    for (int idx = 0; idx < n_off_diag(); ++idx) {
                        root_Mlm_[idx] += tree_[0].Mlm[j][idx];
                    }
                } else {
                    multipole_to_multipole(dr, l_max_,
                        tree_[0].Ml0[j], tree_[0].Mlm[j],
                        root_Ml0_, root_Mlm_);
                }
            }
        }
    }
}

inline void FMMCoulomb::M2L_pass(AtomicSystem& system) {
    // Initialize local expansions
    std::fill(root_Ll0_.begin(), root_Ll0_.end(), 0.0);
    std::fill(root_Llm_.begin(), root_Llm_.end(), std::complex<Scalar>(0.0, 0.0));

    for (int level = 0; level < n_level_; ++level) {
        for (auto& Ll0 : tree_[level].Ll0) {
            std::fill(Ll0.begin(), Ll0.end(), 0.0);
        }
        for (auto& Llm : tree_[level].Llm) {
            std::fill(Llm.begin(), Llm.end(), std::complex<Scalar>(0.0, 0.0));
        }
    }

    // M2L at each level
    // For level 0, interact with cells in interaction list (well-separated)
    for (int level = 0; level < n_level_; ++level) {
        const auto& n = tree_[level].n;
        Vec3 h = upper_ - lower_;
        h[0] /= n[0];
        h[1] /= n[1];
        h[2] /= n[2];

        // For each cell, find well-separated cells and do M2L
        for (int x = 0; x < n[0]; ++x) {
            for (int y = 0; y < n[1]; ++y) {
                for (int z = 0; z < n[2]; ++z) {
                    int i = cell_index(x, y, z, n);

                    if (tree_[level].np[i] == 0) continue;

                    // Well-separated cells: |dx|>1 or |dy|>1 or |dz|>1
                    // but parent cells were neighbors (for level > 0)
                    for (int dx = -3; dx <= 3; ++dx) {
                        for (int dy = -3; dy <= 3; ++dy) {
                            for (int dz = -3; dz <= 3; ++dz) {
                                // Skip near-neighbors
                                if (std::abs(dx) <= 1 && std::abs(dy) <= 1 &&
                                    std::abs(dz) <= 1) continue;

                                // Handle periodicity
                                int sx = x + dx;
                                int sy = y + dy;
                                int sz = z + dz;

                                // For non-periodic, skip out-of-bounds
                                if (!system.pbc()[0] && (sx < 0 || sx >= n[0])) continue;
                                if (!system.pbc()[1] && (sy < 0 || sy >= n[1])) continue;
                                if (!system.pbc()[2] && (sz < 0 || sz >= n[2])) continue;

                                // Wrap for periodic
                                int wx = ((sx % n[0]) + n[0]) % n[0];
                                int wy = ((sy % n[1]) + n[1]) % n[1];
                                int wz = ((sz % n[2]) + n[2]) % n[2];

                                int j = cell_index(wx, wy, wz, n);

                                if (tree_[level].np[j] == 0) continue;

                                // Compute displacement
                                Vec3 dr;
                                dr[0] = -dx * h[0];
                                dr[1] = -dy * h[1];
                                dr[2] = -dz * h[2];

                                multipole_to_local(dr, l_max_,
                                    tree_[level].Ml0[j],
                                    tree_[level].Mlm[j],
                                    tree_[level].Ll0[i],
                                    tree_[level].Llm[i]);
                            }
                        }
                    }
                }
            }
        }
    }
}

inline void FMMCoulomb::L2L_pass() {
    // First, add root local expansion to level 0
    Vec3 r0 = (lower_ + upper_) / 2.0;
    const auto& n0 = tree_[0].n;

    for (int x = 0; x < n0[0]; ++x) {
        for (int y = 0; y < n0[1]; ++y) {
            for (int z = 0; z < n0[2]; ++z) {
                int i = cell_index(x, y, z, n0);

                if (tree_[0].np[i] == 0) continue;

                Vec3 r = cell_center(x, y, z, 0);
                Vec3 dr = r0 - r;

                if (dr.squaredNorm() < 1e-20) {
                    for (int l = 0; l <= l_max_; ++l) {
                        tree_[0].Ll0[i][l] += root_Ll0_[l];
                    }
                    for (int idx = 0; idx < n_off_diag(); ++idx) {
                        tree_[0].Llm[i][idx] += root_Llm_[idx];
                    }
                } else {
                    local_to_local(dr, l_max_, root_Ll0_, root_Llm_,
                                  l_max_, tree_[0].Ll0[i], tree_[0].Llm[i]);
                }
            }
        }
    }

    // Propagate downward
    for (int level = 1; level < n_level_; ++level) {
        const auto& n = tree_[level].n;
        const auto& n_parent = tree_[level - 1].n;

        for (int x = 0; x < n_parent[0]; ++x) {
            for (int y = 0; y < n_parent[1]; ++y) {
                for (int z = 0; z < n_parent[2]; ++z) {
                    int j = cell_index(x, y, z, n_parent);

                    if (tree_[level - 1].np[j] == 0) continue;

                    Vec3 r0 = cell_center(x, y, z, level - 1);

                    // Pass to 8 children
                    for (int dx = 0; dx < 2; ++dx) {
                        for (int dy = 0; dy < 2; ++dy) {
                            for (int dz = 0; dz < 2; ++dz) {
                                int xc = 2 * x + dx;
                                int yc = 2 * y + dy;
                                int zc = 2 * z + dz;
                                int i = cell_index(xc, yc, zc, n);

                                if (tree_[level].np[i] == 0) continue;

                                Vec3 r = cell_center(xc, yc, zc, level);
                                Vec3 dr = r0 - r;

                                if (dr.squaredNorm() < 1e-20) {
                                    for (int l = 0; l <= l_max_; ++l) {
                                        tree_[level].Ll0[i][l] +=
                                            tree_[level - 1].Ll0[j][l];
                                    }
                                    for (int idx = 0; idx < n_off_diag(); ++idx) {
                                        tree_[level].Llm[i][idx] +=
                                            tree_[level - 1].Llm[j][idx];
                                    }
                                } else {
                                    local_to_local(dr, l_max_,
                                        tree_[level - 1].Ll0[j],
                                        tree_[level - 1].Llm[j],
                                        l_max_,
                                        tree_[level].Ll0[i],
                                        tree_[level].Llm[i]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

inline Scalar FMMCoulomb::near_field_and_local(
    AtomicSystem& system,
    bool compute_forces,
    bool compute_virial,
    Mat3& virial) {

    Scalar epot = 0.0;
    int level = n_level_ - 1;
    const auto& n = tree_[level].n;

    Vec3 h = upper_ - lower_;
    h[0] /= n[0];
    h[1] /= n[1];
    h[2] /= n[2];

    for (int x = 0; x < n[0]; ++x) {
        for (int y = 0; y < n[1]; ++y) {
            for (int z = 0; z < n[2]; ++z) {
                int i_cell = cell_index(x, y, z, n);

                if (np_[i_cell] == 0) continue;

                Vec3 r0 = cell_center(x, y, z, level);

                // Process each particle in this cell
                for (int ii = 0; ii < np_[i_cell]; ++ii) {
                    int i = p_in_leaf_[i_cell][ii];
                    Vec3 ri = system.position(i).matrix();
                    Scalar qi = charges_[i];

                    // Contribution from local expansion
                    Vec3 dr_local = r0 - ri;
                    std::vector<Scalar> Ll0_atom(2, 0.0);
                    std::vector<std::complex<Scalar>> Llm_atom(1, 0.0);

                    local_to_local(dr_local, l_max_,
                                  tree_[level].Ll0[i_cell],
                                  tree_[level].Llm[i_cell],
                                  1, Ll0_atom, Llm_atom);

                    // Potential from far-field
                    Scalar phi_far = COULOMB_CONST * Ll0_atom[0];
                    epot += 0.5 * qi * phi_far;

                    // Force from far-field: E = -grad(phi)
                    // E = -(dL/dx, dL/dy, dL/dz) where L is local expansion
                    // For l=1: L_1^0 = E_z, L_1^1 = E_x - i*E_y
                    if (compute_forces) {
                        Vec3 E_far;
                        E_far[0] = -std::real(Llm_atom[0]);
                        E_far[1] = -std::imag(Llm_atom[0]);
                        E_far[2] = Ll0_atom[1];

                        Vec3 force = COULOMB_CONST * qi * E_far;
                        system.forces().col(i) += force.array();
                    }

                    // Self-interaction within cell (j > ii to avoid double counting)
                    for (int jj = ii + 1; jj < np_[i_cell]; ++jj) {
                        int j = p_in_leaf_[i_cell][jj];
                        Vec3 rj = system.position(j).matrix();
                        Scalar qj = charges_[j];

                        Vec3 dr = ri - rj;
                        Scalar r_sq = dr.squaredNorm();

                        if (r_sq < 1e-20) continue;

                        Scalar r = std::sqrt(r_sq);
                        Scalar inv_r = 1.0 / r;

                        Scalar pair_energy = COULOMB_CONST * qi * qj * inv_r;
                        epot += pair_energy;

                        if (compute_forces || compute_virial) {
                            Scalar force_mag = pair_energy * inv_r * inv_r;
                            Vec3 force = force_mag * dr;

                            if (compute_forces) {
                                system.forces().col(i) += force.array();
                                system.forces().col(j) -= force.array();
                            }

                            if (compute_virial) {
                                virial += dr * force.transpose();
                            }
                        }
                    }

                    // Neighbors (only half for non-periodic to avoid double counting)
                    for (int neigh = 0; neigh < NUM_NEIGHBORS / 2; ++neigh) {
                        int tx = x + neighbor_offsets_[neigh][0];
                        int ty = y + neighbor_offsets_[neigh][1];
                        int tz = z + neighbor_offsets_[neigh][2];

                        // Handle boundaries
                        bool valid = true;
                        Vec3 shift = Vec3::Zero();

                        if (system.pbc()[0]) {
                            if (tx < 0) { tx += n[0]; shift[0] = -(upper_[0] - lower_[0]); }
                            if (tx >= n[0]) { tx -= n[0]; shift[0] = upper_[0] - lower_[0]; }
                        } else {
                            if (tx < 0 || tx >= n[0]) valid = false;
                        }

                        if (system.pbc()[1]) {
                            if (ty < 0) { ty += n[1]; shift[1] = -(upper_[1] - lower_[1]); }
                            if (ty >= n[1]) { ty -= n[1]; shift[1] = upper_[1] - lower_[1]; }
                        } else {
                            if (ty < 0 || ty >= n[1]) valid = false;
                        }

                        if (system.pbc()[2]) {
                            if (tz < 0) { tz += n[2]; shift[2] = -(upper_[2] - lower_[2]); }
                            if (tz >= n[2]) { tz -= n[2]; shift[2] = upper_[2] - lower_[2]; }
                        } else {
                            if (tz < 0 || tz >= n[2]) valid = false;
                        }

                        if (!valid) continue;

                        int j_cell = cell_index(tx, ty, tz, n);

                        for (int jj = 0; jj < np_[j_cell]; ++jj) {
                            int j = p_in_leaf_[j_cell][jj];
                            Vec3 rj = system.position(j).matrix() + shift;
                            Scalar qj = charges_[j];

                            Vec3 dr = ri - rj;
                            Scalar r_sq = dr.squaredNorm();

                            if (r_sq < 1e-20) continue;

                            Scalar r = std::sqrt(r_sq);
                            Scalar inv_r = 1.0 / r;

                            Scalar pair_energy = COULOMB_CONST * qi * qj * inv_r;
                            epot += pair_energy;

                            if (compute_forces || compute_virial) {
                                Scalar force_mag = pair_energy * inv_r * inv_r;
                                Vec3 force = force_mag * dr;

                                if (compute_forces) {
                                    system.forces().col(i) += force.array();
                                    system.forces().col(j) -= force.array();
                                }

                                if (compute_virial) {
                                    virial += dr * force.transpose();
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return epot;
}

inline PotentialResults FMMCoulomb::compute_impl(
    AtomicSystem& system,
    NeighborList& /* neighbors */,  // FMM doesn't use neighbor list
    bool compute_forces,
    bool compute_virial) {

    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();

    if (charges_.size() != num_atoms) {
        throw std::runtime_error("FMMCoulomb: charges array size mismatch");
    }

    // Initialize tree if needed
    if (!initialized_) {
        initialize_tree(system);
    }

    // Sort particles into tree
    sort_particles(system);

    // FMM algorithm
    M2M_pass(system);
    M2L_pass(system);
    L2L_pass();

    // Near-field and final evaluation
    Mat3 virial = Mat3::Zero();
    Scalar energy = near_field_and_local(system, compute_forces, compute_virial, virial);

    results.energy = energy;
    results.virial = virial;

    return results;
}

} // namespace atomistica
