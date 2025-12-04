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

#include "../../config.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../potential_base.hpp"
#include "coulomb.hpp"

namespace atomistica {

/**
 * @brief Particle Mesh Ewald (PME) for long-range Coulomb interactions
 *
 * Implements the smooth Particle Mesh Ewald method for periodic systems.
 *
 * The total Coulomb energy is split into:
 *   E_total = E_real + E_reciprocal - E_self
 *
 * where:
 *   E_real = sum_{i<j} q_i*q_j*erfc(alpha*r_ij)/r_ij  (short-range, cutoff)
 *   E_reciprocal = k-space sum via FFT
 *   E_self = sum_i q_i^2 * sqrt(alpha/pi)  (self-interaction correction)
 *
 * References:
 * - Darden, York, Pedersen, J. Chem. Phys. 98, 10089 (1993)
 * - Essmann et al., J. Chem. Phys. 103, 8577 (1995)
 *
 * @note Requires 3D periodic boundary conditions
 */
class PMECoulomb : public PotentialBase<PMECoulomb> {
public:
    /**
     * @brief Construct PME Coulomb solver
     * @param cutoff Real-space cutoff (Angstrom)
     * @param grid_x Grid size in x direction (should be FFT-friendly: 2^a*3^b*5^c)
     * @param grid_y Grid size in y direction
     * @param grid_z Grid size in z direction
     * @param order B-spline interpolation order (4, 6, or 8 recommended)
     * @param alpha Ewald parameter (0 = auto-compute from cutoff)
     */
    explicit PMECoulomb(Scalar cutoff = 10.0,
                        int grid_x = 32, int grid_y = 32, int grid_z = 32,
                        int order = 4, Scalar alpha = 0.0);

    ~PMECoulomb();

    // Disable copy (FFT plans are not copyable)
    PMECoulomb(const PMECoulomb&) = delete;
    PMECoulomb& operator=(const PMECoulomb&) = delete;

    // Move operations
    PMECoulomb(PMECoulomb&& other) noexcept;
    PMECoulomb& operator=(PMECoulomb&& other) noexcept;

    /**
     * @brief Set real-space cutoff
     */
    void set_cutoff(Scalar cutoff);

    /**
     * @brief Set Ewald parameter alpha
     * @param alpha If <= 0, auto-computed from cutoff
     */
    void set_alpha(Scalar alpha);

    /**
     * @brief Set grid dimensions
     */
    void set_grid(int grid_x, int grid_y, int grid_z);

    /**
     * @brief Set B-spline order
     */
    void set_order(int order);

    // Getters
    Scalar alpha() const { return alpha_; }
    int order() const { return order_; }
    std::array<int, 3> grid() const { return grid_; }

    void set_charges(const std::vector<Scalar>& charges) {
        charges_ = charges;
    }

    const std::vector<Scalar>& charges() const { return charges_; }

    // CRTP implementation
    Scalar cutoff_impl() const { return cutoff_; }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial);

private:
    // Initialize/update internal structures
    void initialize();
    void compute_bspline_moduli();

    // B-spline functions
    void fill_bspline(Scalar w, std::vector<Scalar>& theta,
                      std::vector<Scalar>& dtheta) const;

    // Real-space contribution
    void compute_real_space(AtomicSystem& system,
                           NeighborList& neighbors,
                           bool compute_forces,
                           bool compute_virial,
                           Scalar& energy,
                           Mat3& virial);

    // Reciprocal space contribution
    void compute_reciprocal_space(AtomicSystem& system,
                                  bool compute_forces,
                                  bool compute_virial,
                                  Scalar& energy,
                                  Mat3& virial);

    // Self-energy correction
    Scalar compute_self_energy() const;

    // Simple 3D FFT (not optimized, but portable)
    void fft_forward();
    void fft_backward();

    // Parameters
    Scalar cutoff_ = 10.0;
    Scalar cutoff_sq_;
    Scalar alpha_ = 0.0;
    Scalar sqrt_alpha_;
    Scalar sqrt_alpha_pi_;
    std::array<int, 3> grid_ = {32, 32, 32};
    int order_ = 4;

    std::vector<Scalar> charges_;

    // B-spline moduli for structure factor correction
    std::vector<Scalar> bsp_mod_x_;
    std::vector<Scalar> bsp_mod_y_;
    std::vector<Scalar> bsp_mod_z_;

    // Per-atom B-spline coefficients
    std::vector<Scalar> theta_x_, theta_y_, theta_z_;
    std::vector<Scalar> dtheta_x_, dtheta_y_, dtheta_z_;

    // Scaled fractional coordinates
    std::vector<Scalar> fr_x_, fr_y_, fr_z_;

    // Charge grid (complex for FFT)
    std::vector<std::complex<Scalar>> Q_;

    // Initialization flag
    bool initialized_ = false;
};

// ============================================================================
// Implementation
// ============================================================================

inline PMECoulomb::PMECoulomb(Scalar cutoff, int grid_x, int grid_y, int grid_z,
                              int order, Scalar alpha)
    : cutoff_(cutoff)
    , cutoff_sq_(cutoff * cutoff)
    , grid_{grid_x, grid_y, grid_z}
    , order_(order) {
    set_alpha(alpha);
}

inline PMECoulomb::~PMECoulomb() = default;

inline PMECoulomb::PMECoulomb(PMECoulomb&& other) noexcept = default;
inline PMECoulomb& PMECoulomb::operator=(PMECoulomb&& other) noexcept = default;

inline void PMECoulomb::set_cutoff(Scalar cutoff) {
    cutoff_ = cutoff;
    cutoff_sq_ = cutoff * cutoff;
    // Recompute alpha if auto-computed
    if (alpha_ <= 0) {
        set_alpha(0.0);
    }
    initialized_ = false;
}

inline void PMECoulomb::set_alpha(Scalar alpha) {
    if (alpha <= 0) {
        // Auto-compute: ensures erfc(alpha*cutoff) ~ 1e-6
        alpha_ = std::sqrt(12.0 * std::log(10.0)) / cutoff_;
    } else {
        alpha_ = alpha;
    }
    sqrt_alpha_ = std::sqrt(alpha_);
    sqrt_alpha_pi_ = std::sqrt(alpha_ / M_PI);
    initialized_ = false;
}

inline void PMECoulomb::set_grid(int grid_x, int grid_y, int grid_z) {
    grid_ = {grid_x, grid_y, grid_z};
    initialized_ = false;
}

inline void PMECoulomb::set_order(int order) {
    if (order < 2) {
        throw std::invalid_argument("PME order must be >= 2");
    }
    order_ = order;
    initialized_ = false;
}

inline void PMECoulomb::initialize() {
    // Allocate charge grid
    Q_.resize(grid_[0] * grid_[1] * grid_[2]);

    // Compute B-spline moduli
    compute_bspline_moduli();

    initialized_ = true;
}

inline void PMECoulomb::compute_bspline_moduli() {
    // Generate B-spline values at w=0
    std::vector<Scalar> bsp_arr(order_ + 1, 0.0);
    std::vector<Scalar> theta(order_), dtheta(order_);

    fill_bspline(0.0, theta, dtheta);

    // Store in array format
    for (int i = 0; i < order_; ++i) {
        bsp_arr[i + 1] = theta[i];
    }

    // Compute DFT modulus for each dimension
    auto compute_dft_mod = [&](int nfft, std::vector<Scalar>& bsp_mod) {
        bsp_mod.resize(nfft);

        for (int k = 0; k < nfft; ++k) {
            Scalar sum1 = 0.0, sum2 = 0.0;
            for (int j = 0; j <= order_; ++j) {
                Scalar arg = 2.0 * M_PI * k * j / nfft;
                sum1 += bsp_arr[j] * std::cos(arg);
                sum2 += bsp_arr[j] * std::sin(arg);
            }
            bsp_mod[k] = sum1 * sum1 + sum2 * sum2;
        }

        // Smooth near-zero values for numerical stability
        for (int k = 0; k < nfft; ++k) {
            if (bsp_mod[k] < 1e-7) {
                int km1 = (k - 1 + nfft) % nfft;
                int kp1 = (k + 1) % nfft;
                bsp_mod[k] = 0.5 * (bsp_mod[km1] + bsp_mod[kp1]);
            }
        }
    };

    compute_dft_mod(grid_[0], bsp_mod_x_);
    compute_dft_mod(grid_[1], bsp_mod_y_);
    compute_dft_mod(grid_[2], bsp_mod_z_);
}

inline void PMECoulomb::fill_bspline(Scalar w, std::vector<Scalar>& theta,
                                     std::vector<Scalar>& dtheta) const {
    theta.assign(order_, 0.0);
    dtheta.assign(order_, 0.0);

    // Initialize with linear basis
    theta[0] = 1.0 - w;
    theta[1] = w;

    // Cox-de Boor recursion to build higher-order splines
    for (int k = 3; k <= order_; ++k) {
        Scalar div = 1.0 / (k - 1);
        // Work backwards to avoid overwriting needed values
        theta[k - 1] = div * w * theta[k - 2];
        for (int j = k - 2; j >= 1; --j) {
            theta[j] = div * ((w + j) * theta[j - 1] + (k - j - w) * theta[j]);
        }
        theta[0] = div * (1.0 - w) * theta[0];
    }

    // Differentiation: dtheta[i] = theta[i-1] - theta[i]
    dtheta[0] = -theta[0];
    for (int i = 1; i < order_; ++i) {
        dtheta[i] = theta[i - 1] - theta[i];
    }
}

inline void PMECoulomb::compute_real_space(
    AtomicSystem& system,
    NeighborList& neighbors,
    bool compute_forces,
    bool compute_virial,
    Scalar& energy,
    Mat3& virial) {

    const std::size_t num_atoms = system.num_atoms();
    const Mat3& cell = system.cell();

    for (std::size_t i = 0; i < num_atoms; ++i) {
        Scalar qi = charges_[i];
        if (std::abs(qi) < 1e-15) continue;

        Vec3 ri = system.position(i).matrix();

        auto [begin, end] = neighbors.neighbors(i);
        for (auto it = begin; it != end; ++it) {
            const auto& neigh = *it;
            std::size_t j = neigh.index;

            Scalar qj = charges_[j];
            if (std::abs(qj) < 1e-15) continue;

            Vec3 rj = system.position(j).matrix();
            Vec3 dr = rj - ri;
            dr += cell.col(0) * neigh.cell_shift[0];
            dr += cell.col(1) * neigh.cell_shift[1];
            dr += cell.col(2) * neigh.cell_shift[2];

            Scalar r_sq = dr.squaredNorm();
            if (r_sq >= cutoff_sq_ || r_sq < 1e-20) continue;

            Scalar r = std::sqrt(r_sq);
            Scalar inv_r = 1.0 / r;

            // erfc damped energy
            Scalar erfc_val = std::erfc(sqrt_alpha_ * r);
            Scalar pair_energy = COULOMB_CONST * qi * qj * erfc_val * inv_r;

            // Half energy for full neighbor list
            energy += 0.5 * pair_energy;

            if (compute_forces || compute_virial) {
                // Force: includes both erfc and Gaussian derivative
                Scalar exp_val = std::exp(-alpha_ * r_sq);
                Scalar force_factor = COULOMB_CONST * qi * qj *
                    (erfc_val + 2.0 * sqrt_alpha_pi_ * r * exp_val) *
                    inv_r * inv_r * inv_r;

                Vec3 force = force_factor * dr;

                if (compute_forces) {
                    system.forces().col(i) -= force.array();
                }

                if (compute_virial) {
                    virial += 0.5 * dr * force.transpose();
                }
            }
        }
    }
}

inline void PMECoulomb::compute_reciprocal_space(
    AtomicSystem& system,
    bool compute_forces,
    bool compute_virial,
    Scalar& energy,
    Mat3& virial) {

    const std::size_t num_atoms = system.num_atoms();
    const Mat3& cell = system.cell();

    // Compute inverse cell (reciprocal lattice vectors / 2*pi)
    Mat3 inv_cell = cell.inverse();
    Scalar volume = std::abs(cell.determinant());

    // Allocate per-atom arrays
    fr_x_.resize(num_atoms);
    fr_y_.resize(num_atoms);
    fr_z_.resize(num_atoms);

    theta_x_.resize(num_atoms * order_);
    theta_y_.resize(num_atoms * order_);
    theta_z_.resize(num_atoms * order_);
    dtheta_x_.resize(num_atoms * order_);
    dtheta_y_.resize(num_atoms * order_);
    dtheta_z_.resize(num_atoms * order_);

    // Compute fractional coordinates and B-spline coefficients
    for (std::size_t n = 0; n < num_atoms; ++n) {
        Vec3 r = system.position(n).matrix();

        // Fractional coordinates
        Vec3 frac = inv_cell * r;

        // Wrap to [0, 1) and scale by grid
        auto wrap_frac = [](Scalar f) {
            f = f - std::floor(f);
            if (f < 0) f += 1.0;
            if (f >= 1.0) f -= 1.0;
            return f;
        };

        fr_x_[n] = wrap_frac(frac[0]) * grid_[0];
        fr_y_[n] = wrap_frac(frac[1]) * grid_[1];
        fr_z_[n] = wrap_frac(frac[2]) * grid_[2];

        // B-spline coefficients
        Scalar wx = fr_x_[n] - std::floor(fr_x_[n]);
        Scalar wy = fr_y_[n] - std::floor(fr_y_[n]);
        Scalar wz = fr_z_[n] - std::floor(fr_z_[n]);

        std::vector<Scalar> tx(order_), ty(order_), tz(order_);
        std::vector<Scalar> dtx(order_), dty(order_), dtz(order_);

        fill_bspline(wx, tx, dtx);
        fill_bspline(wy, ty, dty);
        fill_bspline(wz, tz, dtz);

        for (int k = 0; k < order_; ++k) {
            theta_x_[n * order_ + k] = tx[k];
            theta_y_[n * order_ + k] = ty[k];
            theta_z_[n * order_ + k] = tz[k];
            dtheta_x_[n * order_ + k] = dtx[k];
            dtheta_y_[n * order_ + k] = dty[k];
            dtheta_z_[n * order_ + k] = dtz[k];
        }
    }

    // Clear and spread charges onto grid
    std::fill(Q_.begin(), Q_.end(), std::complex<Scalar>(0.0, 0.0));

    for (std::size_t n = 0; n < num_atoms; ++n) {
        Scalar q = charges_[n];
        if (std::abs(q) < 1e-15) continue;

        int i0 = static_cast<int>(fr_x_[n]) - order_ + 1;
        int j0 = static_cast<int>(fr_y_[n]) - order_ + 1;
        int k0 = static_cast<int>(fr_z_[n]) - order_ + 1;

        for (int ith3 = 0; ith3 < order_; ++ith3) {
            int k = k0 + ith3;
            while (k < 0) k += grid_[2];
            while (k >= grid_[2]) k -= grid_[2];

            for (int ith2 = 0; ith2 < order_; ++ith2) {
                int j = j0 + ith2;
                while (j < 0) j += grid_[1];
                while (j >= grid_[1]) j -= grid_[1];

                Scalar prod = theta_y_[n * order_ + ith2] *
                              theta_z_[n * order_ + ith3] * q;

                for (int ith1 = 0; ith1 < order_; ++ith1) {
                    int i = i0 + ith1;
                    while (i < 0) i += grid_[0];
                    while (i >= grid_[0]) i -= grid_[0];

                    int idx = i + grid_[0] * (j + grid_[1] * k);
                    Q_[idx] += theta_x_[n * order_ + ith1] * prod;
                }
            }
        }
    }

    // Forward FFT
    fft_forward();

    // K-space energy and modify Q for field calculation
    Scalar fac = M_PI * M_PI / (alpha_);
    Scalar vol_fac = 1.0 / (M_PI * volume);

    for (int k3 = 0; k3 < grid_[2]; ++k3) {
        int m3 = (k3 <= grid_[2] / 2) ? k3 : k3 - grid_[2];

        for (int k2 = 0; k2 < grid_[1]; ++k2) {
            int m2 = (k2 <= grid_[1] / 2) ? k2 : k2 - grid_[1];

            for (int k1 = 0; k1 < grid_[0]; ++k1) {
                // Skip k=0
                if (k1 == 0 && k2 == 0 && k3 == 0) continue;

                int m1 = (k1 <= grid_[0] / 2) ? k1 : k1 - grid_[0];

                // Reciprocal lattice vector (in units of 2*pi/cell)
                Vec3 mvec;
                mvec[0] = inv_cell(0, 0) * m1 + inv_cell(0, 1) * m2 + inv_cell(0, 2) * m3;
                mvec[1] = inv_cell(1, 0) * m1 + inv_cell(1, 1) * m2 + inv_cell(1, 2) * m3;
                mvec[2] = inv_cell(2, 0) * m1 + inv_cell(2, 1) * m2 + inv_cell(2, 2) * m3;

                Scalar msq = mvec.squaredNorm();

                // B-spline modulus correction
                Scalar denom = bsp_mod_x_[k1] * bsp_mod_y_[k2] * bsp_mod_z_[k3] * msq;
                if (denom < 1e-15) continue;

                Scalar eterm = std::exp(-fac * msq) / denom;

                // Structure factor
                int idx = k1 + grid_[0] * (k2 + grid_[1] * k3);
                Scalar struc2 = std::norm(Q_[idx]);

                // Energy
                energy += COULOMB_CONST * vol_fac * eterm * struc2;

                if (compute_virial) {
                    Scalar vterm = 2.0 * (fac * msq + 1.0) / msq;
                    for (int a = 0; a < 3; ++a) {
                        for (int b = 0; b < 3; ++b) {
                            Scalar vab = vterm * mvec[a] * mvec[b];
                            if (a == b) vab -= 1.0;
                            virial(a, b) += COULOMB_CONST * vol_fac * eterm * struc2 * vab;
                        }
                    }
                }

                // Modify Q for field calculation
                Q_[idx] *= eterm;
            }
        }
    }

    // Backward FFT
    fft_backward();

    // Interpolate forces from grid
    if (compute_forces) {
        for (std::size_t n = 0; n < num_atoms; ++n) {
            Scalar q = charges_[n];
            if (std::abs(q) < 1e-15) continue;

            Scalar f1 = 0.0, f2 = 0.0, f3 = 0.0;

            int i0 = static_cast<int>(fr_x_[n]) - order_ + 1;
            int j0 = static_cast<int>(fr_y_[n]) - order_ + 1;
            int k0 = static_cast<int>(fr_z_[n]) - order_ + 1;

            for (int ith3 = 0; ith3 < order_; ++ith3) {
                int k = k0 + ith3;
                while (k < 0) k += grid_[2];
                while (k >= grid_[2]) k -= grid_[2];

                for (int ith2 = 0; ith2 < order_; ++ith2) {
                    int j = j0 + ith2;
                    while (j < 0) j += grid_[1];
                    while (j >= grid_[1]) j -= grid_[1];

                    for (int ith1 = 0; ith1 < order_; ++ith1) {
                        int i = i0 + ith1;
                        while (i < 0) i += grid_[0];
                        while (i >= grid_[0]) i -= grid_[0];

                        int idx = i + grid_[0] * (j + grid_[1] * k);
                        Scalar term = Q_[idx].real();

                        // Gradient in grid coordinates
                        f1 -= grid_[0] * term * dtheta_x_[n * order_ + ith1] *
                              theta_y_[n * order_ + ith2] * theta_z_[n * order_ + ith3];
                        f2 -= grid_[1] * term * theta_x_[n * order_ + ith1] *
                              dtheta_y_[n * order_ + ith2] * theta_z_[n * order_ + ith3];
                        f3 -= grid_[2] * term * theta_x_[n * order_ + ith1] *
                              theta_y_[n * order_ + ith2] * dtheta_z_[n * order_ + ith3];
                    }
                }
            }

            // Convert to Cartesian forces
            Scalar force_fac = COULOMB_CONST * q * vol_fac;
            Vec3 force;
            force[0] = force_fac * (inv_cell(0, 0) * f1 + inv_cell(0, 1) * f2 + inv_cell(0, 2) * f3);
            force[1] = force_fac * (inv_cell(1, 0) * f1 + inv_cell(1, 1) * f2 + inv_cell(1, 2) * f3);
            force[2] = force_fac * (inv_cell(2, 0) * f1 + inv_cell(2, 1) * f2 + inv_cell(2, 2) * f3);

            system.forces().col(n) += force.array();
        }
    }
}

inline Scalar PMECoulomb::compute_self_energy() const {
    Scalar self_energy = 0.0;
    for (const auto& q : charges_) {
        self_energy += q * q;
    }
    return COULOMB_CONST * sqrt_alpha_pi_ * self_energy;
}

// Simple 3D FFT implementation (Cooley-Tukey radix-2 for powers of 2)
// For production use, replace with FFTW3

inline void PMECoulomb::fft_forward() {
    const int nx = grid_[0], ny = grid_[1], nz = grid_[2];

    // Helper for 1D FFT
    auto fft1d = [](std::complex<Scalar>* data, int n, int stride, bool inverse) {
        // Bit-reversal permutation
        for (int i = 0, j = 0; i < n; ++i) {
            if (j > i) {
                std::swap(data[i * stride], data[j * stride]);
            }
            int m = n / 2;
            while (m >= 1 && j >= m) {
                j -= m;
                m /= 2;
            }
            j += m;
        }

        // Cooley-Tukey
        for (int mmax = 1; mmax < n; mmax *= 2) {
            Scalar theta = (inverse ? 1.0 : -1.0) * M_PI / mmax;
            std::complex<Scalar> wp(std::cos(theta) - 1.0, std::sin(theta));
            std::complex<Scalar> w(1.0, 0.0);

            for (int m = 0; m < mmax; ++m) {
                for (int i = m; i < n; i += 2 * mmax) {
                    int j = i + mmax;
                    std::complex<Scalar> temp = w * data[j * stride];
                    data[j * stride] = data[i * stride] - temp;
                    data[i * stride] += temp;
                }
                w += w * wp;
            }
        }

        if (inverse) {
            for (int i = 0; i < n; ++i) {
                data[i * stride] /= n;
            }
        }
    };

    // FFT along x
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            fft1d(&Q_[j * nx + k * nx * ny], nx, 1, false);
        }
    }

    // FFT along y
    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
            fft1d(&Q_[i + k * nx * ny], ny, nx, false);
        }
    }

    // FFT along z
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fft1d(&Q_[i + j * nx], nz, nx * ny, false);
        }
    }
}

inline void PMECoulomb::fft_backward() {
    const int nx = grid_[0], ny = grid_[1], nz = grid_[2];

    auto fft1d = [](std::complex<Scalar>* data, int n, int stride, bool inverse) {
        for (int i = 0, j = 0; i < n; ++i) {
            if (j > i) {
                std::swap(data[i * stride], data[j * stride]);
            }
            int m = n / 2;
            while (m >= 1 && j >= m) {
                j -= m;
                m /= 2;
            }
            j += m;
        }

        for (int mmax = 1; mmax < n; mmax *= 2) {
            Scalar theta = (inverse ? 1.0 : -1.0) * M_PI / mmax;
            std::complex<Scalar> wp(std::cos(theta) - 1.0, std::sin(theta));
            std::complex<Scalar> w(1.0, 0.0);

            for (int m = 0; m < mmax; ++m) {
                for (int i = m; i < n; i += 2 * mmax) {
                    int j_idx = i + mmax;
                    std::complex<Scalar> temp = w * data[j_idx * stride];
                    data[j_idx * stride] = data[i * stride] - temp;
                    data[i * stride] += temp;
                }
                w += w * wp;
            }
        }

        if (inverse) {
            for (int i = 0; i < n; ++i) {
                data[i * stride] /= n;
            }
        }
    };

    // Inverse FFT along z
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            fft1d(&Q_[i + j * nx], nz, nx * ny, true);
        }
    }

    // Inverse FFT along y
    for (int k = 0; k < nz; ++k) {
        for (int i = 0; i < nx; ++i) {
            fft1d(&Q_[i + k * nx * ny], ny, nx, true);
        }
    }

    // Inverse FFT along x
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            fft1d(&Q_[j * nx + k * nx * ny], nx, 1, true);
        }
    }
}

inline PotentialResults PMECoulomb::compute_impl(
    AtomicSystem& system,
    NeighborList& neighbors,
    bool compute_forces,
    bool compute_virial) {

    PotentialResults results;
    const std::size_t num_atoms = system.num_atoms();

    if (charges_.size() != num_atoms) {
        throw std::runtime_error("PMECoulomb: charges array size mismatch");
    }

    // Check PBC
    if (!system.pbc()[0] || !system.pbc()[1] || !system.pbc()[2]) {
        throw std::runtime_error("PMECoulomb requires 3D periodic boundary conditions");
    }

    if (!initialized_) {
        initialize();
    }

    Scalar real_energy = 0.0, recip_energy = 0.0;
    Mat3 real_virial = Mat3::Zero(), recip_virial = Mat3::Zero();

    // Real-space contribution
    compute_real_space(system, neighbors, compute_forces, compute_virial,
                       real_energy, real_virial);

    // Reciprocal space contribution
    compute_reciprocal_space(system, compute_forces, compute_virial,
                             recip_energy, recip_virial);

    // Self-energy correction
    Scalar self_energy = compute_self_energy();

    // Total
    results.energy = real_energy + recip_energy - self_energy;
    results.virial = real_virial + recip_virial;

    return results;
}

} // namespace atomistica
