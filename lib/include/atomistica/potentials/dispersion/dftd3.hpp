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

#include <stdexcept>
#include <vector>

#include "../../config.hpp"
#include "../../core/atomic_system.hpp"
#include "../../core/neighbor_list.hpp"
#include "../potential_base.hpp"

#ifdef HAVE_SDFTD3
#include <s-dftd3.h>
#endif

namespace atomistica {

// Unit conversion constants
static constexpr Scalar DFTD3_BOHR    = 0.529177210903;    // Bohr in Angstrom
static constexpr Scalar DFTD3_HARTREE = 27.211386245988;   // Hartree in eV

/**
 * @brief DFT-D3 dispersion correction potential
 *
 * Grimme's DFT-D3 dispersion potential with Becke-Johnson (BJ) or
 * zero damping. Wraps the s-dftd3 library when available.
 *
 * References:
 *   S. Grimme, J. Antony, S. Ehrlich, H. Krieg,
 *   J. Chem. Phys. 132, 154104 (2010)
 *   S. Grimme, S. Ehrlich, L. Goerigk,
 *   J. Comp. Chem. 32, 1456 (2011)
 *
 * Requires the s-dftd3 library (https://github.com/awvwgk/simple-dftd3).
 * Without it the class is compiled in but compute() raises a runtime error.
 *
 * Default parameters correspond to PBE0-D3(BJ).
 */
class DFTD3Disp : public PotentialBase<DFTD3Disp> {
public:
    // Becke-Johnson damping parameters (default: PBE0-D3(BJ))
    Scalar s6    = 1.0000;
    Scalar s8    = 0.5883;
    Scalar a1    = 0.5719;
    Scalar a2    = 3.6017;

    // Zero-damping parameters
    Scalar sr6   = 0.7461;
    Scalar sr8   = 1.0000;
    Scalar alpha6 = 14.000;

    // Cutoffs in Angstrom
    Scalar cutoff_radius = 80.0;
    Scalar cutoff_cn     = 40.0;

    bool use_bj_damping = true;   // true = BJ damping, false = zero damping

    DFTD3Disp() = default;

    ~DFTD3Disp() {
#ifdef HAVE_SDFTD3
        free_resources();
#endif
    }

    // Non-copyable due to opaque pointer ownership
    DFTD3Disp(const DFTD3Disp&) = delete;
    DFTD3Disp& operator=(const DFTD3Disp&) = delete;
    DFTD3Disp(DFTD3Disp&&) = default;

    Scalar cutoff_impl() const { return cutoff_radius; }

    void bind_to_impl(AtomicSystem& system, NeighborList& neighbors) {
#ifndef HAVE_SDFTD3
        throw std::runtime_error(
            "DFTD3Disp: built without s-dftd3 support. "
            "Install s-dftd3 and rebuild.");
#else
        neighbors.set_cutoff(cutoff_radius);
        setup_calculator(system);
#endif
    }

    PotentialResults compute_impl(AtomicSystem& system,
                                  NeighborList& neighbors,
                                  bool compute_forces,
                                  bool compute_virial) {
#ifndef HAVE_SDFTD3
        throw std::runtime_error(
            "DFTD3Disp: built without s-dftd3 support. "
            "Install s-dftd3 and rebuild.");
        return {};
#else
        PotentialResults results;
        const std::size_t n = system.num_atoms();

        // Convert positions to Bohr, row-major N×3 for s-dftd3
        std::vector<double> positions_bohr(3 * n);
        for (std::size_t i = 0; i < n; ++i) {
            positions_bohr[3*i + 0] = system.positions()(0, i) / DFTD3_BOHR;
            positions_bohr[3*i + 1] = system.positions()(1, i) / DFTD3_BOHR;
            positions_bohr[3*i + 2] = system.positions()(2, i) / DFTD3_BOHR;
        }

        // Convert lattice to Bohr, row-major 3×3 for s-dftd3
        // AtomicSystem stores cell as columns (each column = lattice vector)
        // s-dftd3 wants rows
        double lattice_bohr[9];
        const Mat3& cell = system.cell();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                lattice_bohr[3*i + j] = cell(j, i) / DFTD3_BOHR;

        // PBC flags
        bool periodic[3] = {system.pbc()[0], system.pbc()[1], system.pbc()[2]};

        // Update structure
        dftd3_update_structure(error_, structure_,
                               positions_bohr.data(), lattice_bohr);
        check_error("dftd3_update_structure");

        // Output arrays
        double energy = 0.0;
        std::vector<double> gradient(3 * n, 0.0);
        double sigma[9] = {};

        dftd3_get_dispersion(error_, structure_, model_, param_,
                             &energy,
                             compute_forces ? gradient.data() : nullptr,
                             compute_virial ? sigma : nullptr);
        check_error("dftd3_get_dispersion");

        // Convert energy: Hartree -> eV
        results.energy = energy * DFTD3_HARTREE;

        if (compute_forces) {
            // Convert forces: Hartree/Bohr -> eV/Å; gradient -> forces (sign flip)
            const Scalar fac = -DFTD3_HARTREE / DFTD3_BOHR;
            for (std::size_t i = 0; i < n; ++i) {
                system.forces()(0, i) += gradient[3*i + 0] * fac;
                system.forces()(1, i) += gradient[3*i + 1] * fac;
                system.forces()(2, i) += gradient[3*i + 2] * fac;
            }
        }

        if (compute_virial) {
            // sigma is the stress tensor in Hartree; virial = -sigma * volume
            // s-dftd3 sigma is row-major 3×3
            // Convert: Hartree -> eV (no Bohr factor since it's already volumetric
            // relative to the Bohr^3 volume, but we provide Bohr lattice so sigma
            // comes out in Hartree/Bohr^3 * volume_bohr^3 = Hartree)
            const Scalar vfac = DFTD3_HARTREE;
            Mat3 virial = Mat3::Zero();
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    virial(i, j) = -sigma[3*i + j] * vfac;
            results.virial = virial;
        }

        return results;
#endif
    }

private:
#ifdef HAVE_SDFTD3
    dftd3_error     error_     = nullptr;
    dftd3_structure structure_ = nullptr;
    dftd3_model     model_     = nullptr;
    dftd3_param     param_     = nullptr;

    void free_resources() {
        if (param_)     { dftd3_delete_param(&param_);         param_    = nullptr; }
        if (model_)     { dftd3_delete_model(&model_);         model_    = nullptr; }
        if (structure_) { dftd3_delete_structure(&structure_); structure_= nullptr; }
        if (error_)     { dftd3_delete_error(&error_);         error_    = nullptr; }
    }

    void check_error(const char* where) {
        // s-dftd3 errors are signalled by a non-null message
        char* msg = dftd3_get_error_message(error_);
        if (msg) {
            std::string s(msg);
            dftd3_delete_error_message(&msg);
            throw std::runtime_error(std::string("DFTD3Disp::") + where + ": " + s);
        }
    }

    void setup_calculator(AtomicSystem& system) {
        free_resources();

        const std::size_t n = system.num_atoms();

        error_ = dftd3_new_error();

        // Atomic numbers
        std::vector<int> numbers(n);
        for (std::size_t i = 0; i < n; ++i)
            numbers[i] = static_cast<int>(system.atomic_numbers()(i));

        // Positions in Bohr (row-major N×3)
        std::vector<double> positions_bohr(3 * n);
        for (std::size_t i = 0; i < n; ++i) {
            positions_bohr[3*i + 0] = system.positions()(0, i) / DFTD3_BOHR;
            positions_bohr[3*i + 1] = system.positions()(1, i) / DFTD3_BOHR;
            positions_bohr[3*i + 2] = system.positions()(2, i) / DFTD3_BOHR;
        }

        // Lattice in Bohr (row-major 3×3)
        double lattice_bohr[9];
        const Mat3& cell = system.cell();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                lattice_bohr[3*i + j] = cell(j, i) / DFTD3_BOHR;

        bool periodic[3] = {system.pbc()[0], system.pbc()[1], system.pbc()[2]};

        structure_ = dftd3_new_structure(error_,
                                         static_cast<int>(n),
                                         numbers.data(),
                                         positions_bohr.data(),
                                         lattice_bohr,
                                         periodic);
        check_error("new_structure");

        model_ = dftd3_new_d3_model(error_, structure_);
        check_error("new_d3_model");

        if (use_bj_damping) {
            param_ = dftd3_new_bj_damping(error_, s6, s8, a1, a2);
        } else {
            param_ = dftd3_new_zero_damping(error_, s6, s8, sr6, sr8, alpha6);
        }
        check_error("new_damping");
    }
#endif
};

} // namespace atomistica
