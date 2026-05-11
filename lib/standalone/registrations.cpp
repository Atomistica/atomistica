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

// Single compilation unit that includes every hook/integrator header and
// registers all classes with the Registry.  Headers that contain their own
// REGISTER_* call handle registration when included here; those that don't
// are registered explicitly below.

#include "registry.hpp"

// ---- Potentials (self-register via REGISTER_POTENTIAL / REGISTER_COULOMB) --
#include "potentials/bop_wrappers.hpp"
#include "potentials/dftb_wrapper.hpp"
#include "potentials/eam_wrapper.hpp"
#include "potentials/pair_wrappers.hpp"
#include "potentials/coulomb_wrappers.hpp"

// ---- Integrators ----------------------------------------------------------
#include "integrators/velocity_verlet.hpp"
#include "integrators/fire.hpp"
#include "integrators/no_integration.hpp"
#include "integrators/andersen_p.hpp"

// ---- Thermostats & barostats -----------------------------------------------
#include "hooks/thermostats/berendsen_t.hpp"
#include "hooks/thermostats/langevin_t.hpp"
#include "hooks/thermostats/peters_t.hpp"
#include "hooks/thermostats/remove_rotation.hpp"
#include "hooks/thermostats/berendsen_p.hpp"

// ---- Output ----------------------------------------------------------------
// Each of these headers contains its own REGISTER_HOOK call.
#include "hooks/output/output_energy.hpp"
#include "hooks/output/output_xyz.hpp"
#include "hooks/output/output_cfg.hpp"
#include "hooks/output/output_nc.hpp"

// ---- Analysis --------------------------------------------------------------
#include "hooks/analysis/diffusion_coefficient.hpp"
#include "hooks/analysis/slicing.hpp"
#include "hooks/analysis/heat_flux.hpp"

// ---- Deformation -----------------------------------------------------------
#include "hooks/deformation/constant_strain_rate.hpp"
#include "hooks/deformation/constant_velocity.hpp"
#include "hooks/deformation/constant_force.hpp"
#include "hooks/deformation/harmonic_hook.hpp"
#include "hooks/deformation/confinement.hpp"

// ---- Constraints -----------------------------------------------------------
// settle.hpp also registers FreezeAtoms and RATTLE.
#include "hooks/constraints/settle.hpp"

// ---- Charge equilibration --------------------------------------------------
#include "hooks/charges/charge_equilibration.hpp"

// ---------------------------------------------------------------------------
// Explicit registrations for classes whose headers do not self-register.
// ---------------------------------------------------------------------------

namespace {

// Integrators
const bool _reg_VelocityVerlet =
    atomistica::Registry::instance().register_integrator(
        "VelocityVerlet",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::VelocityVerlet>(cfg);
        });

const bool _reg_FIRE =
    atomistica::Registry::instance().register_integrator(
        "FIRE",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::FIRE>(cfg);
        });

const bool _reg_NoIntegration =
    atomistica::Registry::instance().register_integrator(
        "NoIntegration",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::NoIntegration>(cfg);
        });

const bool _reg_AndersenP =
    atomistica::Registry::instance().register_integrator(
        "AndersenP",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::AndersenP>(cfg);
        });

// Thermostats & barostats (headers have no self-registration)
const bool _reg_BerendsenT =
    atomistica::Registry::instance().register_hook(
        "BerendsenT",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::BerendsenT>(cfg);
        });

const bool _reg_LangevinT =
    atomistica::Registry::instance().register_hook(
        "LangevinT",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::LangevinT>(cfg);
        });

const bool _reg_PetersT =
    atomistica::Registry::instance().register_hook(
        "PetersT",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::PetersT>(cfg);
        });

const bool _reg_RemoveRotation =
    atomistica::Registry::instance().register_hook(
        "RemoveRotation",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::RemoveRotation>(cfg);
        });

const bool _reg_BerendsenP =
    atomistica::Registry::instance().register_hook(
        "BerendsenP",
        [](const atomistica::Config& cfg) {
            return std::make_unique<atomistica::BerendsenP>(cfg);
        });

} // anonymous namespace
