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

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "../include/atomistica/potentials/potential_base.hpp"

#include "config.hpp"
#include "coulomb_solver.hpp"
#include "hook.hpp"
#include "integrator.hpp"

namespace atomistica {

// Runtime class registry — maps section names from md.dat to factory functions.
//
// Each class (VelocityVerlet, BerendsenT, Tersoff, …) self-registers at
// static-initialisation time using one of the REGISTER_* macros at the
// bottom of this file. mdcore.cpp then iterates the top-level sections of
// md.dat and calls the appropriate make_*() for each.
//
// There are four separate factory maps, one per class tier:
//   Potential   — charge-unaware force models (BOPs, EAM, pairs)
//   CoulombSolver — electrostatic solvers (DirectCoulomb, PME, Wolf, …)
//   Integrator  — time-steppers (VelocityVerlet, FIRE, NoIntegration)
//   Hook        — everything called per-step (thermostats, barostats,
//                  output, analysis, constraints, deformation)
//
// Names are stored and matched case-insensitively.
class Registry {
public:
    static Registry& instance();

    using PotentialFactory  = std::function<std::unique_ptr<Potential>   (const Config&)>;
    using CoulombFactory    = std::function<std::unique_ptr<CoulombSolver>(const Config&)>;
    using IntegratorFactory = std::function<std::unique_ptr<Integrator>(const Config&)>;
    using HookFactory       = std::function<std::unique_ptr<Hook>        (const Config&)>;

    // Registration — returns true so the result can drive a static bool init.
    bool register_potential  (const std::string& name, PotentialFactory  f);
    bool register_coulomb    (const std::string& name, CoulombFactory    f);
    bool register_integrator (const std::string& name, IntegratorFactory f);
    bool register_hook       (const std::string& name, HookFactory       f);

    // Query — all comparisons are case-insensitive.
    bool is_potential  (const std::string& name) const;
    bool is_coulomb    (const std::string& name) const;
    bool is_integrator (const std::string& name) const;
    bool is_hook       (const std::string& name) const;
    bool is_known      (const std::string& name) const;

    // Factory — throws std::runtime_error with a helpful message on unknown name.
    std::unique_ptr<Potential>    make_potential  (const std::string& name, const Config& cfg) const;
    std::unique_ptr<CoulombSolver>make_coulomb    (const std::string& name, const Config& cfg) const;
    std::unique_ptr<Integrator>   make_integrator (const std::string& name, const Config& cfg) const;
    std::unique_ptr<Hook>         make_hook       (const std::string& name, const Config& cfg) const;

    // Comma-separated list of all registered names across all tiers.
    // Used in error messages for unknown sections.
    std::string known_names() const;

    // Separate lists per tier (for --help style output).
    std::vector<std::string> potential_names()   const;
    std::vector<std::string> coulomb_names()     const;
    std::vector<std::string> integrator_names()  const;
    std::vector<std::string> hook_names()        const;

private:
    Registry() = default;

    static std::string to_lower(const std::string& s);

    std::map<std::string, PotentialFactory>  potentials_;
    std::map<std::string, CoulombFactory>    coulombs_;
    std::map<std::string, IntegratorFactory> integrators_;
    std::map<std::string, HookFactory>       hooks_;
};

// ---------------------------------------------------------------------------
// Self-registration macros
//
// Place one of these at namespace scope in a .cpp file (or, for header-only
// classes, in an anonymous namespace in the .hpp) to register the class.
//
// The lambda captures nothing; Class must be constructible from const Config&.
//
// Example:
//   REGISTER_POTENTIAL("Tersoff", TersoffPotential)
//   REGISTER_HOOK("BerendsenT", BerendsenT)
//   REGISTER_INTEGRATOR("VelocityVerlet", VelocityVerlet)
//   REGISTER_COULOMB("DirectCoulomb", DirectCoulombSolver)
// ---------------------------------------------------------------------------

#define REGISTER_POTENTIAL(Name, Class)                                      \
    namespace {                                                               \
    const bool _atomistica_reg_potential_##Class =                           \
        ::atomistica::Registry::instance().register_potential(               \
            Name,                                                             \
            [](const ::atomistica::Config& cfg) {                            \
                return std::make_unique<Class>(cfg);                         \
            });                                                               \
    }

#define REGISTER_COULOMB(Name, Class)                                        \
    namespace {                                                               \
    const bool _atomistica_reg_coulomb_##Class =                             \
        ::atomistica::Registry::instance().register_coulomb(                 \
            Name,                                                             \
            [](const ::atomistica::Config& cfg) {                            \
                return std::make_unique<Class>(cfg);                         \
            });                                                               \
    }

#define REGISTER_INTEGRATOR(Name, Class)                                     \
    namespace {                                                               \
    const bool _atomistica_reg_integrator_##Class =                          \
        ::atomistica::Registry::instance().register_integrator(              \
            Name,                                                             \
            [](const ::atomistica::Config& cfg) {                            \
                return std::make_unique<Class>(cfg);                         \
            });                                                               \
    }

#define REGISTER_HOOK(Name, Class)                                           \
    namespace {                                                               \
    const bool _atomistica_reg_hook_##Class =                                \
        ::atomistica::Registry::instance().register_hook(                    \
            Name,                                                             \
            [](const ::atomistica::Config& cfg) {                            \
                return std::make_unique<Class>(cfg);                         \
            });                                                               \
    }

} // namespace atomistica
