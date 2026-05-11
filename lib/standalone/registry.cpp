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

#include <cctype>
#include <stdexcept>

#include "registry.hpp"

namespace atomistica {

// ---------------------------------------------------------------------------
// Singleton
// ---------------------------------------------------------------------------

Registry& Registry::instance() {
    static Registry inst;
    return inst;
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

std::string Registry::to_lower(const std::string& s) {
    std::string r = s;
    for (auto& c : r)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return r;
}

// ---------------------------------------------------------------------------
// Registration
// ---------------------------------------------------------------------------

bool Registry::register_potential(const std::string& name, PotentialFactory f) {
    potentials_[to_lower(name)] = std::move(f);
    return true;
}

bool Registry::register_coulomb(const std::string& name, CoulombFactory f) {
    coulombs_[to_lower(name)] = std::move(f);
    return true;
}

bool Registry::register_integrator(const std::string& name, IntegratorFactory f) {
    integrators_[to_lower(name)] = std::move(f);
    return true;
}

bool Registry::register_hook(const std::string& name, HookFactory f) {
    hooks_[to_lower(name)] = std::move(f);
    return true;
}

// ---------------------------------------------------------------------------
// Query
// ---------------------------------------------------------------------------

bool Registry::is_potential(const std::string& name) const {
    return potentials_.count(to_lower(name)) > 0;
}

bool Registry::is_coulomb(const std::string& name) const {
    return coulombs_.count(to_lower(name)) > 0;
}

bool Registry::is_integrator(const std::string& name) const {
    return integrators_.count(to_lower(name)) > 0;
}

bool Registry::is_hook(const std::string& name) const {
    return hooks_.count(to_lower(name)) > 0;
}

bool Registry::is_known(const std::string& name) const {
    return is_potential(name) || is_coulomb(name)
        || is_integrator(name) || is_hook(name);
}

// ---------------------------------------------------------------------------
// Factory
// ---------------------------------------------------------------------------

std::unique_ptr<Potential>
Registry::make_potential(const std::string& name, const Config& cfg) const {
    auto it = potentials_.find(to_lower(name));
    if (it == potentials_.end())
        throw std::runtime_error(
            "Unknown potential '" + name + "'. Known potentials: "
            + [this]{ std::string s; for (auto& kv : potentials_) { if (!s.empty()) s += ", "; s += kv.first; } return s; }()
        );
    return it->second(cfg);
}

std::unique_ptr<CoulombSolver>
Registry::make_coulomb(const std::string& name, const Config& cfg) const {
    auto it = coulombs_.find(to_lower(name));
    if (it == coulombs_.end())
        throw std::runtime_error(
            "Unknown Coulomb solver '" + name + "'. Known solvers: "
            + [this]{ std::string s; for (auto& kv : coulombs_) { if (!s.empty()) s += ", "; s += kv.first; } return s; }()
        );
    return it->second(cfg);
}

std::unique_ptr<Integrator>
Registry::make_integrator(const std::string& name, const Config& cfg) const {
    auto it = integrators_.find(to_lower(name));
    if (it == integrators_.end())
        throw std::runtime_error(
            "Unknown integrator '" + name + "'. Known integrators: "
            + [this]{ std::string s; for (auto& kv : integrators_) { if (!s.empty()) s += ", "; s += kv.first; } return s; }()
        );
    return it->second(cfg);
}

std::unique_ptr<Hook>
Registry::make_hook(const std::string& name, const Config& cfg) const {
    auto it = hooks_.find(to_lower(name));
    if (it == hooks_.end())
        throw std::runtime_error(
            "Unknown hook '" + name + "'. Known hooks: "
            + [this]{ std::string s; for (auto& kv : hooks_) { if (!s.empty()) s += ", "; s += kv.first; } return s; }()
        );
    return it->second(cfg);
}

// ---------------------------------------------------------------------------
// Enumeration
// ---------------------------------------------------------------------------

std::string Registry::known_names() const {
    std::string s;
    auto append_map = [&](const auto& m) {
        for (const auto& kv : m) {
            if (!s.empty()) s += ", ";
            s += kv.first;
        }
    };
    append_map(integrators_);
    append_map(potentials_);
    append_map(coulombs_);
    append_map(hooks_);
    return s;
}

std::vector<std::string> Registry::potential_names() const {
    std::vector<std::string> v;
    for (const auto& kv : potentials_) v.push_back(kv.first);
    return v;
}

std::vector<std::string> Registry::coulomb_names() const {
    std::vector<std::string> v;
    for (const auto& kv : coulombs_) v.push_back(kv.first);
    return v;
}

std::vector<std::string> Registry::integrator_names() const {
    std::vector<std::string> v;
    for (const auto& kv : integrators_) v.push_back(kv.first);
    return v;
}

std::vector<std::string> Registry::hook_names() const {
    std::vector<std::string> v;
    for (const auto& kv : hooks_) v.push_back(kv.first);
    return v;
}

} // namespace atomistica
