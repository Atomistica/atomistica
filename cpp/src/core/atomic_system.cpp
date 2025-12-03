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

#include "atomistica/core/atomic_system.hpp"

#include <cmath>
#include <stdexcept>

namespace atomistica {

// PropertyMap implementation

bool PropertyMap::has(const std::string& name) const {
    return properties_.find(name) != properties_.end();
}

void PropertyMap::remove(const std::string& name) {
    properties_.erase(name);
}

void PropertyMap::resize(std::size_t n) {
    for (auto& [name, prop] : properties_) {
        std::visit([n](auto& arr) {
            using T = std::decay_t<decltype(arr)>;
            if constexpr (std::is_same_v<T, Array3X>) {
                arr.conservativeResize(3, n);
            } else {
                arr.conservativeResize(n);
            }
        }, prop);
    }
}

// AtomicSystem implementation

AtomicSystem::AtomicSystem(std::size_t num_atoms) {
    resize(num_atoms);
}

void AtomicSystem::resize(std::size_t num_atoms) {
    num_atoms_ = num_atoms;
    positions_.resize(3, num_atoms);
    positions_.setZero();
    atomic_numbers_.resize(num_atoms);
    atomic_numbers_.setZero();
    forces_.resize(3, num_atoms);
    forces_.setZero();
    properties_.resize(num_atoms);
}

void AtomicSystem::set_cell(const Mat3& cell) {
    cell_ = cell;
    ++cell_revision_;
}

void AtomicSystem::zero_forces() {
    forces_.setZero();
}

Scalar AtomicSystem::volume() const {
    // Volume = |det(cell)| = |a . (b x c)|
    return std::abs(cell_.determinant());
}

Mat3 AtomicSystem::inverse_cell() const {
    if (inverse_cell_revision_ != cell_revision_) {
        inverse_cell_ = cell_.inverse();
        inverse_cell_revision_ = cell_revision_;
    }
    return inverse_cell_;
}

Vec3 AtomicSystem::minimum_image(const Vec3& dr) const {
    // Convert to fractional coordinates
    Vec3 s = inverse_cell() * dr;

    // Apply minimum image convention for periodic directions
    for (int i = 0; i < 3; ++i) {
        if (pbc_[i]) {
            s(i) -= std::round(s(i));
        }
    }

    // Convert back to Cartesian
    return cell_ * s;
}

Vec3 AtomicSystem::wrap_position(const Vec3& r) const {
    // Convert to fractional coordinates
    Vec3 s = inverse_cell() * r;

    // Wrap into [0, 1) for periodic directions
    for (int i = 0; i < 3; ++i) {
        if (pbc_[i]) {
            s(i) -= std::floor(s(i));
        }
    }

    // Convert back to Cartesian
    return cell_ * s;
}

} // namespace atomistica
