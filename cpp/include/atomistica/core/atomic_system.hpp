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
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

#include "../config.hpp"

namespace atomistica {

/**
 * @brief Type-erased property storage for per-atom data
 *
 * Allows storing arbitrary per-atom properties (charges, velocities, etc.)
 * with type-safe access.
 */
class PropertyMap {
public:
    using PropertyVariant = std::variant<
        ArrayX,      // Scalar per-atom (e.g., charge, energy)
        Array3X,     // 3-vector per-atom (e.g., velocity, force)
        ArrayXi      // Integer per-atom (e.g., type, molecule ID)
    >;

    template<typename T>
    void add(const std::string& name, std::size_t size);

    template<typename T>
    T& get(const std::string& name);

    template<typename T>
    const T& get(const std::string& name) const;

    bool has(const std::string& name) const;
    void remove(const std::string& name);
    void resize(std::size_t n);

private:
    std::unordered_map<std::string, PropertyVariant> properties_;
};

/**
 * @brief Container for atomic system data
 *
 * Stores positions, atomic numbers, cell vectors, and periodic boundary conditions.
 * Also provides extensible per-atom property storage.
 *
 * This replaces the Fortran particles_t type.
 */
class AtomicSystem {
public:
    AtomicSystem() = default;
    explicit AtomicSystem(std::size_t num_atoms);

    // Resize the system
    void resize(std::size_t num_atoms);

    // Number of atoms
    std::size_t num_atoms() const { return num_atoms_; }

    // Cell vectors (columns are the lattice vectors a, b, c)
    Mat3& cell() { return cell_; }
    const Mat3& cell() const { return cell_; }
    void set_cell(const Mat3& cell);

    // Periodic boundary conditions
    std::array<bool, 3>& pbc() { return pbc_; }
    const std::array<bool, 3>& pbc() const { return pbc_; }

    // Positions (3 x N array, column i is position of atom i)
    Array3X& positions() { return positions_; }
    const Array3X& positions() const { return positions_; }

    // Single atom position access
    auto position(std::size_t i) { return positions_.col(i); }
    auto position(std::size_t i) const { return positions_.col(i); }

    // Atomic numbers
    ArrayXi& atomic_numbers() { return atomic_numbers_; }
    const ArrayXi& atomic_numbers() const { return atomic_numbers_; }

    // Forces (accumulated by potentials)
    Array3X& forces() { return forces_; }
    const Array3X& forces() const { return forces_; }

    // Zero forces
    void zero_forces();

    // Extensible property storage
    PropertyMap& properties() { return properties_; }
    const PropertyMap& properties() const { return properties_; }

    // Cell operations
    Scalar volume() const;
    Mat3 inverse_cell() const;

    // Minimum image convention for distance vectors
    Vec3 minimum_image(const Vec3& dr) const;

    // Wrap position into cell
    Vec3 wrap_position(const Vec3& r) const;

    // Change tracking for lazy updates (e.g., neighbor lists)
    int position_revision() const { return position_revision_; }
    int cell_revision() const { return cell_revision_; }
    void positions_changed() { ++position_revision_; }
    void cell_changed() { ++cell_revision_; }

private:
    std::size_t num_atoms_ = 0;

    Mat3 cell_ = Mat3::Identity();
    std::array<bool, 3> pbc_ = {true, true, true};

    Array3X positions_;
    ArrayXi atomic_numbers_;
    Array3X forces_;

    PropertyMap properties_;

    int position_revision_ = 0;
    int cell_revision_ = 0;

    // Cached inverse cell (lazily computed)
    mutable Mat3 inverse_cell_;
    mutable int inverse_cell_revision_ = -1;
};

// Template implementations

template<typename T>
void PropertyMap::add(const std::string& name, std::size_t size) {
    if constexpr (std::is_same_v<T, ArrayX>) {
        ArrayX arr = ArrayX::Zero(size);
        properties_[name] = std::move(arr);
    } else if constexpr (std::is_same_v<T, Array3X>) {
        Array3X arr = Array3X::Zero(3, size);
        properties_[name] = std::move(arr);
    } else if constexpr (std::is_same_v<T, ArrayXi>) {
        ArrayXi arr = ArrayXi::Zero(size);
        properties_[name] = std::move(arr);
    } else {
        static_assert(sizeof(T) == 0, "Unsupported property type");
    }
}

template<typename T>
T& PropertyMap::get(const std::string& name) {
    return std::get<T>(properties_.at(name));
}

template<typename T>
const T& PropertyMap::get(const std::string& name) const {
    return std::get<T>(properties_.at(name));
}

} // namespace atomistica
