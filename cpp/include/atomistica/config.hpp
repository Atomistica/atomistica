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

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace atomistica {

// Scalar type (can be changed for single precision if needed)
using Scalar = double;

// Vector and matrix types
using Vec3 = Eigen::Matrix<Scalar, 3, 1>;
using Mat3 = Eigen::Matrix<Scalar, 3, 3>;
using VecX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
using MatX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

// Array types for per-atom data (row-major for cache efficiency when iterating atoms)
using Array3X = Eigen::Array<Scalar, 3, Eigen::Dynamic>;
using ArrayX = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
using ArrayXi = Eigen::Array<int, Eigen::Dynamic, 1>;

// Constants
constexpr Scalar PI = 3.14159265358979323846;

} // namespace atomistica
