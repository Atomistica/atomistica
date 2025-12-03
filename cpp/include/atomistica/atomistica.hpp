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

// Configuration and common types
#include "config.hpp"

// Core components
#include "core/atomic_system.hpp"
#include "core/neighbor_list.hpp"

// Math utilities
#include "math/cutoff_functions.hpp"
#include "math/spline.hpp"

// Potential base
#include "potentials/potential_base.hpp"

// Pair potentials
#include "potentials/pair/lj.hpp"

// Bond-order potentials
#include "potentials/bop/bop_base.hpp"
#include "potentials/bop/tersoff.hpp"
