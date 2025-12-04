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

/**
 * @file integrators.hpp
 * @brief Main header for molecular dynamics integrators
 *
 * This header includes all integrator components:
 * - Velocity Verlet integrator
 * - Thermostats (Berendsen, Langevin, Nose-Hoover)
 * - Barostats (Berendsen, Andersen, Parrinello-Rahman)
 */

#include "verlet.hpp"
#include "thermostats.hpp"
#include "barostats.hpp"
