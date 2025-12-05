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

// Include the new BOP kernel implementations
#include "bop_kernel.hpp"

namespace atomistica {

// Backward compatibility: BOPBase is now an alias for the appropriate kernel
// based on the Screening template parameter

/**
 * @brief Base class for Bond-Order Potentials (backward compatibility)
 *
 * This is maintained for backward compatibility. New code should use
 * BOPKernel (unscreened) or ScreenedBOPKernel (screened) directly.
 *
 * @tparam Derived The derived potential class (CRTP)
 * @tparam Screening Whether screening is enabled
 */
template<typename Derived, bool Screening = false>
using BOPBase = std::conditional_t<Screening,
                                   ScreenedBOPKernel<Derived>,
                                   BOPKernel<Derived>>;

} // namespace atomistica
