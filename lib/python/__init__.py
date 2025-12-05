# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
# and others. See the AUTHORS file in the top-level Atomistica directory.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================

"""
Atomistica C++ - Modern C++ implementation of interatomic potentials.
"""

from ._atomistica_cpp import (
    # Core classes
    AtomicSystem,
    NeighborList,
    Neighbor,
    PotentialResults,
    # Potentials
    LJCut,
    LJCutShift,
    # Math utilities
    CubicSpline,
    NonUniformSpline,
)

from .ase_calculator import Atomistica

__all__ = [
    'AtomicSystem',
    'NeighborList',
    'Neighbor',
    'PotentialResults',
    'LJCut',
    'LJCutShift',
    'CubicSpline',
    'NonUniformSpline',
    'Atomistica',
]
