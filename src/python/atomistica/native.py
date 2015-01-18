# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
# See the AUTHORS file in the top-level Atomistica directory.
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
Native MDCore interface.
"""

import numpy as np

from ase.data import atomic_numbers

from _atomistica import *

###

_lees_edwards_info_str = 'shear_dx'

###

def from_atoms(atoms):
    pbc = np.array(atoms.get_pbc())
    particles = Particles()
    particles.allocate(len(atoms))
    particles.set_cell(atoms.get_cell(), pbc)

    Z  = particles.Z
    for i, at in enumerate(atoms):
        Z[i]   = atomic_numbers[at.symbol]

    particles.coordinates[:, :]  = atoms.get_positions()[:, :]
    if _lees_edwards_info_str in atoms.info:
        particles.set_lees_edwards(atoms.info[_lees_edwards_info_str])
    
    # Notify the Particles object of a change
    particles.I_changed_positions()

    particles.update_elements()

    return particles


def neighbor_list(particles, cutoff, avgn=100):
    neighbors = Neighbors(avgn)
    neighbors.request_interaction_range(cutoff)
    neighbors.update(particles)

    return neighbors

