"""Native MDCore interface.
"""

import numpy as np

from ase.data import atomic_numbers

from _atomistica import *

###

_lees_edwards_info_str = 'shear_dx'

###

def from_atoms(atoms):
    pbc = np.array( [ True, True, True ] )
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

