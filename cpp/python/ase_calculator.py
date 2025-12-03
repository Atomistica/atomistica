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
ASE Calculator interface for Atomistica C++.
"""

import numpy as np
from ase.calculators.calculator import Calculator, all_changes

from ._atomistica_cpp import AtomicSystem, NeighborList, LJCut, LJCutShift


class Atomistica(Calculator):
    """
    ASE Calculator for Atomistica C++ potentials.

    Parameters
    ----------
    potential : object
        An Atomistica potential object (e.g., LJCut, LJCutShift).
    """

    implemented_properties = ['energy', 'forces', 'stress']

    def __init__(self, potential, **kwargs):
        super().__init__(**kwargs)
        self._potential = potential
        self._system = AtomicSystem()
        self._neighbors = NeighborList()
        self._neighbors.set_cutoff(potential.cutoff())
        self._neighbors.set_verlet_shell(0.5)  # Default Verlet shell

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)

        # Update system from atoms
        n = len(atoms)
        self._system.resize(n)

        # Set cell (transpose because ASE uses row vectors, we use column vectors)
        self._system.cell = np.array(atoms.cell).T
        self._system.pbc = list(atoms.pbc)

        # Set positions (transpose: ASE is (N, 3), we use (3, N))
        self._system.positions = atoms.positions.T

        # Set atomic numbers
        self._system.atomic_numbers = atoms.numbers

        # Update neighbor list
        self._neighbors.update(self._system)

        # Zero forces
        self._system.zero_forces()

        # Compute
        compute_forces = 'forces' in properties
        compute_virial = 'stress' in properties

        results = self._potential.compute(
            self._system,
            self._neighbors,
            compute_forces,
            compute_virial
        )

        self.results['energy'] = results.energy

        if compute_forces:
            # Transpose back to ASE format (N, 3)
            self.results['forces'] = np.array(self._system.forces).T

        if compute_virial:
            # Convert virial to stress
            # stress = -virial / volume
            volume = atoms.get_volume()
            virial = np.array(results.virial)
            # Convert to Voigt notation: xx, yy, zz, yz, xz, xy
            stress = np.array([
                -virial[0, 0] / volume,
                -virial[1, 1] / volume,
                -virial[2, 2] / volume,
                -virial[1, 2] / volume,
                -virial[0, 2] / volume,
                -virial[0, 1] / volume,
            ])
            self.results['stress'] = stress
