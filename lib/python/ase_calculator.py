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
ASE Calculator interface for Atomistica C++ potentials.
"""

import numpy as np

from ._atomistica_cpp import AtomicSystem, NeighborList

try:
    from ase.calculators.calculator import Calculator, all_changes
    _ASE_AVAILABLE = True
except ImportError:
    _ASE_AVAILABLE = False
    Calculator = object
    all_changes = []


# Coulomb potentials have cutoff() == 0 (or very large) and need external charges.
_COULOMB_TYPES = ('DirectCoulomb', 'CutoffCoulomb', 'WolfCoulomb',
                  'PMECoulomb', 'FMMCoulomb')

# DFTB requires init() before compute().
_DFTB_TYPES = ('DFTB',)


def _is_coulomb(potential):
    return type(potential).__name__ in _COULOMB_TYPES


def _is_dftb(potential):
    return type(potential).__name__ in _DFTB_TYPES


def _make_potential(potential_or_class, param_name=None):
    """
    Instantiate a potential.

    Parameters
    ----------
    potential_or_class :
        Either an already-constructed potential object, or a potential
        class to be instantiated and (optionally) loaded with parameters.
    param_name : str, optional
        Name string passed to ``potential.load_parameters()``.

    Returns
    -------
    potential object
    """
    if isinstance(potential_or_class, type):
        potential = potential_or_class()
    else:
        potential = potential_or_class

    if param_name is not None:
        potential.load_parameters(param_name)

    return potential


class Atomistica(Calculator):
    """
    ASE Calculator for Atomistica C++ potentials.

    Parameters
    ----------
    potential : potential object or potential class
        An instantiated Atomistica potential, or a potential class (in which
        case ``param_name`` must be provided to load parameters).
    param_name : str, optional
        Parameter set name passed to ``potential.load_parameters()``.
        Required when ``potential`` is a class.
    verlet_shell : float, optional
        Verlet shell radius (Å) added to the cutoff for lazy neighbor-list
        updates. Default 0.5.
    charges : array_like, optional
        Per-atom charges (electrons) for Coulomb potentials. Can also be set
        later via ``atoms.arrays['charges']``.
    **kwargs :
        Passed to the ASE ``Calculator`` base class.

    Examples
    --------
    Instantiate with a pre-loaded potential::

        from atomistica_cpp import Tersoff
        pot = Tersoff()
        pot.load_parameters("Tersoff_PRB_39_5566_Si_C")
        calc = Atomistica(pot)

    Instantiate using class + name constant::

        from atomistica_cpp import Tersoff, Tersoff_PRB_39_5566_Si_C
        calc = Atomistica(Tersoff, Tersoff_PRB_39_5566_Si_C)

    DFTB::

        from atomistica_cpp import DFTB
        dftb = DFTB(skf_path='/path/to/skf', enable_scc=True)
        calc = Atomistica(dftb)
    """

    implemented_properties = ['energy', 'free_energy', 'forces', 'stress']

    def __init__(self, potential, param_name=None, verlet_shell=0.5,
                 charges=None, **kwargs):
        super().__init__(**kwargs)
        self._potential = _make_potential(potential, param_name)
        self._verlet_shell = verlet_shell
        self._charges = charges

        self._system = AtomicSystem()
        self._neighbors = NeighborList()

        cutoff = self._potential.cutoff()
        if cutoff > 0:
            self._neighbors.set_cutoff(cutoff)
            self._neighbors.set_verlet_shell(verlet_shell)

        self._dftb_initialized = False

    def set_charges(self, charges):
        """Set per-atom charges for Coulomb potentials."""
        self._charges = np.asarray(charges, dtype=float)
        import numpy as np
        self._potential.set_charges(np.asarray(charges, dtype=float))

    def _update_system(self, atoms):
        """Synchronise AtomicSystem from an ASE Atoms object."""
        n = len(atoms)
        self._system.resize(n)
        self._system.cell = np.array(atoms.cell).T
        self._system.pbc = list(atoms.pbc)
        self._system.positions = atoms.positions.T
        self._system.atomic_numbers = atoms.numbers

    def _update_neighbors(self):
        self._neighbors.update(self._system)

    def _handle_coulomb_charges(self, atoms):
        """Push charges into a Coulomb potential before computing."""
        if self._charges is not None:
            charges = np.asarray(self._charges, dtype=float)
        elif 'charges' in atoms.arrays:
            charges = atoms.get_array('charges').astype(float)
        elif hasattr(atoms, 'get_charges'):
            try:
                charges = atoms.get_charges().astype(float)
            except Exception:
                charges = np.zeros(len(atoms))
        else:
            charges = np.zeros(len(atoms))
        self._potential.set_charges(charges)

    def calculate(self, atoms=None, properties=None, system_changes=all_changes):
        if properties is None:
            properties = self.implemented_properties

        super().calculate(atoms, properties, system_changes)

        self._update_system(atoms)
        self._system.zero_forces()

        # Update neighbor list if the cutoff is finite.
        cutoff = self._potential.cutoff()
        if cutoff > 0:
            self._update_neighbors()
        else:
            # For potentials like DirectCoulomb that have no finite cutoff,
            # use a very large cutoff (whole cell) or no neighbor list.
            # The potential must handle this itself.
            self._neighbors.set_cutoff(1e10)
            self._update_neighbors()

        # Coulomb: push charges
        if _is_coulomb(self._potential):
            self._handle_coulomb_charges(atoms)

        # DFTB: call init() on first use
        if _is_dftb(self._potential) and not self._dftb_initialized:
            self._potential.init(self._system)
            self._dftb_initialized = True

        compute_forces = 'forces' in properties
        compute_virial = 'stress' in properties

        results = self._potential.compute(
            self._system,
            self._neighbors,
            compute_forces,
            compute_virial,
        )

        self.results['energy'] = float(results.energy)
        self.results['free_energy'] = float(results.energy)

        if compute_forces:
            self.results['forces'] = np.array(self._system.forces).T.copy()

        if compute_virial:
            volume = atoms.get_volume()
            virial = np.array(results.virial)
            # Voigt notation: xx, yy, zz, yz, xz, xy
            self.results['stress'] = np.array([
                -virial[0, 0] / volume,
                -virial[1, 1] / volume,
                -virial[2, 2] / volume,
                -virial[1, 2] / volume,
                -virial[0, 2] / volume,
                -virial[0, 1] / volume,
            ])
