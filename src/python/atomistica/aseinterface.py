# ======================================================================
# Atomistica - Interatomic potential library
# https://github.com/pastewka/atomistica
# Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
# See the AUTHORS file in the top-level Atomistica directory.
#
# Copyright (2005-2013) Fraunhofer IWM
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
ASE interface to Atomistica.
"""

import copy
import inspect
from math import sqrt, log

import _atomistica

import numpy as np

from ase.data import atomic_numbers
from ase.units import Hartree, Bohr

from atomistica.parameters import *

###

_warned_about_lees_edwards = False
class Atomistica(object):
    """
    Atomistica ASE calculator.
    """

    CELL_TOL = 1e-16
    POSITIONS_TOL = 1e-16

    potential_class = None
    avgn = 100
    
    _lees_edwards_info_str = 'shear_dx'

    def __init__(self, potentials=None, avgn=None, **kwargs):
        """
        Initialize a potential. *potential* is the a native Atomistica
        potential class, in which case the arguments are given by the
        keywords of this function. Alternatively, it can be a list of tuples
        [ ( pot1, args1 ), ... ] in which case *pot1* is the name of the
        first potential and *args1* a dictionary containing its arguments.
        """

        # List of potential objects
        self.pots = [ ]
        # List of Coulomb solvers
        self.couls = [ ]

        # Loop over all potentials and check whether it is a potntial or a
        # Coulomb solver. Create two separate lists.
        if potentials is not None:
            for pot in potentials:
                if hasattr(pot, 'potential'):
                    self.couls += [ pot ]
                else:
                    self.pots += [ pot ]
        else:
            self.pots = [ self.potential_class(**kwargs) ]
            self.couls = [ ]

        if avgn:
            self.avgn = avgn

        self.particles = None
        self.nl = None

        self.q = None

        self.force_update = True

        self.kwargs = kwargs

        self.lees_edwards_dx = None
        self.lees_edwards_dv = None

        self.compute_epot_per_at = False
        self.compute_epot_per_bond = False
        self.compute_f_per_bond = False
        self.compute_wpot_per_at = False
        self.compute_wpot_per_bond = False

        self.epot_per_at = None
        self.epot_per_bond = None
        self.f_per_bond = None
        self.wpot_per_at = None
        self.wpot_per_bond = None
        

    def initialize(self, atoms):
        # For now, length and elements are fixed
        # FIXME! Override pbc since code exits when atom moves out of the box!
        pbc = np.array( [ True, True, True ] )
        self.particles  = _atomistica.Particles()
        if len(self.couls) > 0:
            self.q = np.zeros(len(atoms))
            self.phi = np.zeros(len(atoms))
            self.E = np.zeros([3,len(atoms)])
            # Coulomb callback should be directed to this wrapper object
            for pot in self.pots:
                pot.set_Coulomb(self)
        for pot in self.pots:
            pot.register_data(self.particles)
        self.particles.allocate(len(atoms))
        self.particles.set_cell(atoms.get_cell(), pbc)
        
        if self._lees_edwards_info_str in atoms.info:
            self._warn_lees_edwards()
            self.lees_edwards_dx = atoms.info[self._lees_edwards_info_str]
            assert self.lees_edwards_dx.shape == (3,)

        if self.lees_edwards_dx is not None:
            self.particles.set_lees_edwards(self.lees_edwards_dx, self.lees_edwards_dv)

        Z  = self.particles.Z

        for i, at in enumerate(atoms):
            Z[i]   = atomic_numbers[at.symbol]

        self.particles.coordinates[:, :]  = atoms.get_positions()[:, :]
        # Notify the Particles object of a change
        self.particles.I_changed_positions()

        self.particles.update_elements()

        # Initialize and set neighbor list
        self.nl = _atomistica.Neighbors(self.avgn)

        # Tell the potential about the new Particles and Neighbors object
        for pot in self.pots:
            pot.bind_to(self.particles, self.nl)
        for coul in self.couls:
            coul.bind_to(self.particles, self.nl)

        # Force re-computation of energies/forces the next time these are
        # requested
        self.force_update = True


    def set_per_at(self, epot=None, wpot=None):
        if epot is not None:
            self.compute_epot_per_at = epot
        if wpot is not None:
            self.compute_wpot_per_at = wpot
        self.force_update = True


    def set_per_bond(self, epot=None, f=None, wpot=None):
        if epot is not None:
            self.compute_epot_per_bond = epot
        if f is not None:
            self.compute_f_per_bond = f
        if wpot is not None:
            self.compute_wpot_per_bond = wpot
        self.force_update = True


    def update(self, atoms, force_update=False):
        if atoms is None:
            return

        # Number of particles changed? -> Reinitialize potential
        if self.particles is None or len(self.particles.Z) != len(atoms):
            self.initialize(atoms)
        # Type of particles changed? -> Reinitialize potential
        elif np.any(self.particles.Z != atoms.get_atomic_numbers()):
            self.initialize(atoms)

        # Cell changed? FIXME! Add PBC,LEBC changed
        cell_chgd  = False
        cell       = self.particles.cell
        pbc        = atoms.get_pbc()
        #if np.any(np.abs(cell - atoms.get_cell()) > self.CELL_TOL):
        if np.any(cell != atoms.get_cell()):
            cell[:, :]  = atoms.get_cell()[:, :]
            self.particles.set_cell(cell, pbc)
            cell_chgd  = True

        pos_chgd   = False
        positions  = self.particles.coordinates
        scaled     = np.linalg.solve(cell.T, positions.T).T
        for i in range(3):
            if pbc[i]:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test
                scaled[:, i] %= 1.0
                scaled[:, i] %= 1.0

        #if np.any(np.abs(scaled - atoms.get_scaled_positions()) > self.POSITIONS_TOL):
        if np.any(scaled != atoms.get_scaled_positions()):
            positions[:, :]  = atoms.get_positions()[:, :]
            pos_chgd         = True
            # Notify the Particles object of a change
            self.particles.I_changed_positions()

        lebc_chgd = False            
        if self._lees_edwards_info_str in atoms.info:
            if self.lees_edwards_dx is None or \
                np.any(atoms.info[self._lees_edwards_info_str] != \
                       self.lees_edwards_dx):

                self._warn_lees_edwards()
                self.lees_edwards_dx = atoms.info[self._lees_edwards_info_str]
                lebc_chgd = True
                assert self.lees_edwards_dx.shape == (3,)

        if self.lees_edwards_dx is not None:
            self.particles.set_lees_edwards(self.lees_edwards_dx,
                                            self.lees_edwards_dv)

        if pos_chgd or cell_chgd or lebc_chgd or self.force_update or \
            force_update:

            self.calculate()
            self.force_update  = False

            if self.q is not None:
               atoms.set_charges(self.q)


    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        return self.epot


    def get_forces(self, atoms):
        self.update(atoms)

        return self.forces.copy()


    def get_stress(self, atoms):
        self.update(atoms)

        st = self.wpot/atoms.get_volume()

        return np.array( [ st[0,0], st[1,1], st[2,2], (st[1,2]+st[2,1])/2,
                           (st[0,2]+st[2,0])/2, (st[0,1]+st[1,0])/2 ] )


    def get_charges(self, atoms=None):
        self.update(atoms)

        return self.q


    def get_electrostatic_potential(self, atoms=None):
        self.update(atoms)

        return self.phi


    def get_electrostatic_field(self, atoms=None):
        if atoms is not None:
            self.update(atoms)

        return self.E


    def calculate(self):
        self.epot = 0.0
        self.forces = np.zeros([len(self.particles),3])
        self.wpot = np.zeros([3,3])

        if self.q is None:
            # No charges? Just call the potentials...
            for pot in self.pots:
                _epot, _forces, _wpot, self.epot_per_at, self.epot_per_bond, \
                    self.f_per_bond, self.wpot_per_at, self.wpot_per_bond  = \
                        pot.energy_and_forces(
                            self.particles, self.nl, forces = self.forces,
                            epot_per_at = self.compute_epot_per_at,
                            epot_per_bond = self.compute_epot_per_bond,
                            f_per_bond = self.compute_f_per_bond,
                            wpot_per_at = self.compute_wpot_per_at,
                            wpot_per_bond = self.compute_wpot_per_bond)

                self.epot += _epot
                self.wpot += _wpot

        else:
            # Charges? Pass charge array to potentials and ...
            for pot in self.pots:
                _epot, _forces, _wpot, self.epot_per_at, self.epot_per_bond, \
                    self.f_per_bond, self.wpot_per_at, self.wpot_per_bond  = \
                        pot.energy_and_forces(
                            self.particles, self.nl, forces = self.forces,
                            charges = self.q,
                            epot_per_at = self.compute_epot_per_at,
                            epot_per_bond = self.compute_epot_per_bond,
                            f_per_bond = self.compute_f_per_bond,
                            wpot_per_at = self.compute_wpot_per_at,
                            wpot_per_bond = self.compute_wpot_per_bond)

                self.epot += _epot
                self.wpot += _wpot

            # ... call Coulomb solvers to get potential and fields
            self.phi = np.zeros(len(self.q))
            epot_coul = 0.0
            self.E = np.zeros([len(self.q),3])
            wpot_coul = 0.0

            for coul in self.couls:
                _phi, _epot, _E, _wpot = \
                    coul.potential_and_field(self.particles, self.nl, self.q,
                                             self.phi, self.E)
                epot_coul += _epot
                wpot_coul += _wpot

            # Convert units
            self.phi *= Hartree * Bohr
            epot_coul *= Hartree * Bohr
            self.E *= Hartree * Bohr
            wpot_coul *= Hartree * Bohr

            self.epot += epot_coul
            self.wpot += wpot_coul

            # Forces are charges time field
            self.forces += self.q.reshape(-1,1) * self.E
            
            
    ### Covenience
    def _warn_lees_edwards(self):
        global _warned_about_lees_edwards
        if not _warned_about_lees_edwards:
            print "Warning: Setting Lees-Edwards boundary conditions from " \
                "information found in the %s entry of the Atoms.info " \
                "dictionary. Is this really what you intended?" % \
                self._lees_edwards_info_str
            _warned_about_lees_edwards = True


    ### Atomistica features
    def set_lees_edwards(self, dx, dv=None):
        self.lees_edwards_dx  = dx
        if dv is None:
            self.lees_edwards_dv  = [ 0.0, 0.0, 0.0 ]
        else:
            self.lees_edwards_dv  = dv

        if self.particles is not None and self.lees_edwards_dx is not None:
            self.particles.set_lees_edwards(self.lees_edwards_dx,
                                            self.lees_edwards_dv)
            self.force_update = True



    def get_atomic_stress(self):
        r = np.zeros( [ 3, 3, len(self.particles) ] )
        for i, a in enumerate(self.particles):
            r[:, :, i] = a.w
        return r


    def get_neighbors(self):
        return self.nl.get_neighbors(self.particles)


    def __str__(self):
        s = 'Atomistica(['
        for pot in self.pots:
            s += pot.__str__()+','
        for coul in self.couls:
            s += pot.__str__()+','
        return s[:-1]+'])'


    ### Coulomb solver callback
    def set_Hubbard_U(self, p, U):
        assert p is self.particles
        for coul in self.couls:
            coul.set_Hubbard_U(p, U)


    def potential_and_field(self, p, nl, q, phi, epot, E, wpot):
        # It appeas this callback is actually not needed. Remove?
        raise NotImplementedError


    def potential(self, p, nl, q, phi):
        assert p is self.particles
        assert nl is self.nl
        _phi = np.zeros_like(phi)
        for coul in self.couls:
            coul.potential(p, nl, q, _phi)
        phi += _phi * Hartree * Bohr
    

### Construct ASE interface wrappers for all potentials in _atomistica

exclude_list = [ 'TightBinding' ]
spec_avgn = dict(Gupta=1000)

for name, cls in inspect.getmembers(_atomistica):
    if hasattr(cls, 'energy_and_forces') and not cls.__name__ in exclude_list:
        avgn = 100
        if cls.__name__ in spec_avgn.keys():
            avgn = spec_avgn[cls.__name__]
        elif cls.__name__.endswith('Scr'):
            avgn = 1000
        globals()[cls.__name__] = type(cls.__name__, (Atomistica,), {
            'potential_class': cls,
            'avgn': avgn
            })

# The tight-binding module needs special attention to make it somewhat easier
# to use.

if hasattr(_atomistica, 'TightBinding'):
    class TightBinding(Atomistica):
        """Non-othogonal tight-binding.
        """

        potential_class = _atomistica.TightBinding

        def __init__(self, width=0.1, kpts=None, kpts_shift=None,
                     database_folder=None):
            
            # Translate from a human friendly format
            d = {
                "SolverLAPACK" : {
                    "electronic_T" : width
                    }
                }

            if database_folder is not None:
                d["database_folder"] = database_folder

            if kpts is not None:
                bz = { }
                bz["nk"] = kpts

                if kpts_shift is not None:
                    bz["shift"] = kpts_shift

                d["BrillouinZone"] = { "MonkhorstPack" : bz }

            d['avgn'] = 1000

            apply(Atomistica.__init__, (self,), d)

