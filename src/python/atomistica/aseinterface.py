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
ASE interface to Atomistica.
"""

import copy
import inspect
from math import sqrt, log

import _atomistica

import numpy as np

from ase.atoms import string2symbols
from ase.data import atomic_numbers
from ase.units import Hartree, Bohr

from atomistica.parameters import *

###

def convpar(p):
    """
    Convert a parameter set from convenient Python dictionary to the format
    expected by the Fortran kernels.
    """

    if 'el' not in p:
        return p

    els = p['el']
    nel = len(els)

    q = { }
    for name, values in p.iteritems():
        
        if isinstance(values, dict):
            # This is a dictionary. We need to first understand what it is and
            # then convert to the appropriate list
            m = 0
            for itype, value in values:
                if itype[0] == '_':
                    continue
                # Ask ASE to turn this into a list of symbols
                syms = string2symbols(itype)
                # Now turn this into a list of indices
                nums = [ els.index(sym) for sym in syms ]
                # We need to know the maximum number of symbols
                m = max(m, nums)

            default = 0.0
            if '__default__' in values:
                default = values['__default__']
            if m == 2:
                # These are pair indices
                new_values = [ default ]*pair_index(nel,nel,nel)
            elif m == 3:
                # These are triplet indices
                new_values = [ default ]*triplet_index(nel,nel,nel,nel)
            else:
                raise RuntimeError("Parameter '%s' appears to involve " \
                                   "interaction between %i atoms. Do not know "\
                                   "how to handle this." % ( name, nel ))

            for itype, value in values:
                if itype[0] == '_':
                    continue
                # Ask ASE to turn this into a list of symbols
                syms = string2symbols(itype)
                # Now turn this into a list of indices
                nums = [ els.index(sym) for sym in syms ]
                if len(nums) == m:
                    # Set value
                    if m == 2:
                        i,j = nums
                        new_values[pair_index(i,j,nel)] = value
                    elif m == 3:
                        i,j,k = nums
                        new_values[triple_index(i,j,k,nel)] = value
                else:
                    # We have fewer values than we need
                    if m == 3 and len(nums) == 1:
                        [k] = nums
                        for i in range(nel):
                            for j in range(nel):
                                # There is six possible permutations
                                new_values[triple_index(i,j,k,nel)] = value
                                new_values[triple_index(i,k,j,nel)] = value
                                new_values[triple_index(k,j,i,nel)] = value
                                new_values[triple_index(j,k,i,nel)] = value
                                new_values[triple_index(k,i,j,nel)] = value
                                new_values[triple_index(i,j,k,nel)] = value
                    elif m == 3 and len(nums) == 2:
                        [i,j] = nums
                        for k in range(nel):
                            # There is six possible permutations
                            new_values[triple_index(i,j,k,nel)] = value
                            new_values[triple_index(i,k,j,nel)] = value
                            new_values[triple_index(k,j,i,nel)] = value
                            new_values[triple_index(j,k,i,nel)] = value
                            new_values[triple_index(k,i,j,nel)] = value
                            new_values[triple_index(i,j,k,nel)] = value
                    else:
                        raise RuntimeError("Parameter '%s' appears to involve " \
                                           "interaction between %i atoms, but " \
                                           "only %i elements were provied. Do " \
                                           "not know how to handle this." \
                                           % ( name, nel, len(nums) ))
            values = new_values

        q[name] = values
    return q

###

_warned_about_lees_edwards = False
class Atomistica(object):
    """
    Atomistica ASE calculator.
    """

    CELL_TOL = 1e-16
    POSITIONS_TOL = 1e-16

    name = None
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
                if isinstance(pot, Atomistica):
                    raise TypeError('Potential passed to Atomistica class is '
                                    'already an Atomistica object. This does '
                                    'not work. You need to pass native '
                                    'potential from atomistica.native or '
                                    '_atomistica.')
                if hasattr(pot, 'potential'):
                    self.couls += [ pot ]
                else:
                    self.pots += [ pot ]
        else:
            pot = self.potential_class(**convpar(kwargs))
            if hasattr(pot, 'potential'):
                self.couls += [ pot ]
            else:
                self.pots += [ pot ]

        if avgn:
            self.avgn = avgn

        self.particles = None
        self.nl = None

        self.q = None

        self.mask = None

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


    def todict(self):
        return self.kwargs


    def check_state(self, atoms):
        # This is a hack to indicate ase.db that state has changed
        return [1]
        

    def initialize(self, atoms):
        if self.mask is not None:
            if len(self.mask) != len(atoms):
                raise RuntimeError('Length of mask array (= {}) does not equal '
                                   'number of atoms (= {}).'
                                   .format(len(self.mask), len(atoms)))

        pbc = atoms.get_pbc()
        self.particles = _atomistica.Particles()

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

        if len(self.couls) > 0:
            if atoms.has('charges'):
                self.q = atoms.get_array('charges')
            else:
                self.q = atoms.get_charges()
            if self.q is None:
                self.q = np.zeros(len(atoms))
            self.E = np.zeros([3,len(atoms)])
            # Coulomb callback should be directed to this wrapper object
            for pot in self.pots:
                pot.set_Coulomb(self)

        for coul in self.couls:
            coul.bind_to(self.particles, self.nl)

        # Force re-computation of energies/forces the next time these are
        # requested
        self.force_update = True


    def set_mask(self, mask):
        self.mask = mask
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
        pbc        = self.particles.pbc
        #if np.any(np.abs(cell - atoms.get_cell()) > self.CELL_TOL):
        if np.any(cell != atoms.get_cell()) or np.any(pbc != atoms.get_pbc()):
            self.particles.set_cell(atoms.get_cell(), atoms.get_pbc())
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

        # Charges changed?
        charges_chgd = False
        if self.q is not None:
            if atoms.has('charges'):
                q = atoms.get_array('charges')
            else:
                q = atoms.get_charges()
            if np.any(self.q != q):
                self.q = q
                charges_chgd = True

        if pos_chgd or cell_chgd or lebc_chgd or charges_chgd or \
                self.force_update or force_update:

            self.calculate()
            self.force_update  = False

        if self.q is not None:
            atoms.set_initial_charges(self.q)

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        return self.epot


    def get_potential_energies(self, atoms):
        self.update(atoms)
        if not self.compute_epot_per_at:
            self.compute_epot_per_at = True
            self.calculate()
        return self.epot_per_at


    def get_forces(self, atoms):
        self.update(atoms)

        return self.forces.copy()


    def get_stress(self, atoms):
        self.update(atoms)

        st = self.wpot/atoms.get_volume()

        return np.array( [ st[0,0], st[1,1], st[2,2], (st[1,2]+st[2,1])/2,
                           (st[0,2]+st[2,0])/2, (st[0,1]+st[1,0])/2 ] )


    def get_stresses(self, atoms):
        self.update(atoms)
        if not self.compute_wpot_per_at:
            self.compute_wpot_per_at = True
            self.calculate()
        return np.transpose([
                self.wpot_per_at[:,0,0],
                self.wpot_per_at[:,1,1],
                self.wpot_per_at[:,2,2],
                (self.wpot_per_at[:,1,2]+self.wpot_per_at[:,2,1])/2,
                (self.wpot_per_at[:,0,2]+self.wpot_per_at[:,2,0])/2,
                (self.wpot_per_at[:,0,1]+self.wpot_per_at[:,1,0])/2])


    def get_charges(self, atoms=None):
        self.update(atoms)

        return self.q


    def get_electrostatic_potential(self, atoms=None):
        self.phi = np.zeros(len(self.particles))
        self.potential(self.particles, self.nl, self.q, self.phi)

        return self.phi


    def get_per_bond_property(self, name):
        for pot in self.pots:
            if hasattr(pot, 'get_per_bond_property'):
                return pot.get_per_bond_property(self.particles, self.nl, name)


    def calculate(self):
        self.epot = 0.0
        self.forces = np.zeros([len(self.particles),3])
        self.wpot = np.zeros([3,3])

        kwargs = dict(epot_per_at = self.compute_epot_per_at,
                      epot_per_bond = self.compute_epot_per_bond,
                      f_per_bond = self.compute_f_per_bond,
                      wpot_per_at = self.compute_wpot_per_at,
                      wpot_per_bond = self.compute_wpot_per_bond)

        if self.mask is not None:
            kwargs['mask'] = self.mask
        if self.q is None:
            # No charges? Just call the potentials...
            for pot in self.pots:
                _epot, _forces, _wpot, self.epot_per_at, self.epot_per_bond, \
                    self.f_per_bond, self.wpot_per_at, self.wpot_per_bond  = \
                        pot.energy_and_forces(self.particles, self.nl,
                                              forces = self.forces,
                                              **kwargs)
                self.epot += _epot
                self.wpot += _wpot

        else:
            # Charges? Pass charge array to potentials and ...
            for pot in self.pots:
                _epot, _forces, _wpot, self.epot_per_at, self.epot_per_bond, \
                    self.f_per_bond, self.wpot_per_at, self.wpot_per_bond  = \
                        pot.energy_and_forces(self.particles, self.nl,
                                              forces = self.forces,
                                              charges = self.q,
                                              **kwargs)
                self.epot += _epot
                self.wpot += _wpot

            # ... call Coulomb solvers to get potential and fields
            epot_coul = 0.0
            forces_coul = np.zeros([len(self.particles),3])
            wpot_coul = np.zeros([3,3])

            for coul in self.couls:
                _epot, _forces, _wpot = \
                    coul.energy_and_forces(self.particles, self.nl, self.q,
                                           forces_coul)
                epot_coul += _epot
                wpot_coul += _wpot

            # Convert units
            epot_coul *= Hartree * Bohr
            forces_coul *= Hartree * Bohr
            wpot_coul *= Hartree * Bohr

            self.epot += epot_coul
            self.wpot += wpot_coul

            # Sum forces
            self.forces += forces_coul
            
            
    ### Convenience
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
            s += coul.__str__()+','
        return s[:-1]+'])'


    ### Coulomb solver callback
    def set_Hubbard_U(self, p, U):
        assert p is self.particles
        for coul in self.couls:
            coul.set_Hubbard_U(p, U)


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
    if hasattr(cls, 'energy_and_forces') and \
            not cls.__name__ in exclude_list:
        avgn = 100
        if cls.__name__ in spec_avgn.keys():
            avgn = spec_avgn[cls.__name__]
        elif cls.__name__.endswith('Scr'):
            avgn = 1000
        globals()[cls.__name__] = type(cls.__name__, (Atomistica,), {
            'name': cls.__name__,
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

        def __init__(self, width=0.1, database_folder=None):
            
            # Translate from a human friendly format
            d = {
                "SolverLAPACK" : {
                    "electronic_T" : width
                    }
                }

            if database_folder is not None:
                d["database_folder"] = database_folder

            d['avgn'] = 1000

            apply(Atomistica.__init__, (self,), d)

