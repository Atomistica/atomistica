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

import numpy as np

class JoinCalculators:
    """Create a joint calculator, i.e. one can augment a DFT calculation
       by adding a classical potential for van-der-Waals interactions.

       Potential energies, forces, etc. are simply summed up.
    """

    def __init__(self, calcs):
        self.calcs = calcs


    def get_forces(self, a):
        """Calculate atomic forces."""
        f = np.zeros( [ len(a), 3 ], dtype=float )
        for c in self.calcs:
            f += c.get_forces(a)
        return f


    def get_magnetic_moments(self, a):
        """Get calculated local magnetic moments."""
        raise NotImplementedError


    def get_potential_energy(self, a):
        """Calculate potential energy."""
        e = 0.0
        for c in self.calcs:
            e += c.get_potential_energy(a)
        return e


    def get_potential_energies(self, a):
        """Calculate the potential energies of all the atoms."""
        raise NotImplementedError


    def get_spin_polarized(self):
        """Get calculated total magnetic moment."""
        raise NotImplementedError


    def get_stress(self, a):
        """Calculate stress tensor."""
        s = np.zeros( 6, dtype=float )
        for c in self.calcs:
            s += c.get_stress(a)
        return s


    def get_stresses(self, a):
        """Calculate the stress-tensor of all the atoms."""
        raise NotImplementedError


    def set_atoms(self, a):
        """Assign an atoms object."""
        for c in self.calcs:
            if hasattr(c, "set_atoms"):
                c.set_atoms(a)


class LinearPotential:
    """ Potential that is linear in some direction, i.e. a constant force
    """
    def __init__(self, force, mask=None):
        self.force = force
        self.mask = mask


    def get_forces(self, a=None):
        """Calculate atomic forces."""
        if a is None:
            a = self.a
        forces = np.zeros([len(a), 3], dtype=float)
        if self.mask is None:
            forces[self.mask] = self.force
        else:
            forces[:] = self.force
        return forces


    def get_potential_energy(self, a=None):
        """Calculate potential energy."""
        if a is None:
            a = self.a
        if self.mask is None:
            return -np.sum(np.dot(a.get_positions()[self.mask], self.force))
        else:
            return -np.sum(np.dot(a.get_positions(), self.force))


    def get_potential_energies(self, a):
        """Calculate the potential energies of all the atoms."""
        raise NotImplementedError


    def get_stress(self, a=None):
        """Calculate stress tensor."""
        return np.zeros(6, dtype=float)


    def set_atoms(self, a):
        """Assign an atoms object."""
        self.a = a
