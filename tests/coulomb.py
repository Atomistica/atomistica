# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
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
#! /usr/bin/env python


import unittest

import numpy as np

import ase
from ase.units import Hartree, Bohr

import atomistica
from atomistica import DirectCoulomb

###

class CoulombTest(unittest.TestCase):

    def test_direct_coulomb(self):
        a = ase.Atoms('NaCl', positions=[[-1.0, 0, 0], [1.0, 0, 0]], pbc=False)
        a.center(vacuum=10.0)

        a.set_initial_charges(np.zeros(len(a)))

        c = DirectCoulomb()
        a.set_calculator(c)

        assert a.get_potential_energy() == 0
        assert (np.abs(c.get_electrostatic_potential()) < 1e-9).all()

        a.set_initial_charges([-1,1])

        c = DirectCoulomb()
        a.set_calculator(c)

        assert abs(a.get_potential_energy()+Hartree*Bohr/2) < 1e-9
        assert (np.abs(c.get_electrostatic_potential()-Hartree*Bohr/2*np.array([1,-1])) < 1e-9).all()
    
###

if __name__ == '__main__':
    unittest.main()
