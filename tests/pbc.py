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
import ase.io
from ase.lattice.cubic import Diamond

from atomistica import Tersoff

###

Jm2 = 1e23/ase.units.kJ

###

class PBCTest(unittest.TestCase):

    def test_pbc(self):
        a = Diamond('Si', latticeconstant=5.432, size=[2,2,2])
        sx, sy, sz = a.get_cell().diagonal()
        a.set_calculator(Tersoff())
        e1 = a.get_potential_energy()

        a.set_pbc([True,True,False])
        e2 = a.get_potential_energy()

        a.set_pbc(True)
        a.set_cell([sx,sy,2*sz])
        e3 = a.get_potential_energy()

        self.assertEqual(e2, e3)

        # This should give the unrelaxed surface energy
        esurf = (e2-e1)/(2*sx*sy) * Jm2
        self.assertTrue(abs(esurf-2.309) < 0.001)

###

if __name__ == '__main__':
    unittest.main()
