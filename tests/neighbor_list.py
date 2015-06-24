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
#! /usr/bin/env python


import unittest

import numpy as np

import ase

import atomistica.io as io
import atomistica.native as native
from atomistica.snippets import mic

###

class NeighborListTest(unittest.TestCase):

    def test_neighbor_list(self):
        a = io.read('aC.cfg')
        an = native.from_atoms(a)
        nl = native.Neighbors(100)
        nl.request_interaction_range(5.0)

        i, j, abs_dr_no_vec = nl.get_neighbors(an)
        i, j, dr, abs_dr = nl.get_neighbors(an, vec=True)

        self.assertTrue(np.all(np.abs(abs_dr_no_vec-abs_dr) < 1e-12))

        r = a.get_positions()
        dr_direct = mic(r[i]-r[j], a.cell)

        abs_dr_from_dr = np.sqrt(np.sum(dr*dr, axis=1))
        abs_dr_direct = np.sqrt(np.sum(dr_direct*dr_direct, axis=1))

        self.assertTrue(np.all(np.abs(abs_dr-abs_dr_from_dr) < 1e-12))
        self.assertTrue(np.all(np.abs(abs_dr-abs_dr_direct) < 1e-12))

        self.assertTrue(np.all(np.abs(dr-dr_direct) < 1e-12))

    def test_pbc(self):
        a = ase.Atoms('CC',
                      positions=[[0.1, 0.5, 0.5],
                                 [0.9, 0.5, 0.5]],
                      cell=[1, 1, 1],
                      pbc=True)
        an = native.from_atoms(a)
        nl = native.Neighbors(100)
        nl.request_interaction_range(0.3)

        # with pbc

        i, j, abs_dr = nl.get_neighbors(an)
        self.assertEqual(len(i), 2)

        a.set_pbc(False)
        an = native.from_atoms(a)
        nl = native.Neighbors(100)
        nl.request_interaction_range(0.3)

        # no pbc

        i, j, abs_dr = nl.get_neighbors(an)
        self.assertEqual(len(i), 0)

        a.set_pbc([False,False,True])
        an = native.from_atoms(a)
        nl = native.Neighbors(100)
        nl.request_interaction_range(0.3)

        # partial pbc

        i, j, abs_dr = nl.get_neighbors(an)
        self.assertEqual(len(i), 0)

        a.set_pbc([True,False,False])
        an = native.from_atoms(a)
        nl = native.Neighbors(100)
        nl.request_interaction_range(0.3)

        # partial pbc

        i, j, abs_dr = nl.get_neighbors(an)
        self.assertEqual(len(i), 2)

###

if __name__ == '__main__':
    unittest.main()
