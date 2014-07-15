#! /usr/bin/env python

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

import unittest

import numpy as np

import atomistica.io as io
import atomistica.native as native
from atomistica.snippets import mic

###

class NeighborListTest(unittest.TestCase):

    def test_neighbor_list(self):
        a = io.read('aC.traj')
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

###

if __name__ == '__main__':
    unittest.main()
