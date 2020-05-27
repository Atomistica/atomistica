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

import atomistica.io as io
import atomistica.native as native
from atomistica import Tersoff
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

    def test_pbc_shift_by_multiple_cells(self):
        a = io.read('aC.cfg')
        a.set_calculator(Tersoff())
        e1 = a.get_potential_energy()
        i1, j1, r1 = a.calc.nl.get_neighbors(a.calc.particles)
        a[100].position += 3*a.cell[0]
        e2 = a.get_potential_energy()
        i2, j2, r2 = a.calc.nl.get_neighbors(a.calc.particles)
        for i in range(len(a)):
            n1 = np.array(sorted(j1[i1==i]))
            n2 = np.array(sorted(j2[i2==i]))
            if np.any(n1 != n2):
                print(i, n1, n2)
        a[100].position += a.cell.T.dot([1,3,-4])
        e3 = a.get_potential_energy()
        self.assertAlmostEqual(e1, e2)
        self.assertAlmostEqual(e1, e3)

    def test_no_pbc_small_cell(self):
        a = io.read('aC.cfg')
        a.set_calculator(Tersoff())
        a.set_pbc(False)
        e1 = a.get_potential_energy()
        i1, j1, r1 = a.calc.nl.get_neighbors(a.calc.particles)
        a.set_cell(a.cell*0.9, scale_atoms=False)
        e2 = a.get_potential_energy()
        self.assertAlmostEqual(e1, e2)
        i2, j2, r2 = a.calc.nl.get_neighbors(a.calc.particles)
        for k in range(len(a)):
            neigh1 = np.array(sorted(j1[i1==k]))
            neigh2 = np.array(sorted(j2[i2==k]))
            self.assertTrue(np.all(neigh1 == neigh2))

    def test_partial_pbc_small_cell(self):
        a = io.read('aC.cfg')
        a.set_cell(a.cell.diagonal(), scale_atoms=True)
        a.set_calculator(Tersoff())
        a.set_pbc([True, False, False])
        e1 = a.get_potential_energy()
        i1, j1, r1 = a.calc.nl.get_neighbors(a.calc.particles)
        a.set_cell(a.cell.diagonal()*np.array([1.0, 0.8, 0.9]), scale_atoms=False)
        e2 = a.get_potential_energy()
        self.assertAlmostEqual(e1, e2)
        i2, j2, r2 = a.calc.nl.get_neighbors(a.calc.particles)
        for k in range(len(a)):
            neigh1 = np.array(sorted(j1[i1==k]))
            neigh2 = np.array(sorted(j2[i2==k]))
            self.assertTrue(np.all(neigh1 == neigh2))

    def test_floating_point_issue(self):
        calc = Tersoff()
        a1 = ase.Atoms('Si4C4', positions=np.array([[-4.41173839e-52,  0.00000000e+00,  0.00000000e+00],
                                                    [-4.41173839e-52,  2.26371743e+00,  2.26371743e+00],
                                                    [ 2.26371743e+00,  0.00000000e+00,  2.26371743e+00],
                                                    [ 2.26371743e+00,  2.26371743e+00,  0.00000000e+00],
                                                    [ 1.13185872e+00,  1.13185872e+00,  1.13185872e+00],
                                                    [ 1.13185872e+00,  3.39557615e+00,  3.39557615e+00],
                                                    [ 3.39557615e+00,  1.13185872e+00,  3.39557615e+00],
                                                    [ 3.39557615e+00,  3.39557615e+00,  1.13185872e+00]]),
                   cell=[4.527434867899659, 4.527434867899659, 4.527434867899659], pbc=True)

        a1.calc = calc
        a1.get_potential_energy()
        self.assertTrue((calc.nl.get_coordination_numbers(calc.particles, 3.0) == 4).all())

        a2 = a1.copy()
        a2.calc = calc
        a2.set_scaled_positions(a2.get_scaled_positions())
        a2.get_potential_energy()
        self.assertTrue((calc.nl.get_coordination_numbers(calc.particles, 3.0) == 4).all())

###

if __name__ == '__main__':
    unittest.main()
