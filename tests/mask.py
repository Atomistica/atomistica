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

import ase.io as io
from ase.lattice.cubic import Diamond, FaceCenteredCubic

import atomistica
from atomistica import LJCut, TabulatedAlloyEAM, Tersoff, TersoffScr

###

class MaskTest(unittest.TestCase):

    def random_mask_test(self, a):
        c = a.calc
        e = a.get_potential_energy()
        f = a.get_forces()
        w = a.get_stress()

        mask = np.random.randint(0, len(a), size=len(a)) < \
            len(a)/2
        imask = np.logical_not(mask)

        c.set_mask(mask)
        e1 = a.get_potential_energy()
        f1 = a.get_forces()
        w1 = a.get_stress()

        c.set_mask(imask)
        e2 = a.get_potential_energy()
        f2 = a.get_forces()
        w2 = a.get_stress()

        c.set_mask(None)
        e3 = a.get_potential_energy()

        self.assertTrue(abs(e-e1-e2) < 1e-6)
        self.assertTrue(abs(e-e3) < 1e-6)
        self.assertTrue(np.max(np.abs(f-f1-f2)) < 1e-6)
        self.assertTrue(np.max(np.abs(w-w1-w2)) < 1e-6)

    def test_mask_decomposition_bop(self):
        a = io.read('aC.cfg')
        for pot in [Tersoff, TersoffScr]:
            c = Tersoff()
            a.set_calculator(c)
            self.random_mask_test(a)

    def test_mask_decomposition_lj_cut(self):
        a = FaceCenteredCubic('Au', size=[2,2,2])
        c = LJCut(el1='Au', el2='Au', epsilon=1.0, sigma=1.0, cutoff=6.0)
        a.set_calculator(c)
        self.random_mask_test(a)

    def test_mask_decomposition_tabulated_alloy_eam(self):
        a = FaceCenteredCubic('Au', size=[2,2,2])
        c = TabulatedAlloyEAM(fn='Au-Grochola-JCP05.eam.alloy')
        a.set_calculator(c)
        self.random_mask_test(a)

###

if __name__ == '__main__':
    unittest.main()
