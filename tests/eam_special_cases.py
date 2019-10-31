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


"""
Test special cases where EAM potential may fail.
"""

from __future__ import print_function

import sys

import unittest

import numpy as np

import ase.io as io

from atomistica import TabulatedAlloyEAM
from atomistica.tests import test_forces

###

dx = 1e-6
tol = 1e-6

###

class TestEAMSpecialCases(unittest.TestCase):

    def test_crash1(self):
        a = io.read('eam_crash1.poscar')
        a.set_calculator(TabulatedAlloyEAM(fn='Cu_mishin1.eam.alloy'))
        a.get_potential_energy()

    def test_dense_forces(self):
        orig_a = io.read('eam_crash2.poscar')
        c = TabulatedAlloyEAM(fn='Cu_mishin1.eam.alloy')
        for fac in [0.2, 0.3, 0.4, 0.5]:
            a = orig_a.copy()
            a.set_cell(fac*a.cell, scale_atoms=True)
            a.set_calculator(c)
            ffd, f0, maxdf = test_forces(a, dx=dx)
            if maxdf > tol:
                nfail += 1
                print("forces .failed.")
                print("max(df)  = %f" % maxdf)
                print("f - from potential")
                for i, f in enumerate(f0):
                    print(i, f)
                print("f - numerically")
                for i, f in enumerate(ffd):
                    print(i, f)
                print("difference between the above")
                for i, f in enumerate(f0-ffd):
                    print(i, f)
            self.assertTrue(maxdf < tol)

###

if __name__ == '__main__':
    unittest.main()
