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

"""
Test computation of stresses from tight binding.
"""

import os
import sys
import unittest

import numpy as np

import ase.io as io
from ase.build import molecule
from ase.optimize import FIRE

import atomistica.native as native
from atomistica import Atomistica

def test_stresses(test=None):
    database_folder = os.getenv('MIO')
    if database_folder is None:
        raise RuntimeError('Please use environment variable MIO to specify path to mio Slater-Koster tables.')

    calc = Atomistica(
        [ native.TightBinding(
        database_folder = database_folder,
        SolverLAPACK = dict(electronic_T=0.001),
        SCC = dict(dq_crit = 1e-4,
                   mixing = 0.2, # 0.2
                   andersen_memory = 3, # 3
                   maximum_iterations = 100,
                   log = True)
        ),
          native.DirectCoulomb(),
          native.SlaterCharges(cutoff=10.0) ],
        avgn = 1000
        )

    a = io.read('aC_small.cfg')
    a.calc = calc
    s = a.get_stress() * a.get_volume()
    s_at = a.get_stresses()

    np.testing.assert_array_almost_equal(s, s_at.sum(axis=0))

###

if __name__ == '__main__':
    test_stresses()
