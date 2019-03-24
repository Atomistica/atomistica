#! /usr/bin/env python

# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2019) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
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
Test the DFTB3 parametrizations.
"""

from __future__ import print_function

import os
import sys
import unittest

import numpy as np

import ase.units
from ase.build import molecule
from ase.io import read
from ase.optimize import FIRE

import atomistica.native as native
from atomistica import Atomistica

###

# Table 5 from Gaus, Cui, Elstner, JCTC 7, 931 (2011)
table5_data = {
    '2H2O': {
        'reference_structures': ['H2O', 'H2O'],
        'structure': 'molecule_database/2H2O.xyz',
        'G3B3': -4.9, # kcal/mol
        'DFTB2': 1.6,
        'DFTB3': 1.5
    },
    '3H2O': {
        'reference_structures': ['H2O', 'H2O', 'H2O'],
        'structure': 'molecule_database/3H2O.xyz',
        'G3B3': -15.1, # kcal/mol
        'DFTB2': 5.5,
        'DFTB3': 5.4
    },
    '4H2O': {
        'reference_structures': ['H2O', 'H2O', 'H2O', 'H2O'],
        'structure': 'molecule_database/4H2O.xyz',
        'G3B3': -27.4, # kcal/mol
        'DFTB2': 9.7,
        'DFTB3': 9.4
    },
    '5H2O': {
        'reference_structures': ['H2O', 'H2O', 'H2O', 'H2O', 'H2O'],
        'structure': 'molecule_database/5H2O.xyz',
        'G3B3': -36.3, # kcal/mol
        'DFTB2': 13.3,
        'DFTB3': 12.5
    }
}

###

def run_dftb3_test(test=None):
    mio_database_folder = os.getenv('MIO')
    if mio_database_folder is None:
        raise RuntimeError('Please use environment variable MIO to specify path to mio Slater-Koster tables.')

    mio_calc = Atomistica(
        [ native.TightBinding(
        database_folder = mio_database_folder,
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

    for name, data in table5_data.items():
        e0 = 0.0
        for structure in data['reference_structures']:
            if os.path.exists(structure):
                a = read(structure)
            else:
                a = molecule(structure)
            a.center(vacuum=10.0)
            a.set_calculator(mio_calc)
            FIRE(a).run(fmax=0.001)
            e0 += a.get_potential_energy()

        eref_G3B3 = data['G3B3']
        eref_DFTB2 = data['DFTB2']
        eref_DFTB3 = data['DFTB3']

        a = read(data['structure'])
        a.center(vacuum=10.0)
        a.set_calculator(mio_calc)
        FIRE(a).run(fmax=0.001)
        e = a.get_potential_energy()

        e = (e - e0)/(ase.units.kcal/ase.units.mol)

        if test is None:
            print('{0:>20} {1:>20.10f} {2:>20.10f} {3:>20.10f}'.format(name, eref_G3B3, e - eref_G3B3, eref_DFTB2))
        else:
            test.assertAlmostEqual(e - eref_G3B3, eref_DFTB2)

###

class TestMIO(unittest.TestCase):

    def test_dftb3(self):
        if os.getenv('MIO') is None:
            print('Skipping DFTB3 test. Specify path to mio Slater-Koster ' \
                  'tables in MIO environment variable if you want to run it.')
        else:
            run_dftb3_test(self)

###

if __name__ == '__main__':
    run_dftb3_test()
