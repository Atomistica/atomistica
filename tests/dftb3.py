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
# mio-1-1 Slater-Koster tables
# zeta = 4.05 for full dftb3 + XH
# zeta = 3.70 for full dftb2 + XH
#
# Table 4 from Gaus, Goez, Elstner, JCTC 9, 338 (2013)
# 3ob Slater-Koster tables
# zeta = 4.00 for full dftb3 + XH
table5_data = {
    '2H2O': {
        'reference_structures': ['H2O', 'H2O'],
        'structure': 'molecule_database/2H2O.xyz',
        'G3B3': -4.9, # kcal/mol
        'DFTB2': 1.6,
        'DFTB2+XH': 0.0,
        'DFTB3': 1.5,
        'DFTB3+XH': 0.0,
        'DFTB3+XH_3ob': 0.3
    },
    '3H2O': {
        'reference_structures': ['H2O', 'H2O', 'H2O'],
        'structure': 'molecule_database/3H2O.xyz',
        'G3B3': -15.1, # kcal/mol
        'DFTB2': 5.5,
        'DFTB2+XH': -0.6,
        'DFTB3': 5.4,
        'DFTB3+XH': -0.3,
        'DFTB3+XH_3ob': 0.8
    },
    '4H2O': {
        'reference_structures': ['H2O', 'H2O', 'H2O', 'H2O'],
        'structure': 'molecule_database/4H2O.xyz',
        'G3B3': -27.4, # kcal/mol
        'DFTB2': 9.7,
        'DFTB2+XH': 0.6,
        'DFTB3': 9.4,
        'DFTB3+XH': 0.8,
        'DFTB3+XH_3ob': 2.8
    },
    '5H2O': {
        'reference_structures': ['H2O', 'H2O', 'H2O', 'H2O', 'H2O'],
        'structure': 'molecule_database/5H2O.xyz',
        'G3B3': -36.3, # kcal/mol
        'DFTB2': 13.3,
        'DFTB2+XH': 1.4,
        'DFTB3': 12.5,
        'DFTB3+XH': 1.3,
        'DFTB3+XH_3ob': 3.7
    }
}

###

def run_dftb3_test(test=None):
    mio_database_folder = os.getenv('MIO')
    if mio_database_folder is None:
        raise RuntimeError('Please use environment variable MIO to specify path to mio Slater-Koster tables.')
    dftb3_database_folder = os.getenv('DFTB3')
    if dftb3_database_folder is None:
        raise RuntimeError('Please use environment variable 3OB to specify path to 3ob Slater-Koster tables.')

    print('mio folder: {}'.format(mio_database_folder))
    print('3ob folder: {}'.format(dftb3_database_folder))

    dftb2_calc = Atomistica(
        [ native.TightBinding(
        database_folder = mio_database_folder,
        SolverLAPACK = dict(electronic_T=0.001),
        SCC = dict(dq_crit = 1e-6,
                   mixing = 0.05, # 0.2
                   andersen_memory = 15, # 3
                   maximum_iterations = 100,
                   log = True)
        ),
          native.DirectCoulomb(),
          native.SlaterCharges(cutoff=10.0) ],
        avgn = 1000
        )

    dftb2_XH_calc = Atomistica(
        [ native.TightBinding(
        database_folder = mio_database_folder,
        SolverLAPACK = dict(electronic_T=0.001),
        SCC = dict(dq_crit = 1e-6,
                   mixing = 0.05, # 0.2
                   andersen_memory = 15, # 3
                   maximum_iterations = 100,
                   log = True)
        ),
          native.DirectCoulomb(),
          native.SlaterCharges(cutoff=10.0, damp_gamma=True, zeta = 3.70) ],
        avgn = 1000
        )

    dftb3_calc = Atomistica(
        [ native.TightBinding(
        database_folder = mio_database_folder,
        SolverLAPACK = dict(electronic_T=0.001),
        SCC = dict(dq_crit = 1e-6,
                   mixing = 0.05, # 0.2
                   andersen_memory = 15, # 3
                   maximum_iterations = 100,
                   log = True)
        ),
          native.DirectCoulomb(),
          native.SlaterCharges(cutoff=10.0, dftb3=True,  
                               HubbardDerivatives=dict(H=-0.1857, O=-0.1575)) ],
        avgn = 1000
        )

    dftb3_XH_calc = Atomistica(
        [ native.TightBinding(
        database_folder = mio_database_folder,
        SolverLAPACK = dict(electronic_T=0.001),
        SCC = dict(dq_crit = 1e-6,
                   mixing = 0.05, # 0.2
                   andersen_memory = 15, # 3
                   maximum_iterations = 100,
                   log = True)
        ),
          native.DirectCoulomb(),
          native.SlaterCharges(cutoff=10.0, dftb3=True, damp_gamma=True, zeta = 4.05, 
                               HubbardDerivatives=dict(H=-0.1857, O=-0.1575)) ],
        avgn = 1000
        )

    dftb3_XH_3ob_calc = Atomistica(
        [ native.TightBinding(
        database_folder = dftb3_database_folder,
        SolverLAPACK = dict(electronic_T=0.001),
        SCC = dict(dq_crit = 1e-6,
                   mixing = 0.05, # 0.2
                   andersen_memory = 15, # 3
                   maximum_iterations = 100,
                   log = True)
        ),
          native.DirectCoulomb(),
          native.SlaterCharges(cutoff=10.0, dftb3=True, damp_gamma=True, zeta = 4.00,
                               HubbardDerivatives=dict(H=-0.1857, O=-0.1575)) ],
        avgn = 1000
        )

    print('    nH2O       G3B3|          DFTB2 (MIO)|       DFTB2+XH (MIO)|          DFTB3 (MIO)|       DFTB3+XH (MIO)|       DFTB3+XH (3OB)|')

    for name, data in table5_data.items():
        e0_DFTB2        = 0.0
        e0_DFTB3        = 0.0
        e0_DFTB2_XH     = 0.0
        e0_DFTB3_XH     = 0.0
        e0_DFTB3_3ob_XH = 0.0
        for structure in data['reference_structures']:
            if os.path.exists(structure):
                a = read(structure)
            else:
                a = molecule(structure)

            a.center(vacuum=10.0)
            a.set_calculator(dftb2_calc)
            FIRE(a).run(fmax=0.001)
            e0_DFTB2 += a.get_potential_energy()

            a.set_calculator(dftb2_XH_calc)
            FIRE(a).run(fmax=0.001)
            e0_DFTB2_XH += a.get_potential_energy()

            a.set_calculator(dftb3_calc)
            FIRE(a).run(fmax=0.001)
            e0_DFTB3 += a.get_potential_energy()

            a.set_calculator(dftb3_XH_calc)
            FIRE(a).run(fmax=0.001)
            e0_DFTB3_XH += a.get_potential_energy()

            a.set_calculator(dftb3_XH_3ob_calc)
            FIRE(a).run(fmax=0.001)
            e0_DFTB3_3ob_XH += a.get_potential_energy()

        eref_G3B3 = data['G3B3']
        eref_DFTB2 = data['DFTB2']
        eref_DFTB2_XH = data['DFTB2+XH']
        eref_DFTB3 = data['DFTB3']
        eref_DFTB3_XH = data['DFTB3+XH']
        eref_DFTB3_XH_3ob = data['DFTB3+XH_3ob']

        a = read(data['structure'])
        a.center(vacuum=10.0)
        a.set_calculator(dftb2_calc)
        FIRE(a).run(fmax=0.001)
        e_DFTB2 = a.get_potential_energy()

        e_DFTB2 = (e_DFTB2 - e0_DFTB2)/(ase.units.kcal/ase.units.mol)

        a.set_calculator(dftb2_XH_calc)
        FIRE(a).run(fmax=0.001)
        e_DFTB2_XH = a.get_potential_energy()

        e_DFTB2_XH = (e_DFTB2_XH - e0_DFTB2_XH)/(ase.units.kcal/ase.units.mol)

        a.set_calculator(dftb3_calc)
        FIRE(a).run(fmax=0.001)
        e_DFTB3 = a.get_potential_energy()

        e_DFTB3 = (e_DFTB3 - e0_DFTB3)/(ase.units.kcal/ase.units.mol)

        a.set_calculator(dftb3_XH_calc)
        FIRE(a).run(fmax=0.001)
        e_DFTB3_XH = a.get_potential_energy()

        e_DFTB3_XH = (e_DFTB3_XH - e0_DFTB3_XH)/(ase.units.kcal/ase.units.mol)

        a.set_calculator(dftb3_XH_3ob_calc)
        FIRE(a).run(fmax=0.001)
        e_DFTB3_3ob_XH = a.get_potential_energy()

        e_DFTB3_3ob_XH = (e_DFTB3_3ob_XH - e0_DFTB3_3ob_XH)/(ase.units.kcal/ase.units.mol)

        if test is None:
            print('{0:>8} {1:>10.5f} {2:>10.5f} {3:>10.5f} {4:>10.5f} {5:>10.5f} {6:>10.5f} {7:>10.5f} {8:>10.5f} {9:>10.5f} {10:>10.5f} {11:>10.5f}'.format(name, eref_G3B3, e_DFTB2 - eref_G3B3, eref_DFTB2, e_DFTB2_XH - eref_G3B3, eref_DFTB2_XH, e_DFTB3 - eref_G3B3, eref_DFTB3, e_DFTB3_XH - eref_G3B3, eref_DFTB3_XH, e_DFTB3_3ob_XH - eref_G3B3, eref_DFTB3_XH_3ob))
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
