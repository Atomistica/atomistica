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
Test the mio parametrization of Frauenheim and co-workers.
"""

from __future__ import print_function

import os
import sys
import unittest

import numpy as np

from ase.build import molecule
from ase.optimize import FIRE

import atomistica.native as native
from atomistica import Atomistica

###

# From Elstner et al., Phys. Rev. B 58, 7260
noscc_db = {
    'C=O': 1.296,
    'C-N': 1.296,
    'N-H': 1.003,
    'C-H': 1.130,
    'OCN': 127.0
    }

scc_db = {
    'C=O': 1.224,
    'C-N': 1.382,
    'N-H': 0.996,
    'C-H': 1.131,
    'OCN': 125.5
    }

db1 = {
    False: noscc_db,
    True: scc_db
    }

# From Kruger et al., J. Chem. Phys. 122, 114110
db2 = {
    'H2': {
        'H-H': ( ( 0, 1 ), 0.750 )
        },
    'C2H2': {
        'C-H': ( ( 1, 2 ), 1.075 ),
        'C-C': ( ( 0, 1 ), 1.203 )
        },
    'C2H4': {
        'C-H': ( ( 0, 2 ), 1.094 ),
        'C-C': ( ( 0, 1 ), 1.328 )
        },
    'C2H6': {
        'C-H': ( ( 0, 3 ), 1.098 ),
        'C-C': ( ( 0, 1 ), 1.501 )
        },
    'HCN': {
        'C-H': ( ( 0, 2 ), 1.078 ),
        'C-N': ( ( 0, 1 ), 1.141 )
        },
    'NH3': {
        'N-H': ( ( 0, 1 ), 1.021 )
        },
    'CH4': {
        'C-H': ( ( 0, 1 ), 1.089 )
        },
    'CO': {
        # This differs from the paper, but I believe it's a typo
        # paper says: 1.200
        'C-O': ( ( 0, 1 ), 1.100 )
        },
    'H2CO': {
        'C-H': ( ( 1, 2 ), 1.143 ),
        'C-O': ( ( 0, 1 ), 1.183 )
        },
    'CH3OH': {
        'O-H': ( ( 1, 3 ), 0.980 ),
        'C-O': ( ( 0, 1 ), 1.422 )
        },
    'H2O': {
        'O-H': ( ( 0, 1 ), 0.968 )
        },
    'N2': {
        # This differs from the paper, but I believe it's a typo
        # paper says: 1.200
        'N-N': ( ( 0, 1 ), 1.113 )
        },
    'N2H4': {
        'N-H': ( ( 0, 2 ), 1.037 ),
        # This differs from the paper, and I don't know why
        # paper says: 1.442
        'N-N': ( ( 0, 1 ), 1.407 )
        },
    'H2O2': {
        'O-H': ( ( 0, 2 ), 0.991 ),
        'O-O': ( ( 0, 1 ), 1.453 )
        },
    'CO2': {
        'C-O': ( ( 0, 1 ), 1.165 )
        }
    }


def check_db(c, db, test=None):
    if test is None:
        print("%10s %10s %10s ( %10s )" \
            % ( "bond", "value", "reference", "error" ))
        print("%10s %10s %10s ( %10s )" \
            % ( "----", "-----", "---------", "-----" ))
    for mol, values in db.items():
        #if mol == 'H2O':
        if 1:
            a = molecule(mol)
            a.center(vacuum=10.0)
            a.set_pbc(False)
            a.set_initial_charges(np.zeros(len(a)))

            a.set_calculator(c)
            FIRE(a, logfile=None).run(fmax=0.001)

            for name, ( ( i1, i2 ), refvalue ) in values.items():
                value = a.get_distance(i1, i2)
                if test is None:
                    print('%10s %10.3f %10.3f ( %10.3f )' % \
                        ( name, value, refvalue, abs(value-refvalue) ))
                else:
                    test.assertTrue(abs(value-refvalue) < 0.01)

###

def run_mio_test(test=None):
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
    check_db(calc, db2, test=test)
    

###

class TestMIO(unittest.TestCase):

    def test_mio(self):
        if os.getenv('MIO') is None:
            print('Skipping MIO test. Specify path to mio Slater-Koster ' \
                  'tables in MIO environment variable if you want to run it.')
        else:
            run_mio_test(self)

###

if __name__ == '__main__':
    run_mio_test()
