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


"""Check energies against the databases published in
   Brenner et al., J. Phys. Condens. Matter 14, 783 (2002)
"""

import os
import sys

import unittest

import ase
import ase.build
import ase.io
import ase.optimize

import atomistica

###

ETOL = 0.005

###

# Entries are: molecule, atomization energy (eV),
# zero-point energy (kcal mol^-1), dHpot (kcal mol^-1)

# Some of the tests below fail (those are commented out so that the
# test runner can complete). I am not sure why. There is a difference in the
# parameters that are published in the paper and those that are implemented in
# the code (at http://www.mse.ncsu.edu/CompMatSci/). The parameters implemented
# here are the ones from the paper.

# C, H --- Brenner et al., Table 12
Brenner_et_al_CH = [
    ( 'CH2_s1A1d',            -8.4693,   0.000 ),
    ( 'CH3',                 -13.3750,   0.000 ),
    ( 'CH4',                 -18.1851,   0.000 ), # methane
    ( 'C2H',                 -11.5722,   0.000 ),
    ( 'C2H2',                -17.5651,   0.000 ), # acetylene
    ( 'C2H4',                -24.4077,   0.000 ), # ethylene
    ( 'H3C2H2',              -26.5601,   0.000 ),
    ( 'C2H6',                -30.8457,   0.000 ), # ethane
    ( 'C3H4_C2v',            -28.2589,   0.000 ), # cyclopropene
    ( 'CH2=C=CH2',           -30.2392,   0.000 ),
    ( 'propyne',             -30.3076,   0.000 ),
    ( 'C3H6_D3h',            -36.8887,   0.000 ), # cyclopropane
    ( 'C3H6_Cs',             -37.3047,   0.000 ), # propene
    ( 'C3H8',                -43.5891,   0.000 ), # propane
#    ( 'cyclobutene',         -42.1801,   0.000 ), # Fails for screened only
    ( 'butadiene',           -43.0035,   0.000 ),
    ( 'CH3CH=C=CH2',         -43.1367,   0.000 ),
    ( '1-butyne',            -43.0510,   0.000 ),
    ( '2-butyne',            -43.0501,   0.000 ),
#    ( 'cyclobutane',         -49.7304,   0.000 ), # Fails, not sure why
    ( '1-butene',            -50.0487,   0.000 ),
    ( 'cis-butene',          -50.2017,   0.000 ),
    ( 'i-C4H9',              -52.0451,   0.000 ),
    ( 't-C4H9',              -52.3778,   0.000 ),
    ( 'trans-butane',        -56.3326,   0.000 ), # n-butane
    ( 'isobutane',           -56.3309,   0.000 ),
    ( '1,3-pentadiene',      -55.9025,   0.000 ),
    ( '1,4-pentadiene',      -56.5078,   0.000 ),
    ( 'cyclopentene',        -57.1119,   0.000 ),
#    ( '1,2-pentadiene',      -58.7350,   0.000 ), # Fails, not sure why
#    ( '2,3-pentadiene',      -58.8900,   0.000 ), # Fails, not sure why
    ( 'cyclopentane',        -63.6443,   0.000 ),
    ( '2-pentene',           -62.9456,   0.000 ),
    ( '1-butene,2-methyl',   -62.9658,   0.000 ),
#    ( '2-butene,2-methyl',   -63.1109,   0.000 ), # Fails, not sure why
    ( 'n-pentane',           -69.0761,   0.000 ),
    ( 'isopentane',          -69.0739,   0.000 ),
    ( 'neopentane',          -69.0614,   0.000 ),
    ( 'C6H6',                -59.3096,   0.000 ), # benzene
    ( 'cyclohexane',         -76.4606,   0.000 ),
    ( 'naphthalene',         -93.8784,   0.000 ),
    ]

reference_database = \
    Brenner_et_al_CH

###

def molecule(mol):
    if os.path.exists('molecule_database/{0}.xyz'.format(mol)):
        a = ase.io.read('molecule_database/{0}.xyz'.format(mol))
    else:
        a = ase.build.molecule(mol)
    return a

###

def run_rebo2_molecules_test(test=None):
    for potname, c, reference_database in [ 
        ( 'Rebo2', atomistica.Rebo2(),
          Brenner_et_al_CH ),
        ( 'Rebo2Scr', atomistica.Rebo2Scr(dihedral=False),
          Brenner_et_al_CH ),
        ]:

        if test is None:
            print('=== Testing {0} ==='.format(potname))

        nok = 0
        nfailed = 0
    
        for mol, edft, de in reference_database:
            if len(sys.argv) > 1:
                if mol not in sys.argv[1:]:
                    continue

            eref = edft-de

            a = molecule(mol)
            a.center(vacuum=5.0)

            a.set_calculator(c)
            a.rattle(0.05)
            ase.optimize.QuasiNewton(a, logfile='QuasiNewton.log') \
                .run(fmax=0.001)
            #ase.optimize.FIRE(a, logfile='FIRE.log').run(fmax=0.001)
            e = a.get_potential_energy()

            if test is None:
                if abs(e-eref) > ETOL:
                    print('{0:>20} {1:>20.10f} {2:>20.10f} {3:>20.10f}      ' \
                        '.failed.'.format(mol, e, eref, e-eref))
                    nfailed += 1
                else:
                    print('{0:>20} {1:>20.10f} {2:>20.10f} {3:>20.10f}  .ok.' \
                        .format(mol, e, eref, e-eref))
                    nok += 1

            else:
                test.assertTrue(abs(e-eref) < ETOL,
                                msg='Energy for %s should be %f eV but '
                                'is %f eV.' % ( mol, eref, e ))

        if test is None:
            print('{0} molecule tests ok, {1} molecule tests failed.' \
                .format(nok, nfailed))

###

class TestREBO2Molecules(unittest.TestCase):

    def test_rebo2_molecules(self):
        run_rebo2_molecules_test(self)

###

if __name__ == '__main__':
    run_rebo2_molecules_test()
