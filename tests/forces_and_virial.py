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


from __future__ import print_function

import math
import sys

import unittest

from numpy.random import randint

import ase
import ase.io as io
from ase.units import mol

from ase.lattice.cubic import Diamond, FaceCenteredCubic, SimpleCubic
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.compounds import B1, B2, B3, L1_2, NaCl

import atomistica.native as native
from atomistica import *
from atomistica.tests import test_forces, test_potential, test_virial

# import ase_ext_io as io

###

sx = 2
dx = 1e-6
tol = 1e-2

###

def random_solid(els, density):
    syms = [ ]
    nat = 0
    for sym, n in els:
        syms += n*[sym]
        nat += n
    r = np.random.rand(nat, 3)
    a = ase.Atoms(syms, positions=r, cell=[1,1,1], pbc=True)

    mass = np.sum(a.get_masses())
    a0 = ( 1e24*mass/(density*mol) )**(1./3)
    a.set_cell([a0,a0,a0], scale_atoms=True)

    return a

def assign_charges(a, els):
    syms = np.array(a.get_chemical_symbols())
    qs = np.zeros(len(a))
    for el, q in els.items():
        qs[syms==el] = q
    if hasattr(a, 'set_initial_charges'):
        a.set_initial_charges(qs)
    else:
        a.set_charges(qs)
    return a

###

# Potential tests
tests  = [
    ( Harmonic, dict(el1='He', el2='He', k=1.0, r0=1.0, cutoff=1.5),
      [ ( "fcc-He", FaceCenteredCubic("He", size=[sx,sx,sx],
                                      latticeconstant=math.sqrt(2.0)) ) ] ),
    ( r6, dict(el1='Si', el2='Si', A=1.0, r0=1.0, cutoff=5.0),
      [ ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( LJCut, dict(el1='He', el2='He', epsilon=10.2, sigma=2.28, cutoff=5.0,
                  shift=True),
      [ dict( name="fcc-He", struct=FaceCenteredCubic("He", size=[sx,sx,sx],
                                                      latticeconstant=3.5),
              mask=True, rattle=0.1 ) ] ),
    ( Brenner, Erhart_PRB_71_035211_SiC,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "a-C", io.read("aC_small.cfg") ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ),
        ( "dia-Si-C", B3( [ "Si", "C" ], latticeconstant=4.3596,
                          size=[sx,sx,sx]) ) ] ),
    ( BrennerScr, Erhart_PRB_71_035211_SiC__Scr,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "a-C", io.read("aC_small.cfg") ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ),
        ( "dia-Si-C", B3( [ "Si", "C" ], latticeconstant=4.3596,
                          size=[sx,sx,sx]) ) ] ),
    ( Brenner, Henriksson_PRB_79_114107_FeC,
      [ dict( name='dia-C', struct=Diamond('C', size=[sx,sx,sx]), mask=True ),
        dict( name="a-C", struct=io.read("aC_small.cfg"), mask=True ),
        dict( name='bcc-Fe',
              struct=BodyCenteredCubic('Fe', size=[sx,sx,sx]), mask=True ),
        dict( name='fcc-Fe',
              struct=FaceCenteredCubic('Fe', size=[sx,sx,sx],
                                       latticeconstant=3.6), mask=True ),
        dict( name='sc-Fe',
              struct=SimpleCubic('Fe', size=[sx,sx,sx], latticeconstant=2.4),
              mask=True ),
        dict( name='B1-Fe-C',
              struct=B1( [ 'Fe', 'C' ], size=[sx,sx,sx], latticeconstant=3.9),
              mask=True ),
        dict( name='B3-Fe-C',
              struct=B3( [ 'Fe', 'C' ], size=[sx,sx,sx], latticeconstant=4.0),
              mask=True ),
        ] ),
    ( Kumagai, Kumagai_CompMaterSci_39_457_Si,
      [ ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( KumagaiScr, Kumagai_CompMaterSci_39_457_Si__Scr,
      [ ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( Tersoff, Tersoff_PRB_39_5566_Si_C,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "a-C", io.read("aC_small.cfg") ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ),
        ( "dia-Si-C", B3( [ "Si", "C" ], latticeconstant=4.3596,
                          size=[sx,sx,sx]) ) ] ),
    ( TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "a-C", io.read("aC_small.cfg") ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ),
        ( "dia-Si-C", B3( [ "Si", "C" ], latticeconstant=4.3596,
                          size=[sx,sx,sx]) ) ] ),
    ( Rebo2, None,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "a-C", io.read("aC_small.cfg") ),
        ( 'random-C-H', random_solid( [('C',50),('H',10)], 3.0 ) ),
        ] ),
    ( Rebo2Scr, None,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "a-C", io.read("aC_small.cfg") ),
        ( 'random-C-H', random_solid( [('C',50),('H',10)], 3.0 ) ),
        ] ),
    ( TabulatedEAM, dict(fn='Au_u3.eam'),
      [ dict( name="fcc-Au", struct=FaceCenteredCubic("Au", size=[sx,sx,sx]),
              rattle=0.1 ) ] ),
    ( TabulatedAlloyEAM, dict(fn='Au-Grochola-JCP05.eam.alloy'),
      [ dict( name="fcc-Au", struct=FaceCenteredCubic("Au", size=[sx,sx,sx]),
              rattle=0.1, mask=True ) ] ),
    ]

# Coulomb potential tests
tests += [
    ( DirectCoulomb, None,
      [ ( "sc-Na-Cl", assign_charges(NaCl(['Na','Cl'], latticeconstant=5.64,
                                          size=[sx,sx,sx]),
                                     dict(Na=1,Cl=-1)) ),
        ( "random-Na-Cl", assign_charges(random_solid([('Na',50),('Cl',50)],
                                                      2.16),
                                         dict(Na=1,Cl=-1)) ),
        ] ),
    ( PME, dict(cutoff=5.0, grid=(10, 10, 10)),
      [ ( "sc-Na-Cl", assign_charges(NaCl(['Na','Cl'], latticeconstant=5.64,
                                          size=[sx,sx,sx]),
                                     dict(Na=1,Cl=-1)) ),
        ( "random-Na-Cl", assign_charges(random_solid([('Na',50),('Cl',50)],
                                                      2.16),
                                         dict(Na=1,Cl=-1)) ),
        ] ),
    ( SlaterCharges, dict(el=['Na','Cl'], U=[1.0,0.5], Z=[0.1,-0.2],
                          cutoff=5.0),
      [ ( "sc-Na-Cl", assign_charges(NaCl(['Na','Cl'], latticeconstant=5.64,
                                          size=[sx,sx,sx]),
                                     dict(Na=1,Cl=-1)) ),
        ( "random-Na-Cl", assign_charges(random_solid([('Na',50),('Cl',50)],
                                                      2.16),
                                         dict(Na=1,Cl=-1)) ),
        ] ),
    # Also test U1==U2
    ( SlaterCharges, dict(el=['Na','Cl'], U=[1.0,1.0], Z=[0.1,-0.2],
                          cutoff=5.0),
      [ ( "sc-Na-Cl", assign_charges(NaCl(['Na','Cl'], latticeconstant=5.64,
                                          size=[sx,sx,sx]),
                                     dict(Na=1,Cl=-1)) ),
        ( "random-Na-Cl", assign_charges(random_solid([('Na',50),('Cl',50)],
                                                      2.16),
                                         dict(Na=1,Cl=-1)) ),
        ] ),
    ( GaussianCharges, dict(el=['Na','Cl'], U=[1.0,0.5],
                            cutoff=5.0),
      [ ( "sc-Na-Cl", assign_charges(NaCl(['Na','Cl'], latticeconstant=5.64,
                                          size=[sx,sx,sx]),
                                     dict(Na=1,Cl=-1)) ),
        ( "random-Na-Cl", assign_charges(random_solid([('Na',50),('Cl',50)],
                                                      2.16),
                                         dict(Na=1,Cl=-1)) ),
        ] ),
    ]

###

def run_forces_and_virial_test(test=None):
    nok = 0
    nfail = 0
    for pot, par, mats in tests:
        if len(sys.argv) > 1:
            found = False
            if par is not None:
                for keyword in sys.argv[1:]:
                    if '__ref__' in par:
                        if par['__ref__'].lower().find(keyword.lower()) != -1:
                            found = True
            try:
                potname = pot.__name__
            except:
                potname = pot.__class__.__name__
            for keyword in sys.argv[1:]:
                if potname.lower().find(keyword.lower()) != -1:
                    found = True
            if not found:
                continue

        try:
            potname = pot.__name__
        except:
            potname = pot.__class__.__name__
        if test is None:
            print("--- %s ---" % potname)
        if par is None:
            c  = pot()
        else:
            c  = pot(**par)
            if test is None and '__ref__' in par:
                print("    %s" % par["__ref__"])

        for imat in mats:
            rattle = 0.5
            mask = False
            if isinstance(imat, tuple):
                name, a = imat
            else:
                name = imat['name']
                a = imat['struct']
                if 'rattle' in imat:
                    rattle = imat['rattle']
                if 'mask' in imat:
                    mask = imat['mask']
            if test is None:
                print("Material:  ", name)
            a.translate([0.1,0.1,0.1])
            a.set_calculator(c)

            masks = [None]
            if mask:
                masks += [randint(0, len(a), size=len(a)) < len(a)/2,
                          randint(0, len(a), size=len(a)) < len(a)/4]

            for dummy in range(2):
                if dummy == 0:
                    errmsg = 'potential: {0}; material: {1}; equilibrium' \
                        .format(potname, name)
                    if test is None:
                        print('=== equilibrium ===')
                else:
                    errmsg = 'potential: {0}; material: {1}; distorted' \
                        .format(potname, name)
                    if test is None:
                        print('=== distorted ===')

                for mask in masks:
                    if test is None and mask is not None:
                        print('--- using random mask ---')
                    c.set_mask(mask)

                    ffd, f0, maxdf = test_forces(a, dx=dx)
        
                    if test is None:
                        if abs(maxdf) < tol:
                            nok += 1
                            print("forces .ok.")
                        else:
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
                    else:
                      test.assertTrue(abs(maxdf) < tol,
                                        msg=errmsg+'; forces')

                    sfd, s0, maxds = test_virial(a, de=dx)

                    if test is None:
                        if abs(maxds) < tol:
                            nok += 1
                            print("virial .ok.")
                        else:
                            nfail += 1
                            print("virial .failed.")
                            print("max(ds)  = %f" % maxds)
                    
                            print("s - from potential")
                            print(s0)
                
                            print("s - numerically")
                            print(sfd)

                            print("difference between the above")
                            print(s0-sfd)
                    else:
                        test.assertTrue(abs(maxds) < tol,
                                        msg=errmsg+'; virial')

                    pfd, p0, maxdp = test_potential(a, dq=dx)

                    if test is None:
                        if abs(maxdp) < tol:
                            nok += 1
                            print("potential .ok.")
                        else:
                            nfail += 1
                            print("potential .failed.")
                            print("max(dp)  = %f" % maxdp)
                    
                            print("p - from potential")
                            print(p0)
                
                            print("p - numerically")
                            print(pfd)

                            print("difference between the above")
                            print(p0-pfd)
                    else:
                        test.assertTrue(abs(maxds) < tol,
                                        msg=errmsg+'; virial')
            
                a.rattle(rattle)
    if test is None:
        print('{0} tests passed, {1} tests failed.'.format(nok, nfail))

###

class TestForcesAndVirial(unittest.TestCase):

    def test_forces_and_virial(self):
        run_forces_and_virial_test(self)

###

if __name__ == '__main__':
    run_forces_and_virial_test()
