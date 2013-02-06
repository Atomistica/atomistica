#! /usr/bin/env python

import sys

import unittest

import ase

from ase.lattice.cubic import Diamond, FaceCenteredCubic, SimpleCubic
from ase.lattice.cubic import BodyCenteredCubic
from ase.lattice.compounds import B1, B2, B3, L1_2

from mdcore import *
from mdcore.tests import test_forces, test_virial

# import ase_ext_io as io

###

sx         = 2
dx         = 1e-6
dev_thres  = 1e-4

###

tests  = [
    ( r6, dict(el1='Si', el2='Si', A=1.0, r0=1.0, cutoff=5.0),
      [ ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( Brenner,   Erhart_PRB_71_035211_SiC,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ),
        ( "dia-Si-C", B3( [ "Si", "C" ], latticeconstant=4.3596,
                          size=[sx,sx,sx]) ) ] ),
    ( BrennerScr, Erhart_PRB_71_035211_SiC__Scr,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ),
        ( "dia-Si-C", B3( [ "Si", "C" ], latticeconstant=4.3596,
                          size=[sx,sx,sx]) ) ] ),
    ( Brenner,   Henriksson_PRB_79_114107_FeC,
      [ dict( name='dia-C', struct=Diamond('C', size=[sx,sx,sx]) ),
        dict( name='bcc-Fe',
              struct=BodyCenteredCubic('Fe', size=[sx,sx,sx]) ),
        dict( name='fcc-Fe',
              struct=FaceCenteredCubic('Fe', size=[sx,sx,sx],
                                       latticeconstant=3.6) ),
        dict( name='sc-Fe',
              struct=SimpleCubic('Fe', size=[sx,sx,sx], latticeconstant=2.4) ),
        dict( name='B1-Fe-C',
              struct=B1( [ 'Fe', 'C' ], size=[sx,sx,sx], latticeconstant=3.9) ),
        dict( name='B3-Fe-C',
              struct=B3( [ 'Fe', 'C' ], size=[sx,sx,sx], latticeconstant=4.0) ),
        ] ),
    ( Kumagai,    Kumagai_CompMaterSci_39_457_Si,
      [ ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( KumagaiScr, Kumagai_CompMaterSci_39_457_Si__Scr,
      [ ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( Tersoff,    Tersoff_PRB_39_5566_Si_C,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ( TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]) ),
        ( "dia-Si", Diamond("Si", size=[sx,sx,sx]) ) ] ),
    ]

###

def run_forces_and_virial_test(test=None):
    for pot, par, mats in tests:
        if len(sys.argv) > 1:
            found = False
            if par is not None:
                for keyword in sys.argv[1:]:
                    if '__ref__' in par:
                        if par['__ref__'].lower().find(keyword.lower()) != -1:
                            found = True
            if not found:
                continue

        try:
            potname = pot.__name__
        except:
            potname = pot.__class__.__name__
        if test is None:
            print "--- %s ---" % potname
        if par is None:
            c  = pot()
        else:
            c  = pot(**par)
            if test is None and '__ref__' in par:
                print "    %s" % par["__ref__"]

        for imat in mats:
            if isinstance(imat, tuple):
                name, a = imat
            else:
                name = imat['name']
                a = imat['struct']
            if test is None:
                print "Material:  ", name
            a.translate([0.1,0.1,0.1])
            a.set_calculator(c)

            for dummy in range(2):
                if dummy == 0:
                    errmsg = 'potential: {0}; material: {1}; equilibrium' \
                        .format(potname, name)
                    if test is None:
                        print "...equilibrium..."
                else:
                    errmsg = 'potential: {0}; material: {1}; distorted' \
                        .format(potname, name)
                    if test is None:
                        print "...distorted..."

                ffd, f0, maxdf  = test_forces(a, dx=dx)
        
                if test is None:
                    if abs(maxdf) < dev_thres:
                        print "forces .ok."
                    else:
                        print "forces .failed."
                        print "max(df)  = %f" % maxdf

                        print "f - from potential"
                        print f0

                        print "f - numerically"
                        print ffd
                else:
                    test.assertTrue(abs(maxdf) < dev_thres,
                                    msg=errmsg+'; forces')

                sfd, s0, maxds  = test_virial(a, de=dx)

                if test is None:
                    if abs(maxds) < dev_thres:
                        print "virial .ok."
                    else:
                        print "virial .failed."
                        print "max(ds)  = %f" % maxds
                    
                        print "s - from potential"
                        print s0
                
                        print "s - numerically"
                        print sfd
                else:
                    test.assertTrue(abs(maxds) < dev_thres,
                                    msg=errmsg+'; virial')
            
                a.rattle(0.5)

###

class TestForcesAndVirial(unittest.TestCase):

    def test_forces_and_virial(self):
        run_forces_and_virial_test(self)

###

if __name__ == '__main__':
    run_forces_and_virial_test()
