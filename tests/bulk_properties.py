#! /usr/bin/env python

#
# Test the bulk properties for a set of potentials
#

import sys

from math import sqrt

import unittest

import numpy as np

import ase
import ase.constraints

from atomistica import *
from atomistica.tests import test_cubic_elastic_constants

from ase.lattice.cubic import Diamond, BodyCenteredCubic
from ase.lattice.cubic import FaceCenteredCubic, SimpleCubic
from ase.lattice.compounds import B1, B2, B3, L1_2

###

sx        = 1
dev_thres = 5

k0 = ase.units.GPa

tests  = [
    ( Brenner,   Brenner_PRB_42_9458_C_I,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]),
          # Ec    a0     C11   C12   C44   B    Cp
          None,   None,  None, None, None, None, None ),
        ] ),
    ( Brenner,   Brenner_PRB_42_9458_C_II,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]),
          # Ec          a0     C11   C12  C44  B    Cp
          7.376-0.0524, 3.558, 621,  415, 383, 484, None ),
        ] ),
    ( Brenner,   Erhart_PRB_71_035211_SiC,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]),
          # Ec    a0     C11   C12  C44  B    Cp
          7.3731, 3.566, 1082, 127, 673, 445, None ),
        dict( name="dia-Si", struct=Diamond("Si", size=[sx,sx,sx]),
              Ec=4.63, a0=5.429, C11=167, C12=65, C440=105, B=99 ),
        dict( name="dia-Si-C", struct=B3( [ "Si", "C" ], latticeconstant=4.3596,
                                          size=[sx,sx,sx]),
              Ec=6.340,a0=4.359, C11=382, C12=145, C440=305, B=224 ) ] ),
    ( BrennerScr,   Erhart_PRB_71_035211_SiC__Scr,
      [ ( "dia-C", Diamond("C", size=[sx,sx,sx]),
          # Ec    a0     C11   C12  C44  B    Cp
          7.3731, 3.566, 1082, 127, 673, 445, None ),
        dict( name="dia-Si", struct=Diamond("Si", size=[sx,sx,sx]),
              Ec=4.63, a0=5.429, C11=167, C12=65, C440=105, B=99 ),
        dict( name="dia-Si-C", struct=B3( [ "Si", "C" ], latticeconstant=4.3596,
                                          size=[sx,sx,sx]),
              Ec=6.340,a0=4.359, C11=382, C12=145, C440=305, B=224 ) ] ),
    ( Kumagai, Kumagai_CompMaterSci_39_457_Si,
      [ dict( name="dia-Si", struct=Diamond("Si", size=[sx,sx,sx]),
              Ec=4.630, a0=5.429, C11=166.4, C12=65.3, C440=120.9 ),
        ] ),
    ( KumagaiScr, Kumagai_CompMaterSci_39_457_Si__Scr,
      [ dict( name="dia-Si", struct=Diamond("Si", size=[sx,sx,sx]),
              Ec=4.630, a0=5.429, C11=166.4, C12=65.3, C440=120.9 ),
        ] ),
    ( Tersoff,    Tersoff_PRB_39_5566_Si_C,
      [ dict( name="dia-C", struct=Diamond("C", size=[sx,sx,sx]),
              Ec=7.396-0.0250, a0=3.566, C11=1067, C12=104, C44=636,
              C440=671 ),
        dict( name="dia-Si", struct=Diamond("Si", size=[sx,sx,sx]),
              Ec=4.63, a0=5.432, C11=143, C12=75, C44=69, C440=119, B=98 ),
        dict( name="dia-Si-C", struct=B3( [ "Si", "C" ], latticeconstant=4.3596,
                                          size=[sx,sx,sx]),
              Ec=6.165, a0=4.321, C11=437, C12=118, C440=311, B=224 ),
        ] ),
    ( TersoffScr,    Tersoff_PRB_39_5566_Si_C__Scr,
      [ dict( name="dia-C", struct=Diamond("C", size=[sx,sx,sx]),
              Ec=7.396-0.0250, a0=3.566, C11=1067, C12=104, C44=636,
              C440=671 ),
        dict( name="dia-Si", struct=Diamond("Si", size=[sx,sx,sx]),
              Ec=4.63, a0=5.432, C11=143, C12=75, C44=69, C440=119, B=98 ),
        dict( name="dia-Si-C", struct=B3( [ "Si", "C" ], latticeconstant=4.3596,
                                          size=[sx,sx,sx]),
              Ec=6.165, a0=4.321, C11=437, C12=118, C440=311, B=224 ),
        ] ),
    ]

###

class TestElasticConstants(unittest.TestCase):

    def test_elastic_constants(self):
        for pot, par, mats in tests:
             test_cubic_elastic_constants(mats, pot, par, sx, dev_thres,
                                         test=self)

###

if __name__ == '__main__':
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

        test_cubic_elastic_constants(mats, pot, par, sx, dev_thres)
