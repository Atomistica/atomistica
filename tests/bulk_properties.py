#! /usr/bin/env python

# ======================================================================
# Atomistica - Interatomic potential library
# https://github.com/pastewka/atomistica
# Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
# See the AUTHORS file in the top-level Atomistica directory.
#
# Copyright (2005-2013) Fraunhofer IWM
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
    ( Rebo2(),   None,
      [ dict( name='dia-C', struct=Diamond('C', size=[sx,sx,sx]),
              Ec=7.370, a0=3.566, C11=1080, C12=130, C44=720 )
        ] ),
    ( Rebo2Scr(),   None,
      [ dict( name='dia-C', struct=Diamond('C', size=[sx,sx,sx]),
              Ec=7.370, a0=3.566, C11=1080, C12=130, C44=720 )
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
    ( Tersoff, Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N,
      [ dict( name="dia-C", struct=Diamond("C", size=[sx,sx,sx]),
              Ec=7.396-0.0250, a0=3.566, C11=1067, C12=104, C44=636,
              C440=671 ),
        dict( name="dia-B-N", struct=B3( [ "B", "N" ], latticeconstant=3.7,
                                         size=[sx,sx,sx]),
              Ec=6.63, a0=3.658, B=385 ),
        ] ),
    ( TersoffScr, Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr,
      [ dict( name="dia-C", struct=Diamond("C", size=[sx,sx,sx]),
              Ec=7.396-0.0250, a0=3.566, C11=1067, C12=104, C44=636,
              C440=671 ),
        dict( name="dia-B-N", struct=B3( [ "B", "N" ], latticeconstant=3.7,
                                         size=[sx,sx,sx]),
              Ec=6.63, a0=3.658, B=385 ),
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
            try:
                potname = pot.__name__
            except:
                potname = pot.__class__.__name__
            for keyword in sys.argv[1:]:
                if potname.lower().find(keyword.lower()) != -1:
                    found = True
            if not found:
                continue

        test_cubic_elastic_constants(mats, pot, par, sx, dev_thres)
