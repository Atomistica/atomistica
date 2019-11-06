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
Test the bulk properties for a set of potentials
"""

from __future__ import print_function

import sys

from math import sqrt

import unittest

import numpy as np

import ase
import ase.constraints
from ase.units import GPa

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
    ( Harmonic, dict(k=1.0, r0=1.0, cutoff=1.3, shift=True),
      [ dict( name="fcc", struct=FaceCenteredCubic("He", size=[sx,sx,sx],
                                                   latticeconstant=sqrt(2)),
              C11=sqrt(2)/GPa, C12=1./sqrt(2)/GPa, C44=1./sqrt(2)/GPa )
        ] ),
    ( DoubleHarmonic, dict(k1=1.0, r1=1.0, k2=1.0, r2=sqrt(2), cutoff=1.6),
      [ dict( name="sc", struct=SimpleCubic("He", size=[sx,sx,sx],
                                            latticeconstant=1.0),
              C11=3./GPa, C12=1./GPa, C44=1./GPa )
        ] ),
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
    ( Juslin, Juslin_JAP_98_123520_WCH,
      [ dict( name= "bcc-W", struct=BodyCenteredCubic("W", size=[sx,sx,sx]),
              Ec=8.89, a0=3.165, C11=542, C12=191, C44=162, B=308 ),
        dict( name="fcc-W", struct=FaceCenteredCubic("W", latticeconstant=4.0,
                                                     size=[sx,sx,sx]),
              Ec=8.89-0.346, a0=4.005 ),
        # Note: The sc test uses 2 2x2x2 unit cell. The lattice constant is
        # therefore twice what is listed in Juslin's paper.
        dict( name="sc-W", struct=SimpleCubic("W", latticeconstant=2.7,
                                              size=[2*sx,2*sx,2*sx]),
               Ec=8.89-1.614, a0=2*2.671 ),
        dict( name="dia-C", struct=Diamond("C", size=[sx,sx,sx]),
              Ec=7.376-0.0524, a0=3.558, C11=621, C12=415, C44=383, B=484 ),
        dict( name='B1-W-C', struct=B1( [ 'W', 'C' ], latticeconstant=4.38,
                                        size=[sx,sx,sx]),
              Ec=(16.68-0.98)/2, a0=4.380, B=433 ),
        dict( name='B2-W-C', struct=B2( [ 'W', 'C' ], latticeconstant=2.7,
                                        size=[sx,sx,sx]),
              Ec=(16.68-2.32)/2, a0=2.704, B=411 ),
        dict( name='B3-W-C', struct=B3( [ 'W', 'C' ], latticeconstant=4.6,
                                        size=[sx,sx,sx]),
              Ec=(16.68-2.12)/2, a0=4.679, B=511 ),
        ] ),
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
    ( TabulatedEAM, dict(fn='Au_u3.eam'),
       [ dict( name='fcc-Au', struct=FaceCenteredCubic('Au', size=[sx,sx,sx]),
               Ec=3.93, a0=4.08, B=167, C11=183, C12=159, C44=45)
         ] ),
    ( TabulatedAlloyEAM, dict(fn='Au-Grochola-JCP05.eam.alloy'),
       [ dict( name='fcc-Au', struct=FaceCenteredCubic('Au', size=[sx,sx,sx]),
               Ec=3.924, a0=4.070, C11=202, C12=170, C44=47, C440=46)
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

        _nok, _nfail = test_cubic_elastic_constants(mats, pot, par, sx, dev_thres)
        nok += _nok
        nfail += _nfail
    print('{0} tests passed, {1} tests failed.'.format(nok, nfail))
