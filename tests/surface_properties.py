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
Test the surface properties for a set of potentials
"""

import sys

from math import sqrt

import unittest

import numpy as np

import ase
import ase.io
try:
    from ase.symbols import string2symbols
except:
    from ase.atoms import string2symbols

from atomistica import *
from atomistica.tests import test_surface_energies

from ase.lattice.cubic import Diamond, FaceCenteredCubic
from ase.lattice.compounds import B1, B2, B3, L1_2

###

nx = 1
nz = 4

###

def dia_111(sym, a0):
    sym = string2symbols(sym)
    if len(sym) == 1:
        a = Diamond(sym[0],
                    size             = [nx, nx, nz],
                    latticeconstant  = a0,
                    directions=[ [1,-1,0], [1,1,-2], [1,1,1] ]
                    )
    else:
        a = B3(sym,
               size             = [nx, nx, nz],
               latticeconstant  = a0,
               directions=[ [1,-1,0], [1,1,-2], [1,1,1] ]
               )
    sx, sy, sz = a.get_cell().diagonal()
    a.translate([sx/(12*nx), sy/(4*nx), sz/(12*nz)])
    return a


def dia_111_glide(sym, a0):
    sym = string2symbols(sym)
    if len(sym) == 1:
        a = Diamond(sym[0],
                    size             = [nx, nx, nz],
                    latticeconstant  = a0,
                    directions=[ [1,-1,0], [1,1,-2], [1,1,1] ]
                    )
    else:
        a = B3(sym,
               size             = [nx, nx, nz],
               latticeconstant  = a0,
               directions=[ [1,-1,0], [1,1,-2], [1,1,1] ]
               )
    sx, sy, sz = a.get_cell().diagonal()
    a.translate([sx/(12*nx), sy/(4*nx), sz/(12*nz)-0.2])
    a.set_scaled_positions(a.get_scaled_positions()%1.0)
    return a


def dia_111_pandey(sym, a0, nx=nx, ny=nx, nz=nz):
    """2x1 Pandey reconstructed (111) surface."""
    sym = string2symbols(sym)
    if len(sym) == 1:
        a = Diamond(sym[0],
                    size             = [nx, ny, nz],
                    latticeconstant  = a0,
                    directions=[ [1,-1,0], [1,1,-2], [1,1,1] ]
                    )
    else:
        a = B3(sym,
               size             = [nx, ny, nz],
               latticeconstant  = a0,
               directions=[ [1,-1,0], [1,1,-2], [1,1,1] ]
               )
    sx, sy, sz = a.get_cell().diagonal()
    a.translate([sx/(12*nx), sy/(4*ny), sz/(6*nz)])
    a.set_scaled_positions(a.get_scaled_positions()%1.0)

    bulk = a.copy()

    bondlen = a0*sqrt(3)/4

    x, y, z = a.positions.T
    mask = np.abs(z-z.max()) < 0.1*a0
    top1, top2 = np.arange(len(a))[mask].reshape(-1, 2).T
    mask = np.logical_and(np.abs(z-z.max()) < bondlen, np.logical_not(mask))
    topA, topB = np.arange(len(a))[mask].reshape(-1, 2).T
    y[topA] += bondlen/3
    y[topB] -= bondlen/3
    y[top1] += bondlen
    x[top1] += a.cell[0,0]/(2*nx)
    x[top2] += a.cell[0,0]/(2*nx)

    mask = np.abs(z-z.min()) < 0.1*a0
    bot1, bot2 = np.arange(len(a))[mask].reshape(-1, 2).T
    mask = np.logical_and(np.abs(z-z.min()) < bondlen, np.logical_not(mask))
    botA, botB = np.arange(len(a))[mask].reshape(-1, 2).T
    y[botA] += bondlen/3
    y[botB] -= bondlen/3
    y[bot2] -= bondlen
    x[bot2] += a.cell[0,0]/(2*nx)
    x[bot1] += a.cell[0,0]/(2*nx)

    a.set_scaled_positions(a.get_scaled_positions()%1.0)

    return bulk, a


def dia_110(sym, a0):
    sym = string2symbols(sym)
    if len(sym) == 1:
        a = Diamond(sym[0],
                    size             = [nx, nx, nz],
                    latticeconstant  = a0,
                    directions=[ [0,0,1], [1,-1,0], [1,1,0] ]
                    )
    else:
        a = B3(sym,
               size             = [nx, nx, nz],
               latticeconstant  = a0,
               directions=[ [0,0,1], [1,-1,0], [1,1,0] ]
               )
    sx, sy, sz = a.get_cell().diagonal()
    a.translate([sx/(4*nx), sy/(4*nx), sz/(8*nz)])
    return a


def dia_100(sym, a0):
    sym = string2symbols(sym)
    if len(sym) == 1:
        a = Diamond(sym[0],
                    size             = [nx, nx, nz],
                    latticeconstant  = a0,
                    directions=[ [1,0,0], [0,1,0], [0,0,1] ]
                    )
    else:
        a = B3(sym,
               size             = [nx, nx, nz],
               latticeconstant  = a0,
               directions=[ [1,0,0], [0,1,0], [0,0,1] ]
               )
    sx, sy, sz = a.get_cell().diagonal()
    a.translate([sx/(8*nx), sy/(8*nx), sz/(8*nz)])
    return a


def dia_100_2x1(sym, a0):
    sym = string2symbols(sym)
    if len(sym) == 1:
        a = Diamond(sym[0],
                    size             = [2*nx, nx, nz],
                    latticeconstant  = a0,
                    directions=[ [1,-1,0], [1,1,0], [0,0,1] ]
                    )
    else:
        a = B3(sym,
               size             = [2*nx, nx, nz],
               latticeconstant  = a0,
               directions=[ [1,-1,0], [1,1,0], [0,0,1] ]
               )
    sx, sy, sz = a.get_cell().diagonal()
    a.translate([sx/(8*nx), sy/(8*nx), sz/(8*nz)])

    bulk = a.copy()
    
    for i in a:
        if i.z < sz/(4*nz) or i.z > sz-sz/(4*nz):
            if i.x < sx/2:
                i.x = i.x+0.5
            else:
                i.x = i.x-0.5

    return bulk, a

###

vacuum = 10.0

tests = [
    ( Brenner, Erhart_PRB_71_035211_SiC,
      [ dict( name="dia-C-111", struct=dia_111('C', 3.566), r_Jm2=2.06 ),
        dict( name="dia-C-111-pandey", struct=dia_111_pandey('C', 3.566) ),
        dict( name="dia-C-110", struct=dia_110('C', 3.566), r_Jm2=2.96 ),
        dict( name="dia-C-100", struct=dia_100('C', 3.566), r_Jm2=5.59 ),
        dict( name="dia-C-100-2x1", struct=dia_100_2x1('C', 3.566),
              r_Jm2=5.65 ), 
        dict( name="dia-Si-111", struct=dia_111('Si', 5.432), r_Jm2=0.999 ),
        dict( name="dia-Si-111-pandey", struct=dia_111_pandey('Si', 5.432) ),
        dict( name="dia-Si-110", struct=dia_110('Si', 5.432), r_Jm2=1.23 ),
        dict( name="dia-Si-100", struct=dia_100('Si', 5.432), r_Jm2=1.95 ),
        dict( name="dia-Si-100-2x1", struct=dia_100_2x1('Si', 5.432),
              r_Jm2=1.13 ),
        dict( name="dia-SiC-111", struct=dia_111('SiC', 4.321), r_Jm2=1.67 ),
        dict( name="dia-SiC-110", struct=dia_110('SiC', 4.321), r_Jm2=2.29 ),
        dict( name="dia-SiC-100", struct=dia_100('SiC', 4.321), r_Jm2=3.93 ),
        dict( name="dia-SiC-100-2x1", struct=dia_100_2x1('SiC', 4.321),
              r_Jm2=2.85 )
        ] ),
    ( BrennerScr, Erhart_PRB_71_035211_SiC__Scr,
      [ dict( name="dia-C-111", struct=dia_111('C', 3.566), r_Jm2=2.06 ),
        dict( name="dia-C-111-pandey", struct=dia_111_pandey('C', 3.566) ),
        dict( name="dia-C-110", struct=dia_110('C', 3.566), r_Jm2=2.96 ),
        dict( name="dia-C-100", struct=dia_100('C', 3.566), r_Jm2=5.88 ),
        dict( name="dia-C-100-2x1", struct=dia_100_2x1('C', 3.566),
              r_Jm2=5.89 ), 
        dict( name="dia-Si-111", struct=dia_111('Si', 5.432), r_Jm2=0.999 ),
        dict( name="dia-Si-111-pandey", struct=dia_111_pandey('Si', 5.432) ),
        dict( name="dia-Si-110", struct=dia_110('Si', 5.432), r_Jm2=1.23 ),
        dict( name="dia-Si-100", struct=dia_100('Si', 5.432), r_Jm2=1.90 ),
        dict( name="dia-Si-100-2x1", struct=dia_100_2x1('Si', 5.432),
              r_Jm2=1.13 ),
        dict( name="dia-SiC-111", struct=dia_111('SiC', 4.321), r_Jm2=1.67 ),
        dict( name="dia-SiC-110", struct=dia_110('SiC', 4.321), r_Jm2=2.29 ),
        dict( name="dia-SiC-100", struct=dia_100('SiC', 4.321), r_Jm2=3.87 ),
        dict( name="dia-SiC-100-2x1", struct=dia_100_2x1('SiC', 4.321),
              r_Jm2=2.91 )
        ] ),
    ( Kumagai, Kumagai_CompMaterSci_39_457_Si,
      [ dict( name="dia-Si-111", struct=dia_111('Si', 5.429) ),
        dict( name="dia-Si-111-pandey", struct=dia_111_pandey('Si', 5.429) ),
        dict( name="dia-Si-110", struct=dia_110('Si', 5.429) ),
        dict( name="dia-Si-100", struct=dia_100('Si', 5.429) ),
        dict( name="dia-Si-100-2x1", struct=dia_100_2x1('Si', 5.429) ),
        ] ),
    ( Rebo2, {},
      [ dict( name="dia-C-111", struct=dia_111('C', 3.566) ),
        dict( name="dia-C-111-glide", struct=dia_111_glide('C', 3.566) ),
        dict( name="dia-C-111-pandey", struct=dia_111_pandey('C', 3.566) ),
        dict( name="dia-C-110", struct=dia_110('C', 3.566) ),
        dict( name="dia-C-100", struct=dia_100('C', 3.566) ),
        dict( name="dia-C-100-2x1", struct=dia_100_2x1('C', 3.566) ),
        ] ),
    ( Rebo2Scr, {},
      [ dict( name="dia-C-111", struct=dia_111('C', 3.566) ),
        dict( name="dia-C-111-glide", struct=dia_111_glide('C', 3.566) ),
        dict( name="dia-C-111-pandey", struct=dia_111_pandey('C', 3.566) ),
        dict( name="dia-C-110", struct=dia_110('C', 3.566) ),
        dict( name="dia-C-100", struct=dia_100('C', 3.566) ),
        dict( name="dia-C-100-2x1", struct=dia_100_2x1('C', 3.566) ),
        ] ),
    ( Tersoff, Tersoff_PRB_39_5566_Si_C,
      [ dict( name="dia-C-111", struct=dia_111('C', 3.566) ),
        dict( name="dia-C-111-glide", struct=dia_111_glide('C', 3.566) ),
        dict( name="dia-C-111-pandey", struct=dia_111_pandey('C', 3.566) ),
        dict( name="dia-C-110", struct=dia_110('C', 3.566) ),
        dict( name="dia-C-100", struct=dia_100('C', 3.566) ),
        dict( name="dia-C-100-2x1", struct=dia_100_2x1('C', 3.566) ),
        dict( name="dia-Si-111", struct=dia_111('Si', 5.432) ),
        dict( name="dia-Si-111-pandey", struct=dia_111_pandey('Si', 5.432) ),
        dict( name="dia-Si-110", struct=dia_110('Si', 5.432) ),
        dict( name="dia-Si-100", struct=dia_100('Si', 5.432) ),
        dict( name="dia-Si-100-2x1", struct=dia_100_2x1('Si', 5.432) ),
        dict( name="dia-SiC-111", struct=dia_111('SiC', 4.321) ),
        dict( name="dia-SiC-110", struct=dia_110('SiC', 4.321) ),
        dict( name="dia-SiC-100", struct=dia_100('SiC', 4.321) ),
        dict( name="dia-SiC-100-2x1", struct=dia_100_2x1('SiC', 4.321) )
        ] ),
    ]

###

class TestSurfaceEnergies(unittest.TestCase):

    def test_surface_energies(self):
        for pot, par, mats in tests:
            test_surface_energies(mats, pot, par, nx, vacuum, test=self)

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

        test_surface_energies(mats, pot, par, nx, vacuum, dump=True)
