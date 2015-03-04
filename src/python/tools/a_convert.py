#! /usr/bin/env python

# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
# See the AUTHORS file in the top-level Atomistica directory.
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
Convenience command line tool for conversion between file formats. Similar to
ASE's ag, but with support for AMBER-style NetCDF files.
"""

import os
import sys

from optparse import OptionParser

import atomistica.io as io

import ase
from ase.io.bader import attach_charges

###

parser  = OptionParser()
parser.add_option("--toA", action="store_true", dest="toA", default=False)
parser.add_option("--toBohr", action="store_true", dest="toBohr", default=False)
parser.add_option("--supercell", action="store", dest="supercell")
parser.add_option("--wrap-to-cell", action="store_true", dest="wrap_to_cell", default=False)
parser.add_option("--cell", action="store", dest="cell")
parser.add_option("--cellfn", action="store", dest="cellfn")
parser.add_option("--center", action="store_true", dest="center", default=False)
parser.add_option("--acf", action="store", dest="acf")
parser.add_option("--clear-velocities", action="store_true", dest="clear_velocities", default=False)
parser.add_option("--specorder", action="store", dest="specorder")
parser.add_option("--format", action="store", dest="format")

(opt, args)  = parser.parse_args()

###

convlen = None

if opt.toA:
    convlen  = ase.units.Bohr
elif opt.toBohr:
    convlen  = 1.0/ase.units.Bohr

infn    = args[0]
outfn   = args[1]

a  = io.read(infn)

if 'shear_dx' in a.info:
    cx, cy, cz = a.cell
    assert abs(cx[1]) < 1e-12
    assert abs(cx[2]) < 1e-12
    assert abs(cy[0]) < 1e-12
    assert abs(cy[2]) < 1e-12
    assert abs(cz[0]) < 1e-12
    assert abs(cz[1]) < 1e-12
    dx, dy, dz = a.info['shear_dx']
    sx, sy, sz = a.cell.diagonal()
    a.set_cell([[sx,0,0],[0,sy,0],[dx,dy,sz]], scale_atoms=False)
    
if opt.cell is not None:
    cell = map(float, opt.cell.split(','))
    a.set_cell(cell)

if opt.cellfn is not None:
    io.read_cyc(a, opt.cellfn)

if opt.center is not None:
    if opt.center:
        a.center()

if opt.acf is not None:
    attach_charges(a, opt.acf)

if convlen is not None:
    a.set_cell(a.get_cell()*convlen, scale_atoms=True)

if opt.wrap_to_cell:
    a.set_scaled_positions(a.get_scaled_positions())

if opt.clear_velocities:
    a.set_momenta(None)

if opt.supercell is not None:
    supercell  = map(int, opt.supercell.split(','))
    a *= supercell

d = { }
if opt.specorder is not None:
    d['specorder'] = opt.specorder.split(',')

if opt.format is not None:
    d['format'] = opt.format

io.write(outfn, a, **d)
