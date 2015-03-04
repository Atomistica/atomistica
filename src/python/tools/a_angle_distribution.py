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
Compute angle distribution.
"""

import sys
from math import pi

import numpy as np

import ase
from atomistica.io import read

import atomistica.native as native

###

# Default parameters
nbins = 100
cutoff = 5.0
avgn = 100
outfn = 'angle_distribution.out'

import getopt
optlist, args = getopt.getopt(sys.argv[1:], '',
                              [ 'nbins=', 'cutoff=',
                                'avgn=', 'out=' ])

assert len(args) == 1
fn = args[0]
for key, value in optlist:
    if key == '--nbins':
        nbins = int(value)
    elif key == '--cutoff':
        cutoff = float(value)
    elif key == '--avgn':
        avgn = int(value)
    elif key == '--out':
        outfn = value

###

print '# fn = ', fn
print '# nbins = ', nbins
print '# cutoff = ', cutoff
print '# avgn = ', avgn
print '# outfn = ', outfn

###

a = read(fn)

print '{0} atoms.'.format(len(a))

p = native.from_atoms(a)
nl = native.neighbor_list(p, cutoff, avgn=avgn)

i, j, dr, abs_dr = nl.get_neighbors(p, vec=True)
pavg, pvar = native.angle_distribution(i, j, dr, nbins, cutoff)
r = np.linspace(0.0, 2*pi, nbins+1)
r = (r[1:]+r[:-1])/2

np.savetxt(outfn, np.transpose([r, pavg, pvar]))

