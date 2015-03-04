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
Run Voronoi analysis of a trajectory, and store the results into that
trajectory.
"""

import os
import sys

from ase.io import NetCDFTrajectory

from atomistica.analysis import voropp

###

traj = NetCDFTrajectory(sys.argv[1], 'a')

for i, a in enumerate(traj):
    sys.stdout.write('=== {0}/{1} ===\r'.format(i+1, len(traj)))

    vol, area = voropp(a, q='%v %F', fast=True)

    a.set_array('voronoi_volume', vol)
    a.set_array('voronoi_surface_area', area)

    traj.write_arrays(a, i, ['voronoi_volume', 'voronoi_surface_area'])
traj.close()
