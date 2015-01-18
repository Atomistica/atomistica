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
Checkpointing class.
"""

import os
import atomistica.io as io

###

class Checkpoint(object):
    """
    Call function only if checkpoint file does not exist, otherwise read atoms
    object from checkpoint file. This allows to have multiple consecutive
    operations in a single script and restart from the latest one that
    completed.

    Example:

    a = Checkpoint(optimize_geometry, 'geometry_optimized.traj')(a)
    a = Checkpoint(compute_charge, 'charge_computed.traj')(a)
    """

    def __init__(self, func, fn):
        self.func = func
        self.fn = fn

    def __call__(self, *args):
        if os.path.exists(self.fn):
            print 'Reading configuration from checkpoint file %s...' % self.fn
            a = io.read(self.fn)
        else:
            a = self.func(*args)
            io.write(self.fn, a)
        return a

