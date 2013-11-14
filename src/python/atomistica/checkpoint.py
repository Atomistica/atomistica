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

"""
Checkpointing class.
"""

class Checkpoint(object):
    """
    Call function only if checkpoint file does not exist, otherwise read atoms
    object from checkpoint file. This allows to have multiple consecutive
    operations in a single script and restart from the latest one that
    completed.

    Example:

    a = Checkpoint(optimize_geometry, 'optimized_geometry.traj')(a)
    a = Checkpoint(compute_charge, 'charged_computed.traj')(a)
    """

    def __init__(self, func, fn):
        self.func = func
        self.fn = fn

    def __call__(self, *args):
        if os.path.exists(self.fn):
            a = io.read(self.fn)
            a.set_pbc(True)
        else:
            a = self.func(*args)
            io.write(self.fn, a)
        return a

