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


'''
Basic test for Tersoff potential with ASE interface
'''

import unittest

from ase.build import bulk

from atomistica import Tersoff


class TestTersoff(unittest.TestCase):
    """Test Tersoff potential basic functionality"""

    def test_tersoff_silicon(self):
        """Test Tersoff potential energy calculation for silicon"""
        si = bulk('Si')
        t = Tersoff()
        si.calc = t
        energy = si.get_potential_energy()

        # Check that energy is a finite number
        self.assertTrue(abs(energy) < 1e10)

        # Check that we can get forces (tests array writability)
        forces = si.get_forces()
        self.assertEqual(forces.shape, (len(si), 3))


if __name__ == '__main__':
    unittest.main()
