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
Test if distance-energy and distance-force functions of dimers are continuous
for Rebo2, Rebo2Scr, Kumagai, KumagaiScr
'''

import unittest

import numpy as np

from ase import Atoms
from atomistica import Rebo2, Rebo2Scr, Kumagai, KumagaiScr


class DimerTest(unittest.TestCase):

    def test_C2(self):
        vac = 4
        dist_min = 1.2
        dist_max = 3.1

        a = Atoms('CC', positions=[[0, 0, 0], [dist_min, 0, 0]])
        a.center(vacuum=vac)
        distances = np.linspace(dist_min, dist_max, 1000)

        for potential in [Rebo2(), Rebo2Scr()]:
            a.set_calculator(potential)
            energies = []
            forces = []
            for dist in distances:
                a[1].position[0] = dist + vac
                forces += [a.get_forces()[0][0]]
                energies += [a.get_potential_energy()]

            forces = np.array(forces)
            energies = np.array(energies)
            en_differences = np.abs(energies[1:] - energies[:-1])
            self.assertTrue(np.max(en_differences) < 0.05)
            self.assertTrue(np.max(np.abs(forces[1:] - forces[:-1])) < 0.5)

    def test_H2(self):
        vac = 4
        dist_min = 0.6
        dist_max = 1.8

        a = Atoms('HH', positions=[[0, 0, 0], [dist_min, 0, 0]])
        a.center(vacuum=vac)
        distances = np.linspace(dist_min, dist_max, 1000)

        for potential in [Rebo2(), Rebo2Scr()]:
            a.set_calculator(potential)
            energies = []
            forces = []
            for dist in distances:
                a[1].position[0] = dist + vac
                forces += [a.get_forces()[0][0]]
                energies += [a.get_potential_energy()]

            forces = np.array(forces)
            energies = np.array(energies)

            en_differences = np.abs(energies[1:] - energies[:-1])
            self.assertTrue(np.max(en_differences) < 0.02)
            self.assertTrue(np.max(np.abs(forces[1:] - forces[:-1])) < 0.2)

    def test_CH(self):
        vac = 4
        dist_min = 0.8
        dist_max = 1.9

        a = Atoms('CH', positions=[[0, 0, 0], [dist_min, 0, 0]])
        a.center(vacuum=vac)
        distances = np.linspace(dist_min, dist_max, 1000)

        for potential in [Rebo2(), Rebo2Scr()]:
            a.set_calculator(potential)
            energies = []
            forces = []
            for dist in distances:
                a[1].position[0] = dist + vac
                forces += [a.get_forces()[0][0]]
                energies += [a.get_potential_energy()]

            forces = np.array(forces)
            energies = np.array(energies)
            en_differences = np.abs(energies[1:] - energies[:-1])
            self.assertTrue(np.max(en_differences) < 0.03)
            self.assertTrue(np.max(np.abs(forces[1:] - forces[:-1])) < 0.3)

    def test_Si2(self):
        vac = 4
        dist_min = 1.8
        dist_max = 6.2

        a = Atoms('Si2', positions=[[0, 0, 0], [dist_min, 0, 0]])
        a.center(vacuum=vac)
        distances = np.linspace(dist_min, dist_max, 1000)

        for potential in [Kumagai(), KumagaiScr()]:
            a.set_calculator(potential)
            energies = []
            forces = []
            for dist in distances:
                a[1].position[0] = dist + vac
                forces += [a.get_forces()[0][0]]
                energies += [a.get_potential_energy()]

            forces = np.array(forces)
            energies = np.array(energies)
            en_differences = np.abs(energies[1:] - energies[:-1])
            self.assertTrue(np.max(en_differences) < 0.08)
            self.assertTrue(np.max(np.abs(forces[1:] - forces[:-1])) < 0.4)


if __name__ == '__main__':
    unittest.main()
