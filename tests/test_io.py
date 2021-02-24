'''
Test for IO in Rebo2, Rebo2Scr
'''

import unittest

import numpy as np

from ase import Atoms, io
from atomistica import Rebo2, Rebo2Scr


class IOTest(unittest.TestCase):
    def test_energy(self):
        vac = 8
        dist_min = 1.2

        atoms = Atoms('CC', positions=[[0, 0, 0], [dist_min, 0, 0]])
        atoms.center(vacuum=vac)

        for calc in [Rebo2(), Rebo2Scr()]:
            atoms.calc = calc
            energy = atoms.get_potential_energy()
            forces_ac = atoms.get_forces()
            stress = atoms.get_stress()

            fname = 'structure.traj'
            atoms.write(fname)
            atoms = io.read(fname)

            self.assertTrue(np.abs(energy - atoms.get_potential_energy())
                            < 1e-10)
            self.assertTrue((np.abs(forces_ac - atoms.get_forces())
                             < 1e-10).all())
            self.assertTrue((np.abs(stress - atoms.get_stress())
                             < 1e-10).all())


if __name__ == '__main__':
    unittest.main()
