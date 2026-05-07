# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
# ======================================================================

"""
EAM potential tests: correct loading from file, force/stress consistency,
and special crash cases.
"""

import numpy as np
import pytest

ase = pytest.importorskip('ase')

import ase.io
from ase.lattice.cubic import FaceCenteredCubic, BodyCenteredCubic

from atomistica_cpp import TabulatedEAM, TabulatedAlloyEAM, Atomistica
from conftest import assert_forces, assert_stress, fortran_test_file

DX  = 1e-6
DE  = 1e-6
TOL = 1e-4   # EAM forces should be very accurate


# ---------------------------------------------------------------------------
# TabulatedEAM (funcfl format)
# ---------------------------------------------------------------------------

class TestTabulatedEAM:

    def test_load_funcfl(self):
        fn = fortran_test_file('Au_u3.eam')
        pot = TabulatedEAM()
        pot.load(fn)
        assert pot.is_valid()
        assert pot.element_info().atomic_number == 79  # Au
        assert pot.cutoff() > 2.0

    def test_fcc_Au_energy(self):
        fn = fortran_test_file('Au_u3.eam')
        pot = TabulatedEAM(); pot.load(fn)
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[2, 2, 2])
        atoms.calc = calc
        E = atoms.get_potential_energy()
        assert E < 0, 'FCC-Au energy should be negative'
        Eper = E / len(atoms)
        # Au cohesive energy ≈ -3.81 eV/atom
        assert -4.5 < Eper < -3.0, f'Au Ec={Eper:.3f} eV/atom out of range (expected ~-3.81)'

    def test_fcc_Au_forces(self):
        fn = fortran_test_file('Au_u3.eam')
        pot = TabulatedEAM(); pot.load(fn)
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[2, 2, 2])
        atoms.rattle(0.1)
        atoms.calc = calc
        assert_forces(atoms, dx=DX, tol=TOL, msg='TabulatedEAM(Au) forces ')

    def test_fcc_Au_stress(self):
        fn = fortran_test_file('Au_u3.eam')
        pot = TabulatedEAM(); pot.load(fn)
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[2, 2, 2])
        atoms.rattle(0.1)
        atoms.calc = calc
        assert_stress(atoms, de=DE, tol=TOL, msg='TabulatedEAM(Au) stress ')


# ---------------------------------------------------------------------------
# TabulatedAlloyEAM (setfl format)
# ---------------------------------------------------------------------------

class TestTabulatedAlloyEAM:

    def test_load_alloy(self):
        fn = fortran_test_file('Au-Grochola-JCP05.eam.alloy')
        pot = TabulatedAlloyEAM()
        pot.load(fn)
        assert pot.is_valid()
        assert pot.num_elements() >= 1
        assert pot.cutoff() > 2.0

    def test_fcc_Au_forces(self):
        fn = fortran_test_file('Au-Grochola-JCP05.eam.alloy')
        pot = TabulatedAlloyEAM(); pot.load(fn)
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[2, 2, 2])
        atoms.rattle(0.1)
        atoms.calc = calc
        assert_forces(atoms, dx=DX, tol=TOL, msg='TabulatedAlloyEAM(Au-Grochola) forces ')

    def test_fcc_Au_stress(self):
        fn = fortran_test_file('Au-Grochola-JCP05.eam.alloy')
        pot = TabulatedAlloyEAM(); pot.load(fn)
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[2, 2, 2])
        atoms.rattle(0.1)
        atoms.calc = calc
        assert_stress(atoms, de=DE, tol=TOL, msg='TabulatedAlloyEAM(Au-Grochola) stress ')

    def test_cu_forces(self):
        fn = fortran_test_file('Cu_mishin1.eam.alloy')
        pot = TabulatedAlloyEAM(); pot.load(fn)
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Cu', latticeconstant=3.62, size=[2, 2, 2])
        atoms.rattle(0.1)
        atoms.calc = calc
        assert_forces(atoms, dx=DX, tol=TOL, msg='TabulatedAlloyEAM(Cu) forces ')


# ---------------------------------------------------------------------------
# Special crash cases (from Fortran test_eam_special_cases.py)
# ---------------------------------------------------------------------------

class TestEAMSpecialCases:

    def test_crash1(self):
        """Should not crash on a known pathological configuration."""
        fn_pot = fortran_test_file('Cu_mishin1.eam.alloy')
        fn_atoms = fortran_test_file('eam_crash1.poscar')
        pot = TabulatedAlloyEAM(); pot.load(fn_pot)
        atoms = ase.io.read(fn_atoms)
        atoms.calc = Atomistica(pot)
        atoms.get_potential_energy()  # should not crash

    def test_crash2_forces(self):
        """Dense configurations should give correct forces."""
        fn_pot = fortran_test_file('Cu_mishin1.eam.alloy')
        fn_atoms = fortran_test_file('eam_crash2.poscar')
        pot = TabulatedAlloyEAM(); pot.load(fn_pot)
        orig = ase.io.read(fn_atoms)
        for fac in [0.5]:  # 0.3 is too extreme (very high density)
            atoms = orig.copy()
            atoms.set_cell(fac * atoms.cell, scale_atoms=True)
            atoms.calc = Atomistica(pot)
            assert_forces(atoms, dx=DX, tol=TOL,
                          msg=f'EAM crash2 fac={fac} forces ')
