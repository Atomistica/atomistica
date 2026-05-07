# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
# ======================================================================

"""
Periodic boundary condition tests for the C++ backend.

Tests that energies and forces are consistent when atoms are moved across
cell boundaries (wrapping) and that forces sum to zero in periodic systems.
"""

import numpy as np
import pytest

ase = pytest.importorskip('ase')

import ase
from ase.lattice.cubic import Diamond, FaceCenteredCubic

import atomistica_cpp as a
from conftest import make_calc, assert_forces, assert_stress

DX  = 1e-6
DE  = 1e-6
TOL = 1e-2


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tersoff_calc():
    return make_calc(a.Tersoff, a.Tersoff_PRB_39_5566_Si_C)


def _brenner_calc():
    return make_calc(a.Brenner, a.Erhart_PRB_71_035211_SiC)


# ---------------------------------------------------------------------------
# Force sum (Newton's 3rd law)
# ---------------------------------------------------------------------------

class TestForceSum:
    """Forces must sum to zero in a periodic system (no net force on COM)."""

    def test_tersoff_Si(self):
        atoms = Diamond('Si', size=[2, 2, 2])
        atoms.calc = _tersoff_calc()
        forces = atoms.get_forces()
        assert np.max(np.abs(forces.sum(axis=0))) < 1e-10, \
            'Tersoff Si forces do not sum to zero'

    def test_brenner_C(self):
        atoms = Diamond('C', size=[2, 2, 2])
        atoms.calc = _brenner_calc()
        forces = atoms.get_forces()
        assert np.max(np.abs(forces.sum(axis=0))) < 1e-10, \
            'Brenner C forces do not sum to zero'

    def test_tersoff_rattled(self):
        atoms = Diamond('Si', size=[2, 2, 2])
        atoms.rattle(0.1)
        atoms.calc = _tersoff_calc()
        forces = atoms.get_forces()
        assert np.max(np.abs(forces.sum(axis=0))) < 1e-10, \
            'Tersoff Si (rattled) forces do not sum to zero'


# ---------------------------------------------------------------------------
# Energy invariance under cell-boundary wrapping
# ---------------------------------------------------------------------------

class TestWrapping:
    """Translating all atoms by a full lattice vector must not change energy."""

    def test_translate_full_cell(self):
        atoms = Diamond('Si', size=[2, 2, 2])
        atoms.calc = _tersoff_calc()
        E0 = atoms.get_potential_energy()

        cell = atoms.get_cell()
        atoms.translate(cell[0])
        E1 = atoms.get_potential_energy()
        assert abs(E1 - E0) < 1e-10, \
            f'Energy changed by {E1-E0:.2e} after full-cell translation'

    def test_translate_half_cell(self):
        """Translating by half a cell and back should give the same energy."""
        atoms = Diamond('C', size=[2, 2, 2])
        atoms.calc = _brenner_calc()
        E0 = atoms.get_potential_energy()

        cell = atoms.get_cell()
        atoms.translate(0.5 * cell[0])
        E1 = atoms.get_potential_energy()
        atoms.translate(-0.5 * cell[0])
        E2 = atoms.get_potential_energy()

        assert abs(E2 - E0) < 1e-10, \
            f'Energy changed after translate-untranslate: {E2-E0:.2e}'


# ---------------------------------------------------------------------------
# PBC vs. non-PBC
# ---------------------------------------------------------------------------

class TestPBCvsNoPBC:
    """A pair in PBC and its clone in a big non-PBC box should have same energy."""

    def test_dimer_pbc_vs_nopbc(self):
        r = 2.35
        # Periodic dimer (only nearest-neighbour interaction)
        a_pbc = ase.Atoms('SiSi',
                          positions=[[0, 0, 0], [r, 0, 0]],
                          cell=[20, 20, 20], pbc=True)
        a_pbc.calc = _tersoff_calc()
        E_pbc = a_pbc.get_potential_energy()

        # Non-periodic dimer in same big box
        a_nopbc = a_pbc.copy()
        a_nopbc.pbc = False
        a_nopbc.calc = _tersoff_calc()
        E_nopbc = a_nopbc.get_potential_energy()

        assert abs(E_pbc - E_nopbc) < 1e-8, \
            f'PBC ({E_pbc:.6f}) and no-PBC ({E_nopbc:.6f}) differ for dimer'


# ---------------------------------------------------------------------------
# Force/stress consistency under PBC
# ---------------------------------------------------------------------------

class TestForcesStressPBC:

    def test_tersoff_Si_forces(self):
        atoms = Diamond('Si', size=[2, 2, 2])
        atoms.translate([0.1, 0.1, 0.1])   # avoid exact-boundary atoms
        atoms.calc = _tersoff_calc()
        assert_forces(atoms, dx=DX, tol=TOL, msg='Tersoff Si PBC forces ')

    def test_tersoff_Si_stress(self):
        atoms = Diamond('Si', size=[2, 2, 2])
        atoms.translate([0.1, 0.1, 0.1])
        atoms.calc = _tersoff_calc()
        assert_stress(atoms, de=DE, tol=TOL, msg='Tersoff Si PBC stress ')

    def test_brenner_C_forces(self):
        atoms = Diamond('C', size=[2, 2, 2])
        atoms.translate([0.1, 0.1, 0.1])
        atoms.calc = _brenner_calc()
        assert_forces(atoms, dx=DX, tol=TOL, msg='Brenner C PBC forces ')

    def test_no_pbc_forces(self):
        """Forces are consistent with no PBC (molecule-like)."""
        atoms = ase.Atoms('SiSi',
                          positions=[[0, 0, 0], [2.35, 0, 0]],
                          cell=[20, 20, 20], pbc=False)
        atoms.calc = _tersoff_calc()
        assert_forces(atoms, dx=DX, tol=TOL, msg='Tersoff no-PBC forces ')
