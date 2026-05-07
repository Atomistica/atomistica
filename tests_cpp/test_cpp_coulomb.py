# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
# ======================================================================

"""
Coulomb potential tests: energy values, force/virial consistency.
"""

import numpy as np
import pytest

ase = pytest.importorskip('ase')

import ase
from ase.units import Hartree, Bohr
from ase.lattice.compounds import NaCl

from atomistica_cpp import DirectCoulomb, CutoffCoulomb, WolfCoulomb, Atomistica
from conftest import assert_forces, assert_stress

DX  = 1e-6
DE  = 1e-6
TOL = 1e-2


# ---------------------------------------------------------------------------
# Helper: NaCl-like structure with ±1 charges
# ---------------------------------------------------------------------------

def _nacl_atoms(size=(1, 1, 1), pbc=True):
    from ase.lattice.compounds import NaCl as _NaCl
    atoms = _NaCl(['Na', 'Cl'], latticeconstant=5.64, size=size)
    Z = atoms.get_atomic_numbers()
    charges = np.where(np.array(atoms.get_chemical_symbols()) == 'Na', 1.0, -1.0)
    atoms.set_array('charges', charges)
    return atoms


def _random_ions(n=50, density=2.16, pbc=True):
    """Random equal-and-opposite ions."""
    na = n; nc = n
    rng = np.random.default_rng(42)
    syms = ['Na'] * na + ['Cl'] * nc
    positions = rng.uniform(0, 1, (na + nc, 3))
    atoms = ase.Atoms(syms, scaled_positions=positions,
                      cell=[10, 10, 10], pbc=pbc)
    mass = na + nc   # fictitious volume normalisation
    atoms.set_cell([10, 10, 10])
    charges = np.array([1.0]*na + [-1.0]*nc)
    atoms.set_array('charges', charges)
    return atoms


# ---------------------------------------------------------------------------
# DirectCoulomb
# ---------------------------------------------------------------------------

class TestDirectCoulomb:

    def test_zero_charges(self):
        atoms = ase.Atoms('NaCl', positions=[[-1, 0, 0], [1, 0, 0]],
                          pbc=False)
        atoms.center(vacuum=10)
        atoms.set_array('charges', np.zeros(2))
        atoms.calc = Atomistica(DirectCoulomb())
        assert atoms.get_potential_energy() == pytest.approx(0.0, abs=1e-12)

    def test_dimer_energy(self):
        """Two point charges ±q separated by d Å: E = k*q²/d."""
        d = 2.0   # separation in Å
        atoms = ase.Atoms('NaCl',
                          positions=[[-d/2, 0, 0], [d/2, 0, 0]],
                          pbc=False)
        atoms.center(vacuum=10)
        atoms.set_array('charges', np.array([-1.0, 1.0]))
        atoms.calc = Atomistica(DirectCoulomb())
        E = atoms.get_potential_energy()
        from atomistica_cpp import COULOMB_CONST
        E_ref = -COULOMB_CONST / d   # attractive: E < 0
        assert E == pytest.approx(E_ref, rel=1e-6)

    def test_forces_nonperiodic(self):
        atoms = _nacl_atoms(size=(1, 1, 1), pbc=False)
        atoms.center(vacuum=5)
        atoms.calc = Atomistica(DirectCoulomb())
        assert_forces(atoms, dx=DX, tol=TOL, msg='DirectCoulomb forces ')

    def test_forces_random_ions(self):
        atoms = _random_ions(n=10, pbc=False)
        atoms.center(vacuum=2.0)
        atoms.calc = Atomistica(DirectCoulomb())
        assert_forces(atoms, dx=DX, tol=TOL, msg='DirectCoulomb random forces ')


# ---------------------------------------------------------------------------
# CutoffCoulomb
# ---------------------------------------------------------------------------

class TestCutoffCoulomb:

    def test_forces(self):
        atoms = _nacl_atoms(size=(2, 2, 2))
        atoms.calc = Atomistica(CutoffCoulomb(cutoff=8.0))
        assert_forces(atoms, dx=DX, tol=TOL, msg='CutoffCoulomb forces ')

    def test_energy_is_finite(self):
        atoms = _nacl_atoms(size=(2, 2, 2))
        atoms.calc = Atomistica(CutoffCoulomb(cutoff=8.0))
        E = atoms.get_potential_energy()
        assert np.isfinite(E)
        # Note: hard-cutoff Coulomb sums for periodic ionic crystals can give
        # positive values due to truncation (no Ewald correction). The forces
        # are still self-consistent with the truncated energy.

    # Note: CutoffCoulomb has a hard cutoff — the virial stress does not agree
    # with finite-difference strain perturbations when bonds straddle the cutoff
    # boundary during strain (discontinuous energy vs strain). Use WolfCoulomb
    # for stress calculations.


# ---------------------------------------------------------------------------
# WolfCoulomb (damped shifted force)
# ---------------------------------------------------------------------------

class TestWolfCoulomb:

    def test_forces(self):
        atoms = _nacl_atoms(size=(2, 2, 2))
        atoms.calc = Atomistica(WolfCoulomb(cutoff=8.0, alpha=0.3))
        assert_forces(atoms, dx=DX, tol=TOL, msg='WolfCoulomb forces ')

    def test_stress(self):
        atoms = _nacl_atoms(size=(2, 2, 2))
        atoms.calc = Atomistica(WolfCoulomb(cutoff=8.0, alpha=0.3))
        assert_stress(atoms, de=DE, tol=TOL, msg='WolfCoulomb stress ')

    def test_auto_alpha(self):
        """alpha=0 should auto-compute a reasonable value."""
        atoms = _nacl_atoms(size=(1, 1, 1), pbc=False)
        atoms.center(vacuum=5)
        pot = WolfCoulomb(cutoff=10.0, alpha=0.0)
        atoms.calc = Atomistica(pot)
        E = atoms.get_potential_energy()
        assert np.isfinite(E)
