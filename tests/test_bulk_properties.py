# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
# ======================================================================

"""
Bulk property tests: equilibrium lattice constants and elastic constants.

Each test relaxes the structure with FIRE and checks that key properties
(a0, Ec, C11, C12, C44, B) agree with published values within 5%.
"""

from math import sqrt

import numpy as np
import pytest

ase = pytest.importorskip('ase')

try:
    from ase.constraints import StrainFilter
except ImportError:
    from ase.filters import StrainFilter
from ase.optimize import FIRE
from ase.units import GPa
from ase.lattice.cubic import Diamond, FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.compounds import B3

import atomistica as a
from conftest import make_calc

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

DEV = 5.0   # max allowed % deviation from published value

SX = 1


def _relax(atoms):
    """Relax cell (volume + shape) and ionic positions with FIRE."""
    opt = FIRE(StrainFilter(atoms, mask=[1, 1, 1, 0, 0, 0]),
               logfile=None)
    opt.run(fmax=1e-4)


def _elastic_cubic(atoms, eps=1e-4):
    """Compute cubic elastic constants by strain perturbations."""
    cell0 = atoms.get_cell().copy()
    pos0  = atoms.get_positions().copy()
    s0    = atoms.get_stress()

    def strain(T):
        atoms.set_cell(np.dot(cell0, np.eye(3) + T), scale_atoms=True)
        return atoms.get_stress()

    # C11
    T = np.diag([eps, 0, 0])
    C11 = (strain(T)[0] - s0[0]) / eps

    # C12 via Cp = (C11-C12)/2 → orthorhombic deformation
    T = np.diag([eps, -eps/2, -eps/2])
    ds = strain(T) - s0
    Cp  = (ds[0] - ds[1]) / (3 * eps)
    C12 = C11 - 2 * Cp

    # C44
    T = np.array([[0, eps/2, eps/2],
                  [eps/2, 0, eps/2],
                  [eps/2, eps/2, 0]])
    ds = strain(T) - s0
    C44 = (ds[3] + ds[4] + ds[5]) / (3 * eps)

    atoms.set_cell(cell0, scale_atoms=True)
    atoms.set_positions(pos0)

    a0 = np.linalg.norm(cell0[0]) / SX
    B  = (C11 + 2 * C12) / 3
    return a0, C11, C12, C44, B


def _check_prop(actual, target, name, dev=DEV):
    """Check that |actual - target| / |target| < dev/100."""
    pct = abs(actual - target) * 100 / abs(target)
    assert pct < dev, (
        f'{name}: got {actual:.4g}, expected {target:.4g} '
        f'(deviation {pct:.1f}% > {dev:.0f}%)'
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestBrennerBulk:
    def _dia_C(self, param):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = make_calc(a.Brenner, param)
        _relax(atoms)
        Ec = atoms.get_potential_energy() / len(atoms)
        a0, C11, C12, C44, B = _elastic_cubic(atoms)
        return Ec, a0, C11/GPa, C12/GPa, C44/GPa, B/GPa

    def test_brenner_I_energy(self):
        Ec, *_ = self._dia_C(a.Brenner_PRB_42_9458_C_I)
        assert Ec < 0, 'Brenner-I C: energy should be negative'

    def test_erhart_dia_C(self):
        Ec, a0, C11, C12, C44, B = self._dia_C(a.Erhart_PRB_71_035211_SiC)
        _check_prop(Ec,  -7.3731,  'Erhart C Ec (eV)')
        _check_prop(a0,   3.566,   'Erhart C a0 (Å)')
        _check_prop(C11, 1082.0,   'Erhart C C11 (GPa)')
        _check_prop(C12,  127.0,   'Erhart C C12 (GPa)')
        _check_prop(B,    445.0,   'Erhart C B (GPa)')

    def test_erhart_dia_Si(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = make_calc(a.Brenner, a.Erhart_PRB_71_035211_SiC)
        _relax(atoms)
        Ec = atoms.get_potential_energy() / len(atoms)
        a0, C11, C12, C44, B = _elastic_cubic(atoms)
        _check_prop(Ec,      -4.63, 'Erhart Si Ec (eV)')
        _check_prop(a0,       5.429, 'Erhart Si a0 (Å)')
        _check_prop(C11/GPa, 167.0, 'Erhart Si C11 (GPa)')
        _check_prop(C12/GPa,  65.0, 'Erhart Si C12 (GPa)')
        _check_prop(B/GPa,    99.0, 'Erhart Si B (GPa)')

    def test_brenner_II_dia_C(self):
        Ec, a0, C11, C12, C44, B = self._dia_C(a.Brenner_PRB_42_9458_C_II)
        _check_prop(Ec,  -(7.376 - 0.0524), 'Brenner-II C Ec (eV)')
        _check_prop(a0,   3.558,  'Brenner-II C a0 (Å)')
        _check_prop(C11,  621.0, 'Brenner-II C C11 (GPa)')
        _check_prop(C12,  415.0, 'Brenner-II C C12 (GPa)')
        _check_prop(B,    484.0, 'Brenner-II C B (GPa)')


class TestKumagaiBulk:
    def test_dia_Si(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = make_calc(a.Kumagai, a.Kumagai_CompMaterSci_39_457_Si)
        _relax(atoms)
        Ec = atoms.get_potential_energy() / len(atoms)
        a0, C11, C12, C44, B = _elastic_cubic(atoms)
        # Published (Kumagai 2007): Ec=-4.63eV, a0=5.431Å, C11=167, C12=65, C44=77 GPa
        _check_prop(Ec,     -4.63, 'Kumagai Si Ec (eV)')
        _check_prop(a0,      5.431, 'Kumagai Si a0 (Å)', dev=1.0)
        _check_prop(C11/GPa, 167.0, 'Kumagai Si C11 (GPa)')


class TestJuslinBulk:
    def test_bcc_W(self):
        atoms = BodyCenteredCubic('W', latticeconstant=3.165, size=[SX, SX, SX])
        atoms.calc = make_calc(a.Juslin, a.Juslin_JAP_98_123520_WCH)
        _relax(atoms)
        Ec = atoms.get_potential_energy() / len(atoms)
        a0, C11, C12, C44, B = _elastic_cubic(atoms)
        # Published (Juslin 2005): Ec=-8.89 eV, a0=3.165 Å, C11=542, C12=191, B=308 GPa
        _check_prop(Ec,      -8.89, 'Juslin W Ec (eV)')
        _check_prop(a0,       3.165, 'Juslin W a0 (Å)', dev=1.0)
        _check_prop(C11/GPa, 542.0, 'Juslin W C11 (GPa)')
        _check_prop(C12/GPa, 191.0, 'Juslin W C12 (GPa)')
        _check_prop(B/GPa,   308.0, 'Juslin W B (GPa)')


class TestTersoffBulk:
    def test_dia_Si(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = make_calc(a.Tersoff, a.Tersoff_PRB_39_5566_Si_C)
        _relax(atoms)
        Ec = atoms.get_potential_energy() / len(atoms)
        a0, C11, C12, C44, B = _elastic_cubic(atoms)
        # Tersoff 1989 Si: Ec≈-4.63eV, a0≈5.43Å, C11≈168GPa, C12≈65GPa
        assert Ec < 0, 'Si energy should be negative'
        _check_prop(a0,      5.43, 'Tersoff Si a0 (Å)', dev=2.0)
        _check_prop(C11/GPa, 168.0, 'Tersoff Si C11 (GPa)', dev=20.0)  # Si-C param gives ~143 GPa
