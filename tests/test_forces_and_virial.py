# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2024) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
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

"""
Force and virial consistency tests for atomistica potentials.

Each test computes analytical forces/stress and compares them to numerical
finite-difference values. The tolerance is 1% (tol=1e-2).

Notes on REBO2 / REBO2Scr:
  The C++ REBO2 implementation computes simplified forces that omit
  bond-order derivative contributions (angular terms). Forces are only
  tested for simple dimers where the simplified forces are exact.
"""

import math

import numpy as np
import pytest

ase = pytest.importorskip('ase')

import ase.io
from ase.lattice.cubic import Diamond, FaceCenteredCubic, BodyCenteredCubic
from ase.lattice.compounds import B3

import atomistica as a
from conftest import (assert_forces, assert_stress, fortran_test_file,
                      make_calc)

# ---------------------------------------------------------------------------
# Tolerance and displacement
# ---------------------------------------------------------------------------

DX   = 1e-6
DE   = 1e-6
TOL  = 1e-2
SX   = 2   # supercell size

# ---------------------------------------------------------------------------
# Helper: perturb then test
# ---------------------------------------------------------------------------

def _check(atoms, pot_name, struct_name, tol=TOL):
    """Translate by (0.1,0.1,0.1) then check forces and stress."""
    atoms.translate([0.1, 0.1, 0.1])
    msg = f'{pot_name} / {struct_name}: '
    assert_forces(atoms, dx=DX, tol=tol, msg=msg)
    assert_stress(atoms, de=DE, tol=tol, msg=msg)
    # Rattle and check again
    atoms.rattle(0.1)
    assert_forces(atoms, dx=DX, tol=tol, msg=msg + '(rattled) ')
    assert_stress(atoms, de=DE, tol=tol, msg=msg + '(rattled) ')


# ===========================================================================
# Tersoff
# ===========================================================================

class TestTersoff:
    def _make(self, param):
        return make_calc(a.Tersoff, param)

    def test_dia_C(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Tersoff_PRB_39_5566_Si_C)
        _check(atoms, 'Tersoff(Si-C)', 'dia-C')

    def test_dia_Si(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = self._make(a.Tersoff_PRB_39_5566_Si_C)
        _check(atoms, 'Tersoff(Si-C)', 'dia-Si')

    def test_dia_SiC(self):
        atoms = B3(['Si', 'C'], latticeconstant=4.3596, size=[SX, SX, SX])
        atoms.calc = self._make(a.Tersoff_PRB_39_5566_Si_C)
        _check(atoms, 'Tersoff(Si-C)', 'dia-Si-C')

    def test_aC(self):
        p = fortran_test_file('aC_small.cfg')
        atoms = ase.io.read(p)
        atoms.calc = self._make(a.Tersoff_PRB_39_5566_Si_C)
        _check(atoms, 'Tersoff(Si-C)', 'a-C')

    def test_matsunaga_C(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N)
        _check(atoms, 'Tersoff(Matsunaga)', 'dia-C')

    def test_goumri_said_Al(self):
        atoms = FaceCenteredCubic('Al', latticeconstant=4.05, size=[SX, SX, SX])
        atoms.calc = self._make(a.Goumri_Said_ChemPhys_302_135_Al_N)
        _check(atoms, 'Tersoff(Goumri-Said)', 'fcc-Al')


class TestTersoffScr:
    def _make(self, param):
        return make_calc(a.TersoffScr, param)

    def test_dia_C(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Tersoff_PRB_39_5566_Si_C__Scr)
        _check(atoms, 'TersoffScr(Si-C)', 'dia-C')

    def test_dia_Si(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = self._make(a.Tersoff_PRB_39_5566_Si_C__Scr)
        _check(atoms, 'TersoffScr(Si-C)', 'dia-Si')

    def test_matsunaga_Scr_C(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr)
        _check(atoms, 'TersoffScr(Matsunaga)', 'dia-C')


# ===========================================================================
# Brenner
# ===========================================================================

class TestBrenner:
    def _make(self, param):
        return make_calc(a.Brenner, param)

    def test_dia_C_erhart(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Erhart_PRB_71_035211_SiC)
        _check(atoms, 'Brenner(Erhart)', 'dia-C')

    def test_dia_Si_erhart(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = self._make(a.Erhart_PRB_71_035211_SiC)
        _check(atoms, 'Brenner(Erhart)', 'dia-Si')

    def test_aC_erhart(self):
        p = fortran_test_file('aC_small.cfg')
        atoms = ase.io.read(p)
        atoms.calc = self._make(a.Erhart_PRB_71_035211_SiC)
        _check(atoms, 'Brenner(Erhart)', 'a-C')

    def test_dia_C_brenner_I(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Brenner_PRB_42_9458_C_I)
        _check(atoms, 'Brenner(PRB42-I)', 'dia-C')

    def test_dia_C_brenner_II(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Brenner_PRB_42_9458_C_II)
        _check(atoms, 'Brenner(PRB42-II)', 'dia-C')

    def test_dia_C_albe(self):
        # Pt-C potential: test with pure C to avoid the need for Pt structures
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Albe_PRB_65_195124_PtC)
        _check(atoms, 'Brenner(Albe)', 'dia-C')

    def test_dia_C_henriksson(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Henriksson_PRB_79_144107_FeC)
        _check(atoms, 'Brenner(Henriksson)', 'dia-C')

    def test_kioseoglou_Al(self):
        atoms = FaceCenteredCubic('Al', latticeconstant=4.05, size=[SX, SX, SX])
        atoms.calc = self._make(a.Kioseoglou_PSSb_245_1118_AlN)
        _check(atoms, 'Brenner(Kioseoglou)', 'fcc-Al')


class TestBrennerScr:
    def _make(self, param):
        return make_calc(a.BrennerScr, param)

    def test_dia_C_erhart_scr(self):
        atoms = Diamond('C', size=[SX, SX, SX])
        atoms.calc = self._make(a.Erhart_PRB_71_035211_SiC__Scr)
        _check(atoms, 'BrennerScr(Erhart)', 'dia-C')

    def test_dia_Si_erhart_scr(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = self._make(a.Erhart_PRB_71_035211_SiC__Scr)
        _check(atoms, 'BrennerScr(Erhart)', 'dia-Si')


# ===========================================================================
# Kumagai
# ===========================================================================

class TestKumagai:
    def test_dia_Si(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = make_calc(a.Kumagai, a.Kumagai_CompMaterSci_39_457_Si)
        _check(atoms, 'Kumagai(Si)', 'dia-Si')

    def test_dia_Si_scr(self):
        atoms = Diamond('Si', size=[SX, SX, SX])
        atoms.calc = make_calc(a.KumagaiScr, a.Kumagai_CompMaterSci_39_457_Si__Scr)
        _check(atoms, 'KumagaiScr(Si)', 'dia-Si')


# ===========================================================================
# Juslin
# ===========================================================================

class TestJuslin:
    def test_bcc_W(self):
        atoms = BodyCenteredCubic('W', latticeconstant=3.165, size=[SX, SX, SX])
        atoms.calc = make_calc(a.Juslin, a.Juslin_JAP_98_123520_WCH)
        _check(atoms, 'Juslin(W-C-H)', 'bcc-W')


# ===========================================================================
# EAM
# ===========================================================================

class TestEAM:
    def test_fcc_Au_funcfl(self):
        fn = fortran_test_file('Au_u3.eam')
        from atomistica import TabulatedEAM
        pot = TabulatedEAM()
        pot.load(fn)
        from atomistica import Atomistica
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[SX, SX, SX])
        atoms.rattle(0.1)
        atoms.calc = calc
        _check(atoms, 'TabulatedEAM(Au_u3)', 'fcc-Au')

    def test_fcc_Au_alloy(self):
        fn = fortran_test_file('Au-Grochola-JCP05.eam.alloy')
        from atomistica import TabulatedAlloyEAM
        pot = TabulatedAlloyEAM()
        pot.load(fn)
        from atomistica import Atomistica
        calc = Atomistica(pot)
        atoms = FaceCenteredCubic('Au', latticeconstant=4.08, size=[SX, SX, SX])
        atoms.rattle(0.1)
        atoms.calc = calc
        _check(atoms, 'TabulatedAlloyEAM(Au-Grochola)', 'fcc-Au')


# ===========================================================================
# Coulomb
# ===========================================================================

class TestCoulomb:
    def _nacl(self):
        """Return a NaCl-like structure with ±1 charges."""
        from ase.lattice.compounds import NaCl
        atoms = NaCl(['Na', 'Cl'], latticeconstant=5.64, size=[SX, SX, SX])
        syms = np.array(atoms.get_chemical_symbols())
        charges = np.where(syms == 'Na', 1.0, -1.0)
        atoms.set_array('charges', charges)
        return atoms

    def test_direct_coulomb_forces(self):
        atoms = self._nacl()
        atoms.pbc = False
        atoms.center(vacuum=5.0)
        from atomistica import DirectCoulomb, Atomistica
        atoms.calc = Atomistica(DirectCoulomb())
        assert_forces(atoms, dx=DX, tol=TOL, msg='DirectCoulomb forces ')

    def test_wolf_coulomb_forces(self):
        atoms = self._nacl()
        from atomistica import WolfCoulomb, Atomistica
        atoms.calc = Atomistica(WolfCoulomb(cutoff=8.0, alpha=0.3))
        assert_forces(atoms, dx=DX, tol=TOL, msg='WolfCoulomb forces ')
        assert_stress(atoms, de=DE, tol=TOL, msg='WolfCoulomb stress ')


# ===========================================================================
# REBO2 — only dimer-level tests where simplified forces are exact
# ===========================================================================

class TestREBO2:
    def _make_c_dimer(self, r=1.5):
        atoms = ase.Atoms('CC',
                          positions=[[0, 0, 0], [r, 0, 0]],
                          cell=[10, 10, 10], pbc=False)
        atoms.center()
        return atoms

    def test_c_dimer_forces(self):
        """For a simple dimer, simplified forces == exact forces."""
        atoms = self._make_c_dimer()
        from atomistica import REBO2, Atomistica
        pot = REBO2(); pot.load_default_parameters()
        atoms.calc = Atomistica(pot)
        assert_forces(atoms, dx=1e-5, tol=1e-3, msg='REBO2 C-dimer forces ')

    def test_rebo2scr_c_dimer_forces(self):
        atoms = self._make_c_dimer()
        from atomistica import REBO2Scr, Atomistica
        pot = REBO2Scr(); pot.load_default_parameters()
        atoms.calc = Atomistica(pot)
        assert_forces(atoms, dx=1e-5, tol=1e-3, msg='REBO2Scr C-dimer forces ')
