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
Shared fixtures and utilities for the atomistica_cpp Python test suite.

Run from the tests_cpp/ directory or with:
    pytest atomistica_cpp/tests_cpp/
"""

from pathlib import Path
import numpy as np
import pytest

# ---------------------------------------------------------------------------
# Test data location: look for Fortran test fixtures next to this repo
# ---------------------------------------------------------------------------

_HERE = Path(__file__).parent
# Test data lives in tests/ alongside tests_cpp/
_TEST_DATA = _HERE.parent / 'tests'
# Backward compat: also accept data from the sibling Fortran repo if present
_FORTRAN_TESTS = _HERE.parent.parent / 'atomistica_fortran' / 'tests'

def fortran_test_file(name):
    """Return path to a test data file.

    Checks the local tests/ directory first, then the sibling Fortran repo.
    Raises pytest.skip if the file cannot be found in either location.
    """
    # Primary: local tests/ directory
    p = _TEST_DATA / name
    if p.exists():
        return str(p)
    # Fallback: sibling Fortran repo (multi-repo layout)
    p2 = _FORTRAN_TESTS / name
    if p2.exists():
        return str(p2)
    pytest.skip(f'Test data file not found: {name} (looked in {_TEST_DATA} and {_FORTRAN_TESTS})')


# ---------------------------------------------------------------------------
# Core numerical test utilities
# ---------------------------------------------------------------------------

def numerical_forces(atoms, dx=1e-6):
    """Compute forces by central finite differences.

    Returns (f_numerical, f_analytical, max_rms_error).
    """
    f0 = atoms.get_forces().copy()
    ffd = np.zeros_like(f0)

    pos0 = atoms.get_positions().copy()
    for iatom in range(len(atoms)):
        for c in range(3):
            pos = pos0.copy()
            pos[iatom, c] -= dx
            atoms.set_positions(pos)
            e1 = atoms.get_potential_energy()

            pos[iatom, c] += 2 * dx
            atoms.set_positions(pos)
            e2 = atoms.get_potential_energy()

            ffd[iatom, c] = -(e2 - e1) / (2 * dx)

    atoms.set_positions(pos0)

    df = ffd - f0
    max_err = np.sqrt(np.max(np.sum(df * df, axis=1)))
    return ffd, f0, max_err


def numerical_stress(atoms, de=1e-6):
    """Compute stress by central finite differences of strain.

    Returns (stress_numerical, stress_analytical, max_abs_error) in Voigt
    notation (xx, yy, zz, yz, xz, xy) in units of eV/Å³.
    """
    s0 = atoms.get_stress().copy()
    V0 = atoms.get_volume()
    cell0 = atoms.get_cell().copy()

    smat = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            eps = np.eye(3)
            eps[i, j] -= de
            atoms.set_cell(np.dot(cell0, eps), scale_atoms=True)
            e1 = atoms.get_potential_energy()

            eps[i, j] += 2 * de
            atoms.set_cell(np.dot(cell0, eps), scale_atoms=True)
            e2 = atoms.get_potential_energy()

            smat[i, j] = (e2 - e1) / (2 * de)

    atoms.set_cell(cell0, scale_atoms=True)

    sfd = np.array([
        smat[0, 0], smat[1, 1], smat[2, 2],
        (smat[1, 2] + smat[2, 1]) / 2,
        (smat[0, 2] + smat[2, 0]) / 2,
        (smat[0, 1] + smat[1, 0]) / 2,
    ]) / V0

    max_err = np.max(np.abs(sfd - s0))
    return sfd, s0, max_err


def assert_forces(atoms, dx=1e-6, tol=1e-2, msg=''):
    """Assert analytical forces agree with finite-difference forces."""
    ffd, f0, max_err = numerical_forces(atoms, dx=dx)
    assert max_err < tol, (
        f'{msg}Force mismatch: max_rms={max_err:.2e} (tol={tol:.2e})\n'
        f'  analytical: {f0}\n  numerical: {ffd}'
    )


def assert_stress(atoms, de=1e-6, tol=1e-2, msg=''):
    """Assert analytical stress agrees with finite-difference stress."""
    sfd, s0, max_err = numerical_stress(atoms, de=de)
    assert max_err < tol, (
        f'{msg}Stress mismatch: max_abs={max_err:.2e} (tol={tol:.2e})\n'
        f'  analytical: {s0}\n  numerical: {sfd}'
    )


# ---------------------------------------------------------------------------
# Calculator factory
# ---------------------------------------------------------------------------

def make_calc(PotClass, param=None):
    """Create an atomistica_cpp.Atomistica calculator.

    Parameters
    ----------
    PotClass : potential class
    param : str or None
        Parameter set name passed to ``pot.load_parameters()``.

    Returns
    -------
    Atomistica calculator
    """
    from atomistica_cpp import Atomistica
    return Atomistica(PotClass, param)
