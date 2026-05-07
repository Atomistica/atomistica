# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
# ======================================================================

"""
Tests for the C++ NeighborList class: correctness under PBC, cell changes,
minimum image convention, and Verlet shell behaviour.
"""

import numpy as np
import pytest

ase = pytest.importorskip('ase')

import ase
import ase.io
from ase.lattice.cubic import Diamond

import atomistica as a
from conftest import fortran_test_file


# ---------------------------------------------------------------------------
# Helper: build AtomicSystem from an ASE Atoms object
# ---------------------------------------------------------------------------

def _ase_to_cpp(atoms):
    """Create and return (AtomicSystem, NeighborList)."""
    sys = a.AtomicSystem(len(atoms))
    sys.cell = np.array(atoms.cell).T
    sys.pbc  = list(atoms.pbc)
    sys.positions = atoms.positions.T
    sys.atomic_numbers = atoms.numbers
    return sys


# ---------------------------------------------------------------------------
# Basic correctness
# ---------------------------------------------------------------------------

class TestNeighborListCorrectness:

    def test_distances_consistent(self):
        """All distances reported by the NL must match explicit calculations."""
        p = fortran_test_file('aC_small.cfg')
        ase_a = ase.io.read(p)
        sys = _ase_to_cpp(ase_a)

        nl = a.NeighborList()
        cutoff = 5.0
        nl.set_cutoff(cutoff)
        nl.update(sys)

        pos = ase_a.get_positions()
        cell = np.array(ase_a.cell)

        for i in range(sys.num_atoms):
            for nb in nl.neighbors(i):
                j = nb.index
                cs = nb.cell_shift
                dr = pos[j] - pos[i] + cs[0]*cell[0] + cs[1]*cell[1] + cs[2]*cell[2]
                r  = np.linalg.norm(dr)
                assert r <= cutoff + 1e-10, (
                    f'Atom {i}→{j}: r={r:.4f} > cutoff={cutoff}')

    def test_symmetry(self):
        """Every bond (i→j) should also appear as (j→i)."""
        atoms = Diamond('C', size=[2, 2, 2])
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(2.0)
        nl.update(sys)

        bonds_fwd = set()
        for i in range(sys.num_atoms):
            for nb in nl.neighbors(i):
                bonds_fwd.add((i, nb.index))

        # Every (i,j) should have a matching (j,i)
        for (i, j) in bonds_fwd:
            assert (j, i) in bonds_fwd, f'Bond {i}→{j} has no reverse {j}→{i}'

    def test_num_pairs(self):
        """num_pairs should equal sum of all neighbor counts."""
        atoms = Diamond('Si', size=[2, 2, 2])
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(3.0)
        nl.update(sys)

        total = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))
        assert nl.num_pairs == total


# ---------------------------------------------------------------------------
# PBC handling
# ---------------------------------------------------------------------------

class TestPBC:

    def test_pbc_bond_across_boundary(self):
        """Two atoms on opposite sides of the cell should see each other."""
        atoms = ase.Atoms('CC',
                          positions=[[0.1, 0.5, 0.5],
                                     [0.9, 0.5, 0.5]],
                          cell=[1, 1, 1], pbc=True)
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(0.3)
        nl.update(sys)

        n_bonds = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))
        assert n_bonds == 2, f'Expected 2 bonds (PBC), got {n_bonds}'

    def test_no_pbc_no_bond_across_boundary(self):
        """Without PBC, the same pair should NOT see each other."""
        atoms = ase.Atoms('CC',
                          positions=[[0.1, 0.5, 0.5],
                                     [0.9, 0.5, 0.5]],
                          cell=[1, 1, 1], pbc=False)
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(0.3)
        nl.update(sys)

        n_bonds = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))
        assert n_bonds == 0, f'Expected 0 bonds (no PBC), got {n_bonds}'

    def test_partial_pbc(self):
        """Atoms should bond through the periodic directions only."""
        atoms = ase.Atoms('CC',
                          positions=[[0.1, 0.5, 0.5],
                                     [0.9, 0.5, 0.5]],
                          cell=[1, 1, 1], pbc=[True, False, False])
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(0.3)
        nl.update(sys)

        n_bonds = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))
        assert n_bonds == 2, f'Expected 2 bonds (x-PBC only), got {n_bonds}'


# ---------------------------------------------------------------------------
# Invalidation and update
# ---------------------------------------------------------------------------

class TestNeighborListUpdate:

    def test_invalidation_after_position_change(self):
        """After positions_changed(), nl must be updated before use."""
        atoms = Diamond('Si', size=[1, 1, 1])
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(3.0)
        nl.update(sys)
        n0 = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))

        # Move atoms far apart → fewer neighbours
        sys.positions = sys.positions * 3.0
        sys.positions_changed()
        nl.update(sys)
        n1 = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))

        assert n1 <= n0, 'Expected fewer neighbours after moving atoms apart'

    def test_verlet_shell_avoids_rebuild(self):
        """With a Verlet shell the neighbour list preserves bond connectivity
        for small displacements that keep all pairs inside the extended cutoff."""
        atoms = Diamond('Si', size=[2, 2, 2])
        atoms.translate([0.5, 0.5, 0.5])  # avoid boundary atoms
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        cutoff = 2.4  # strictly below nearest-neighbour distance ~2.35 Å
        nl.set_cutoff(cutoff)
        nl.set_verlet_shell(0.5)
        nl.update(sys)
        n0 = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))

        # Tiny rattle (0.001 Å) — within Verlet shell, no bonds cross cutoff
        np.random.seed(42)
        sys.positions = sys.positions + np.random.randn(*sys.positions.shape) * 0.001
        sys.positions_changed()
        nl.update(sys)
        n1 = sum(nl.num_neighbors(i) for i in range(sys.num_atoms))
        assert n1 == n0, (
            f'Neighbour count changed ({n0} → {n1}) with 0.001 Å rattle'
        )


# ---------------------------------------------------------------------------
# Self-image bond tests
# ---------------------------------------------------------------------------

class TestSelfImageBonds:
    """Atoms must see their own periodic images when the cell is smaller than the cutoff."""

    def test_single_atom_sees_self_images(self):
        """A single atom in a small PBC cell should have self as neighbor."""
        atoms = ase.Atoms('Si', positions=[[0,0,0]], cell=[5,5,5], pbc=True)
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(7.0)
        nl.update(sys)

        nbs = list(nl.neighbors(0))
        # Must find self-images
        self_images = [nb for nb in nbs if nb.index == 0]
        assert len(self_images) > 0, 'No self-image bonds found for single-atom cell'
        # 6 nearest self-images along ±x, ±y, ±z at distance 5 Å
        assert len(self_images) >= 6

    def test_bcc_2atom_cell_correct_coordination(self):
        """2-atom BCC cell: atom 0 must find 8 1st NN (atom 1) + 6 2nd NN (self)."""
        a_bcc = 3.165
        atoms = ase.Atoms(
            'WW',
            positions=[[0,0,0], [a_bcc/2]*3],
            cell=[[a_bcc,0,0],[0,a_bcc,0],[0,0,a_bcc]],
            pbc=True
        )
        sys = _ase_to_cpp(atoms)
        nl = a.NeighborList()
        nl.set_cutoff(4.0)
        nl.update(sys)

        cell = np.array(atoms.cell)
        count_1nn, count_2nn_self = 0, 0
        for nb in nl.neighbors(0):
            dr = atoms.positions[nb.index] - atoms.positions[0]
            dr += nb.cell_shift[0]*cell[0] + nb.cell_shift[1]*cell[1] + nb.cell_shift[2]*cell[2]
            r = np.linalg.norm(dr)
            if nb.index == 1 and abs(r - a_bcc*3**0.5/2) < 0.01:
                count_1nn += 1
            elif nb.index == 0 and abs(r - a_bcc) < 0.01:
                count_2nn_self += 1

        assert count_1nn == 8, f'Expected 8 1st NN, got {count_1nn}'
        assert count_2nn_self == 6, f'Expected 6 self-image 2nd NN, got {count_2nn_self}'

    def test_juslin_w_bcc_1x1x1(self):
        """BCC-W 1×1×1 cell (2 atoms) gives correct energy with self-image fix."""
        import atomistica as ac
        from ase.lattice.cubic import BodyCenteredCubic

        atoms = BodyCenteredCubic('W', latticeconstant=3.165, size=[1, 1, 1])
        atoms.calc = ac.Atomistica(ac.Juslin, ac.Juslin_JAP_98_123520_WCH)
        E_per_atom = atoms.get_potential_energy() / len(atoms)
        assert abs(E_per_atom - (-8.89)) < 0.05, (
            f'Juslin BCC-W 1x1x1 Ec={E_per_atom:.3f} eV/atom, expected -8.89'
        )
