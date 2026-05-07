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
Atomistica C++ — Modern C++ implementation of interatomic potentials.

Quick start::

    from atomistica_cpp import Tersoff, Tersoff_PRB_39_5566_Si_C
    from ase.lattice.cubic import Diamond

    atoms = Diamond('Si', latticeconstant=5.43)
    calc = Atomistica(Tersoff, Tersoff_PRB_39_5566_Si_C)
    atoms.calc = calc
    print(atoms.get_potential_energy())
"""

from ._atomistica_cpp import (
    # Core
    AtomicSystem,
    NeighborList,
    Neighbor,
    PotentialResults,

    # Math utilities
    CubicSpline,
    NonUniformSpline,
    CutoffResult,
    TrigOffCutoff,
    TrigOnCutoff,
    ExpCutoff,

    # Simple pair potentials
    LJCut,
    LJCutShift,
    BornMayer,
    Harmonic,
    DoubleHarmonic,
    R6,

    # Bond-order potentials — Tersoff
    TersoffElementParams,
    TersoffPairParams,
    ScreeningParams,
    Tersoff,
    TersoffScr,

    # Bond-order potentials — Brenner
    BrennerElementParams,
    BrennerPairParams,
    Brenner,
    BrennerScr,

    # Bond-order potentials — Kumagai
    KumagaiElementParams,
    KumagaiPairParams,
    Kumagai,
    KumagaiScr,

    # Bond-order potentials — Juslin
    Juslin,
    JuslinScr,

    # REBO2
    REBO2,
    REBO2Scr,
    REBO2_C_C,
    REBO2_C_H,
    REBO2_H_H,
    REBO2_C,
    REBO2_H,

    # EAM
    EAMElementInfo,
    TabulatedEAM,
    TabulatedAlloyEAM,

    # Coulomb
    COULOMB_CONST,
    DirectCoulomb,
    CutoffCoulomb,
    WolfCoulomb,
    PMECoulomb,
    FMMCoulomb,

    # Dispersion
    DFTD3Disp,

    # Tight-binding / DFTB
    TBElementParams,
    SCCParams,
    SolverParams,
    DenseHamiltonian,
    MaterialsDatabase,
    DFTB,

    # Parameter set discovery
    available_tersoff_parameters,
    available_brenner_parameters,
    available_kumagai_parameters,
    available_juslin_parameters,
)

try:
    from .ase_calculator import Atomistica
except ImportError:
    pass  # ASE not available; Atomistica calculator not imported

from .parameters import (
    # Tersoff
    Tersoff_PRB_39_5566_Si_C,
    Tersoff_PRB_39_5566_Si_C__Scr,
    Goumri_Said_ChemPhys_302_135_Al_N,
    Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N,
    Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr,
    # Brenner
    Erhart_PRB_71_035211_SiC,
    Erhart_PRB_71_035211_SiC__Scr,
    Albe_PRB_65_195124_PtC,
    Henriksson_PRB_79_144107_FeC,
    Kioseoglou_PSSb_245_1118_AlN,
    Brenner_PRB_42_9458_C_I,
    Brenner_PRB_42_9458_C_II,
    # Kumagai
    Kumagai_CompMaterSci_39_457_Si,
    Kumagai_CompMaterSci_39_457_Si__Scr,
    # Juslin
    Juslin_JAP_98_123520_WCH,
)

__all__ = [
    # Core
    'AtomicSystem', 'NeighborList', 'Neighbor', 'PotentialResults',
    # Math
    'CubicSpline', 'NonUniformSpline',
    'CutoffResult', 'TrigOffCutoff', 'TrigOnCutoff', 'ExpCutoff',
    # Simple pair potentials
    'LJCut', 'LJCutShift', 'BornMayer', 'Harmonic', 'DoubleHarmonic', 'R6',
    # BOPs
    'TersoffElementParams', 'TersoffPairParams', 'ScreeningParams',
    'Tersoff', 'TersoffScr',
    'BrennerElementParams', 'BrennerPairParams', 'Brenner', 'BrennerScr',
    'KumagaiElementParams', 'KumagaiPairParams', 'Kumagai', 'KumagaiScr',
    'Juslin', 'JuslinScr',
    'REBO2', 'REBO2Scr',
    'REBO2_C_C', 'REBO2_C_H', 'REBO2_H_H', 'REBO2_C', 'REBO2_H',
    # EAM
    'EAMElementInfo', 'TabulatedEAM', 'TabulatedAlloyEAM',
    # Coulomb
    'COULOMB_CONST', 'DirectCoulomb', 'CutoffCoulomb', 'WolfCoulomb',
    'PMECoulomb', 'FMMCoulomb',
    # Dispersion
    'DFTD3Disp',
    # TB/DFTB
    'TBElementParams', 'SCCParams', 'SolverParams',
    'DenseHamiltonian', 'MaterialsDatabase', 'DFTB',
    # Discovery
    'available_tersoff_parameters', 'available_brenner_parameters',
    'available_kumagai_parameters', 'available_juslin_parameters',
    # ASE calculator
    'Atomistica',
    # Parameter name constants
    'Tersoff_PRB_39_5566_Si_C', 'Tersoff_PRB_39_5566_Si_C__Scr',
    'Goumri_Said_ChemPhys_302_135_Al_N',
    'Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N',
    'Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr',
    'Erhart_PRB_71_035211_SiC', 'Erhart_PRB_71_035211_SiC__Scr',
    'Albe_PRB_65_195124_PtC', 'Henriksson_PRB_79_144107_FeC',
    'Kioseoglou_PSSb_245_1118_AlN',
    'Brenner_PRB_42_9458_C_I', 'Brenner_PRB_42_9458_C_II',
    'Kumagai_CompMaterSci_39_457_Si', 'Kumagai_CompMaterSci_39_457_Si__Scr',
    'Juslin_JAP_98_123520_WCH',
]
