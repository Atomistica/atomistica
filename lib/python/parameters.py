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
Parameter set name constants for Atomistica C++ potentials.

Each constant is the string passed to potential.load_parameters().
Screened variants (``__Scr`` suffix) are used with the screened potential
class (e.g. ``TersoffScr``).

Usage::

    from atomistica_cpp import TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr
    calc = TersoffScr()
    calc.load_parameters(Tersoff_PRB_39_5566_Si_C__Scr)
"""

# ---------------------------------------------------------------------------
# Tersoff potential parameter sets
# ---------------------------------------------------------------------------

#: Tersoff Si-C, Phys. Rev. B 39, 5566 (1989)
Tersoff_PRB_39_5566_Si_C = "Tersoff_PRB_39_5566_Si_C"

#: Screened Tersoff Si-C (Pastewka et al., PRB 87, 205410, 2013)
#: Use with TersoffScr — same underlying string, different potential class.
Tersoff_PRB_39_5566_Si_C__Scr = "Tersoff_PRB_39_5566_Si_C"

#: Goumri-Said Al-N, Chem. Phys. 302, 135 (2004)
Goumri_Said_ChemPhys_302_135_Al_N = "Goumri_Said_ChemPhys_302_135_Al_N"

#: Matsunaga B-C-N, Jpn. J. Appl. Phys. 39, 48 (2000)
Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N = (
    "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N"
)

#: Screened Matsunaga B-C-N — use with TersoffScr.
Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr = (
    "Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N__Scr"
)

# ---------------------------------------------------------------------------
# Brenner potential parameter sets
# ---------------------------------------------------------------------------

#: Erhart-Albe Si-C, Phys. Rev. B 71, 035211 (2005)
Erhart_PRB_71_035211_SiC = "Erhart_PRB_71_035211_SiC"

#: Screened Erhart-Albe Si-C — use with BrennerScr.
Erhart_PRB_71_035211_SiC__Scr = "Erhart_PRB_71_035211_SiC"

#: Albe Pt-C, Phys. Rev. B 65, 195124 (2002)
Albe_PRB_65_195124_PtC = "Albe_PRB_65_195124_PtC"

#: Henriksson Fe-C, Phys. Rev. B 79, 144107 (2009)
Henriksson_PRB_79_144107_FeC = "Henriksson_PRB_79_144107_FeC"

#: Kioseoglou Al-N, Phys. Stat. Sol. (b) 245, 1118 (2008)
Kioseoglou_PSSb_245_1118_AlN = "Kioseoglou_PSSb_245_1118_AlN"

#: Original Brenner C, Phys. Rev. B 42, 9458 (1990), parameter set I
Brenner_PRB_42_9458_C_I = "Brenner_PRB_42_9458_C_I"

#: Original Brenner C, Phys. Rev. B 42, 9458 (1990), parameter set II
Brenner_PRB_42_9458_C_II = "Brenner_PRB_42_9458_C_II"

# ---------------------------------------------------------------------------
# Kumagai potential parameter sets
# ---------------------------------------------------------------------------

#: Kumagai Si, Comp. Mater. Sci. 39, 457 (2007)
Kumagai_CompMaterSci_39_457_Si = "Kumagai_CompMaterSci_39_457_Si"

#: Screened Kumagai Si — use with KumagaiScr.
Kumagai_CompMaterSci_39_457_Si__Scr = "Kumagai_CompMaterSci_39_457_Si"

# ---------------------------------------------------------------------------
# Juslin potential parameter sets
# ---------------------------------------------------------------------------

#: Juslin W-C-H, J. Appl. Phys. 98, 123520 (2005)
Juslin_JAP_98_123520_WCH = "Juslin_JAP_98_123520_WCH"
