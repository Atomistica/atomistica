Interatomic potentials
======================

This file contains a list of implemented interatomic potentials and parameters
set provided. The code fragments show Python code necessary to instantiate an
ASE calculator object for the respective potential/parameterization. Note that
for all potentials listed below there is a default parameter set that is used
if instatiated without and explicit database dictionary. This default is the 
first potential listed below.

Empirical bond-order potentials
-------------------------------

Potentials that employ the Tersoff functional form:

*   J. Tersoff  
    "Modeling solid-state chemistry: Interatomic potentials for multicomponent systems"  
    Phys. Rev. B 39, 5566 (1989) - http://dx.doi.org/10.1103/PhysRevB.39.5566  

        from atomistica import Tersoff, Tersoff_PRB_39_5566_Si_C  
        calc = Tersoff(**Tersoff_PRB_39_5566_Si_C)

*   S. Goumri-Said, M.B. Kanoun, A.E. Merad, G. Merad, H. Aourag  
    "Prediction of structural and thermodynamic properties of zinc-blende AlN: molecular dynamics simulation"
    Chem. Phys. 302, 135 (2004) - http://dx.doi.org/10.1016/j.chemphys.2004.03.030

        from atomistica import Tersoff, Goumri_Said_ChemPhys_302_135_Al_N  
        calc = Tersoff(**Goumri_Said_ChemPhys_302_135_Al_N)

*   Katsuyuki Matsunaga, Craig Fisher, Hideaki Matsubara  
    "Tersoff potential parameters for simulating cubic boron carbonitrides"  
    Jpn. J. Appl. Phys. 39, 48 (2000) - http://dx.doi.org/10.1143/JJAP.39.L48  

        from atomistica import Tersoff, Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N  
        calc = Tersoff(**Matsunaga_Fisher_Matsubara_Jpn_J_Appl_Phys_39_48_B_C_N)

Potentials that employ Brenner's functional form:

*   Paul Erhart, Karsten Albe  
    "Analytical potential for atomistic simulations of silicon, carbon, and silicon carbide"  
    Phys. Rev. B 71, 035211 (2005) - http://dx.doi.org/10.1103/PhysRevB.71.035211  

        from atomistica import Brenner, Erhart_PRB_71_035211_SiC  
        calc = Brenner(**Erhart_PRB_71_035211_SiC)

*   Donald Brenner  
    "Empirical potential for hydrocarbons for use in simulating the chemical vapor deposition of diamond films"  
    Phys. Rev. B 42, 9458 (1990) - http://dx.doi.org/Phys. Rev. B 42, 9458 (1990)  
    Warning: The original formulation uses look-up tables to correct the
    energies for radicals and hydrocarbons. This potential is implemented
    without these lookup tables.  
    Note: The paper describes two parameter sets, denoted as I and II.  
    Parameter set I:

        from atomistica import Brenner, Brenner_PRB_42_9458_C_I  
        calc = Brenner(**Brenner_PRB_42_9458_C_I)

    Parameter set II:

        from atomistica import Brenner, Brenner_PRB_42_9458_C_II  
        calc = Brenner(**Brenner_PRB_42_9458_C_II)

*   Karsten Albe, Kai Nordlund, Robert S. Averback  
    "Analytical bond-order potential for platinum-carbon"  
    Phys. Rev. B 65, 195124 (2002) - http://dx.doi.org/10.1103/PhysRevB.65.195124  

        from atomistica import Brenner, Albe_PRB_65_195124_PtC  
        calc = Brenner(**Albe_PRB_65_195124_PtC)

*   K.O.E. Henriksson, Kai Nordlund  
    "Simulations of cementite: An analytical potential for the Fe-C system"  
    Phys. Rev. B 79, 144107 (2009) - http://dx.doi.org/10.1103/PhysRevB.79.144107  

        from atomistica import Brenner, Henriksson_PRB_79_114107_FeC  
        calc = Brenner(**Henriksson_PRB_79_114107_FeC)

*   J. Kioseoglou, Ph. Komninou, Th. Karakostas  
    "Interatomic potential calculations of III (Al, In)-N planar defects with a III‐species environment approach"  
    Phys. Stat. Sol. (b) 245, 1118 (2008) - http://dx.doi.org/10.1002/pssb.200844122  

        from atomistica import Brenner, Kioseoglou_PSSb_245_1118_AlN  
        calc = Brenner(**Kioseoglou_PSSb_245_1118_AlN)

Potentials that employ the REBO2 functional form:

*   Donald W. Brenner, Olga A. Shenderova, Judith A. Harrison, Steven J. Stuart, Boris Ni, Susan B. Sinnott   
    "A second-generation reactive empirical bond order (REBO) potential energy expression for hydrocarbons"   
    J. Phys.: Condens. Matter 14, 783 (2002) - http://dx.doi.org/10.1088/0953-8984/14/4/312  

        from atomistica import Rebo2  
        calc = Rebo2()

Special potentials:

*   N. Juslin, Paul Erhart, P. Traskelin, J. Nord, Krister O.E. Henriksson, Kai Nordlund, E. Salonen, Karsten Albe  
    "Analytical interatomic potential for modeling nonequilibrium processes in the W–C–H system"  
    J. Appl. Phys. 98, 123520 (2005) - http://dx.doi.org/10.1063/1.2149492  

        from atomistica import Juslin, Juslin_JAP_98_123520_WCH  
        calc = Juslin(**Juslin_JAP_98_123520_WCH)

*   T. Kumagai, S. Izumi, S. Hara, S. Sakai  
    "Development of bond-order potentials that can reproduce the elastic constants and melting point of silicon for classical molecular dynamics simulation"  
    Comp. Mater. Sci. 39, 457 (2007) - http://dx.doi.org/10.1016/j.commatsci.2006.07.013  

        from atomistica import Kumagai, Kumagai_CompMaterSci_39_457_Si  
        calc = Kumagai(**Kumagai_CompMaterSci_39_457_Si)

Screened empirical bond-order potentials:

*   Lars Pastewka, Pablo Pou, Ruben Perez, Peter Gumbsch, Michael Moseler   
    "Describing bond-breaking processes by reactive potentials: Importance of an environment-dependent interaction range"  
    Phys. Rev. B 78, 161402(R) (2008) - http://dx.doi.org/10.1103/PhysRevB.78.161402  

        from atomistica import Rebo2Scr  
        calc = Rebo2Scr()

*   Lars Pastewka, Andreas Klemenz, Peter Gumbsch, Michael Moseler  
    "Screened empirical bond-order potential for Si-C"  
    Phys. Rev. B 87, 205410 (2013) - http://dx.doi.org/10.1103/PhysRevB.87.205410  
    arXiv:1301.2142 - http://arxiv.org/abs/1301.2142  
    This paper describes screened versions of the Tersoff, Erhart & Albe and
    Kumagai potentials.

        from atomistica import TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr  
        calc = TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr)  
        from atomistica import BrennerScr, Erhart_PRB_71_035211_SiC__Scr  
        calc = BrennerScr(**Erhart_PRB_71_035211_SiC__Scr)  
        from atomistica import KumagaiScr, Kumagai_CompMaterSci_39_457_Si__Scr  
        calc = KumagaiScr(**Kumagai_CompMaterSci_39_457_Si__Scr)

*   A general overview on bond-order potentials can be found in   
    Lars Pastewka, Matous Mrovec, Michael Moseler, Peter Gumbsch   
    "Bond order potentials for fracture, wear, and plasticity"   
    MRS Bulletin 37, 493 (2012) - http://dx.doi.org/10.1557/mrs.2012.94

Embedded-atom method potentials:

*   General EAM potentials tabulated in the DYNAMO 'setfl' format.  
    See http://www.ctcms.nist.gov/~cbecker/ for a database of potential
    datafiles.

