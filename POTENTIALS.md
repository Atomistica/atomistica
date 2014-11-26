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
    ```python
    from atomistica import Tersoff, Tersoff_PRB_39_5566_Si_C
    calc = Tersoff(**Tersoff_PRB_39_5566_Si_C)
    ```

Potentials that employ Brenner's functional form:

*   Paul Erhart, Karsten Albe  
    "Analytical potential for atomistic simulations of silicon, carbon, and silicon carbide"  
    Phys. Rev. B 71, 035211 (2005) - http://dx.doi.org/10.1103/PhysRevB.71.035211  
    ```python
    from atomistica import Brenner
    calc = Brenner(**Erhart_PRB_71_035211_SiC)
    ```

Potentials that employ the REBO2 functional form:

*   Donald W. Brenner, Olga A. Shenderova, Judith A. Harrison, Steven J. Stuart, Boris Ni, Susan B. Sinnott   
    "A second-generation reactive empirical bond order (REBO) potential energy expression for hydrocarbons"   
    J. Phys.: Condens. Matter 14, 783 (2002) - http://dx.doi.org/10.1088/0953-8984/14/4/312  
    ```python
    from atomistica import Rebo2
    calc = Rebo2()
    ````

Special potentials:

*   T. Kumagai, S. Izumi, S. Hara, S. Sakai  
    "Development of bond-order potentials that can reproduce the elastic constants and melting point of silicon for classical molecular dynamics simulation"  
    Comp. Mater. Sci. 39, 457 (2007) - http://dx.doi.org/10.1016/j.commatsci.2006.07.013  
    ```python
    from atomistica import Kumagai, Kumagai_CompMaterSci_39_457_Si
    calc = Kumagai(**Kumagai_CompMaterSci_39_457_Si)
    ```

Screened empirical bond-order potentials:

*   Lars Pastewka, Pablo Pou, Ruben Perez, Peter Gumbsch, Michael Moseler   
    "Describing bond-breaking processes by reactive potentials: Importance of an environment-dependent interaction range"  
    Phys. Rev. B 78, 161402(R) (2008) - http://dx.doi.org/10.1103/PhysRevB.78.161402  
    ```python
    from atomistica import Rebo2Scr
    calc = Rebo2Scr()
    ```

*   Lars Pastewka, Andreas Klemenz, Peter Gumbsch, Michael Moseler  
    "Screened empirical bond-order potential for Si-C"  
    Phys. Rev. B 87, 205410 (2013) - http://dx.doi.org/10.1103/PhysRevB.87.205410  
    arXiv:1301.2142 - http://arxiv.org/abs/1301.2142  
    This paper describes screened version of the Tersoff, Erhart & Albe and Kumagai
    potentials.
    ```python
    from atomistica import TersoffScr, Tersoff_PRB_39_5566_Si_C__Scr
    calc = TersoffScr(**Tersoff_PRB_39_5566_Si_C__Scr)
    from atomistica import BrennerScr, Erhart_PRB_71_035211_SiC__Scr
    calc = BrennerScr(**Erhart_PRB_71_035211_SiC__Scr)
    from atomistica import KumagaiScr, Kumagai_CompMaterSci_39_457_Si__Scr
    calc = KumagaiScr(**Kumagai_CompMaterSci_39_457_Si__Scr)
    ```

*   A general overview on bond-order potentials can be found in   
    Lars Pastewka, Matous Mrovec, Michael Moseler, Peter Gumbsch   
    "Bond order potentials for fracture, wear, and plasticity"   
    MRS Bulletin 37, 493 (2012) - http://dx.doi.org/10.1557/mrs.2012.94

Embedded-atom method potentials:

*   General EAM potentials tabulated in the DYNAMO 'setfl' format.  
    See http://www.ctcms.nist.gov/~cbecker/ for a database of potential datafiles.

