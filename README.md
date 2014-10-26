Atomistica
==========

Atomistica is a library of interatomic potentials written in Fortran. It is
designed to be used with different computational engines. Currently supported
are ASE and LAMMPS (see below).


Supported interatomic potentials
--------------------------------

Atomistica implements the potentials that are documented in the following
publications. Others that employ identical functional forms can be used by
adjusting the parameter set of one of the below.

Empirical bond-order potentials:

*   J. Tersoff  
    "Modeling solid-state chemistry: Interatomic potentials for multicomponent systems"  
    Phys. Rev. B 39, 5566 (1989) - http://dx.doi.org/10.1103/PhysRevB.39.5566

*   Donald W. Brenner, Olga A. Shenderova, Judith A. Harrison, Steven J. Stuart, Boris Ni, Susan B. Sinnott   
    "A second-generation reactive empirical bond order (REBO) potential energy expression for hydrocarbons"   
    J. Phys.: Condens. Matter 14, 783 (2002) - http://dx.doi.org/10.1088/0953-8984/14/4/312

*   Paul Erhart, Karsten Albe  
    "Analytical potential for atomistic simulations of silicon, carbon, and silicon carbide"  
    Phys. Rev. B 71, 035211 (2005) - http://dx.doi.org/10.1103/PhysRevB.71.035211

*   T. Kumagai, S. Izumi, S. Hara, S. Sakai  
    "Development of bond-order potentials that can reproduce the elastic constants and melting point of silicon for classical molecular dynamics simulation"  
    Comp. Mater. Sci. 39, 457 (2007) - http://dx.doi.org/10.1016/j.commatsci.2006.07.013

Screened empirical bond-order potentials:

*   Lars Pastewka, Pablo Pou, Ruben Perez, Peter Gumbsch, Michael Moseler   
    "Describing bond-breaking processes by reactive potentials: Importance of an environment-dependent interaction range"   
    Phys. Rev. B 78, 161402(R) (2008) - http://dx.doi.org/10.1103/PhysRevB.78.161402

*   Lars Pastewka, Andreas Klemenz, Peter Gumbsch, Michael Moseler  
    "Screened empirical bond-order potential for Si-C"  
    Phys. Rev. B 87, 205410 (2013) - http://dx.doi.org/10.1103/PhysRevB.87.205410  
    arXiv:1301.2142 - http://arxiv.org/abs/1301.2142

*   A general overview on bond-order potentials can be found in   
    Lars Pastewka, Matous Mrovec, Michael Moseler, Peter Gumbsch   
    "Bond order potentials for fracture, wear, and plasticity"   
    MRS Bulletin 37, 493 (2012) - http://dx.doi.org/10.1557/mrs.2012.94

Embedded-atom method potentials:

*   General EAM potentials tabulated in the DYNAMO 'setfl' format.  
    See http://www.ctcms.nist.gov/~cbecker/ for a database of potential datafiles.


Usage
-----

Atomistica can be used in two distinct manners. It is recommended to compile the
ASE interface first and run the tests in the "tests" subdirectory. (See ASE
intructions below.)

Currently supported engines are...

1.  ...the Atomistic Simulation Environment   
    (ASE - see https://wiki.fysik.dtu.dk/ase/)
    * Build instructions are in README.PYTHON
    * Examples are in examples/ASE
    * Tests are in tests

2.  ...the Large-scale Atomic/Molecular Massively Parallel Simulator   
    (LAMMPS - see https://lammps.sandia.gov/)
    * Build instructions are in README.LAMMPS
    * Examples are in examples/LAMMPS

You need the following packages:

* Python 2.4.0 or greater (even if you do not compile the Python interface)


Known issues
------------

We have experienced issues with icc/ifort 11.1 20100414 and do not recommend
usage of this particular compiler.


Contact
-------

This software is mainly developed at [Fraunhofer IWM](http://www.en.iwm.fraunhofer.de/business-units/tribology/multiscale-modeling-and-tribosimulation/) and [Karlsruhe Institute of Technology](http://www.iam.kit.edu/zbs/english/). Please write to Lars Pastewka (lars.pastewka@kit.edu) for questions and suggestions. A detailed list of contributors can be found [here](AUTHORS.md).
