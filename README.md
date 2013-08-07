Atomistica
==========

Atomistica is a library of interatomic potentials. It currently supports the
potentials described in the following publication:

*   J. Tersoff  
    "Modeling solid-state chemistry: Interatomic potentials for multicomponent systems"  
    Phys. Rev. B 39, 5566 (1989) - http://dx.doi.org/10.1103/PhysRevB.39.5566

*   Paul Erhart, Karsten Albe  
    "Analytical potential for atomistic simulations of silicon, carbon, and silicon carbide"  
    Phys. Rev. B 71, 035211 (2005) - http://dx.doi.org/10.1103/PhysRevB.71.035211

*   T. Kumagai, S. Izumi, S. Hara, S. Sakai  
    "Development of bond-order potentials that can reproduce the elastic constants and melting point of silicon for classical molecular dynamics simulation"  
    Comp. Mater. Sci. 39, 457 (2007) - http://dx.doi.org/10.1016/j.commatsci.2006.07.013

*   Lars Pastewka, Andreas Klemenz, Peter Gumbsch, Michael Moseler  
    "Screened empirical bond-order potential for Si-C"  
    Phys. Rev. B 87, 205410 (2013) - http://dx.doi.org/10.1103/PhysRevB.87.205410  
    arXiv:1301.2142 - http://arxiv.org/abs/1301.2142

Atomistica can be used in two distinct manners:

1.  From Python using the Atomistic Simulation Environment
    (ASE - see https://wiki.fysik.dtu.dk/ase/)
    * Build instructions are in README.PYTHON
    * Examples are in examples/ASE

2.  From within LAMMPS as a separate pair style
    (LAMMPS - see https://lammps.sandia.gov/)
    * Build instructions are in README.LAMMPS
    * Examples are in examples/LAMMPS

You need the following packages:

* Python 2.4.0 or greater (even if you do not compile the Python interface)
  (Some parts of the Python interface could rely on features from Python 2.6.0)

