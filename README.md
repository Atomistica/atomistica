Atomistica
==========

Atomistica is a library of interatomic potentials written in Fortran. It is
designed to be used with different computational engines. Currently supported
are ASE and LAMMPS (see below). A list of supported interatomic potentials can
be found [here](POTENTIALS.md).


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
