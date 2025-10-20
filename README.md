![Atomistica](https://github.com/Atomistica/atomistica/blob/master/images/logo.png)

Atomistica is a library of interatomic potentials, including empirical potentials and non-orthogonal tight-binding.
It is designed to be plugged into different simulation environments. We currently support the
[Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/) and
the [Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS)](http://lammps.sandia.gov/).
A list of interatomic potentials can be found [here](doc/potentials.md).

Build status
------------

[![Build Status](https://github.com/Atomistica/atomistica/actions/workflows/run_tests.yml/badge.svg)](https://github.com/Atomistica/atomistica/actions)

Usage
-----

Atomistica can be used in two distinct manners. It is recommended to compile the
ASE interface first and run the tests in the `tests` subdirectory. (See ASE
intructions below.)

Currently supported simulation environments are...

1.  ...the Atomistic Simulation Environment
    (ASE - see https://wiki.fysik.dtu.dk/ase/)
    * Build instructions are in [doc/install.python.md](doc/install.python.md)
    * Examples are in [examples/ASE](examples/ASE)
    * Tests are in tests
    * Atomistica supports Python 3.8+ (Python 3.12+ with NumPy 2.x recommended)

2.  ...the Large-scale Atomic/Molecular Massively Parallel Simulator   
    (LAMMPS - see https://lammps.sandia.gov/)
    * Build instructions are in [doc/install.lammps.md](doc/install.lammps.md)
    * Examples are in [examples/LAMMPS](examples/LAMMPS)

3.  ...MDCORE, the standalone molecular dynamics code of Atomistica.
    Note that there is _no documentation_ for the standalone version.
    * Build instructions are in [doc/install.standalone.md](doc/install.standalone.md)
    * A brief documentation is in [doc/standalone.md](doc/standalone.md)
    * Examples are in [examples/STANDALONE](examples/STANDALONE)

You need the following packages:

* Python 3.8.0 or greater (Python 3.12+ recommended with NumPy 2.x support)
* Meson >= 1.1.0 and meson-python >= 0.15.0 (for Python interface)
* NumPy >= 1.21.0 (NumPy 2.x supported)
* LAPACK library
* A Fortran compiler (gfortran or ifort), C compiler, and C++ compiler


Contact
-------

This software is developed at the
[University of Freiburg](http://www.imtek.uni-freiburg.de/laboratories/simulation)
and
[Fraunhofer IWM](http://www.en.iwm.fraunhofer.de/business-units/tribology/multiscale-modeling-and-tribosimulation/).
Please write to Lars Pastewka (lars.pastewka@imtek.uni-freiburg.de) for questions and suggestions.
A complete list of contributors can be found [here](AUTHORS.md).
