![Atomistica](https://github.com/Atomistica/atomistica/blob/master/images/logo.png)

Atomistica is a library of interatomic potentials. It is designed to be plugged
into different simulation environments. We currently support the
[Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/) and
the [Large-scale Atomic/Molecular Massively Parallel Simulator (LAMMPS)](http://lammps.sandia.gov/).
A list of interatomic potentials can be found [here](POTENTIALS.md).

Build status
------------

[![Build Status](https://travis-ci.org/Atomistica/atomistica.svg?branch=master)](https://travis-ci.org/Atomistica/atomistica)

Usage
-----

Atomistica can be used in two distinct manners. It is recommended to compile the
ASE interface first and run the tests in the `tests` subdirectory. (See ASE
intructions below.)

Currently supported simulation environments are...

1.  ...the Atomistic Simulation Environment   
    (ASE - see https://wiki.fysik.dtu.dk/ase/)
    * Build instructions are in [README.PYTHON.md](README.PYTHON.md)
    * Examples are in [examples/ASE](examples/ASE)
    * Tests are in tests
    * Atomistica supports Python 2 and Python 3

2.  ...the Large-scale Atomic/Molecular Massively Parallel Simulator   
    (LAMMPS - see https://lammps.sandia.gov/)
    * Build instructions are in [README.LAMMPS.md](README.LAMMPS.md)
    * Examples are in [examples/LAMMPS](examples/LAMMPS)

3.  ...MDCORE, the standalone molecular dynamics code of Atomistica.
    Note that there is _no documentation_ for the standalone version.
    * Build instruction are in [README.STANDALONE.md](README.STANDALONE.md)
    * Examples are in [examples/STANDALONE](examples/STANDALONE)

You need the following packages:

* Python 2.6.0 or greater (Python is needed even if you do not compile the
  Python interface to auto-generate parts of the source code)


Contact
-------

This software is developed at the
[University of Freiburg](http://www.imtek.uni-freiburg.de/laboratories/simulation)
and
[Fraunhofer IWM](http://www.en.iwm.fraunhofer.de/business-units/tribology/multiscale-modeling-and-tribosimulation/).
Please write to Lars Pastewka (lars.pastewka@imtek.uni-freiburg.de) for questions and suggestions.
A complete list of contributors can be found [here](AUTHORS.md).
