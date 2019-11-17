Installation instruction for the standalone code
================================================

IMPORTANT: The standalone code is as of now experimental and significantly
less stable than either the Python or LAMMPS interface.

To compile the Atomistica standalone code (MDCORE) do the following:

1. In the build_standalone directory copy one of the Makefile.* to Makefile and
   edit appropriately. Example makefiles for the Intel and GNU compilers are
   provided. Then type

        make factories

   This will create a set of *_factory_* and *_dispatch.f90" files. These files
   are responsible for instantianting and calling instances of potentials,
   Coulomb solvers, integrators and callables. (Callables are the equivalent to
   LAMMPS fixes.)


2. Compile the code with

        make mdcore
