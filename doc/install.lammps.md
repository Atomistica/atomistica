_NOTE_: If you are using LAMMPS < 07Sep16, you will need to replace `list->ghost`
by `list->ghostflag` in `pair_atomistica.cpp`.

To compile the LAMMPS-Atomistica interface do the following:

1. Extract the current version of Atomistica to a directory _outside_ of your
   LAMMPS source directory. Alternative: clone the git-repository directly via

     ```
     git clone https://github.com/Atomistica/atomistica.git
     ```

   Cloning the git-repository will download the latest development version, not
   the latest release. It may be safer to download the release tarball.

   This should give the following directory hierachy

     ```
     atomistica/
        build/
        build_lammps/
        build_standalone/
        build_unittests/
        ...
     ```

2. Change to the `atomistica/build_lammps` directory. This directory contains
   makefiles for compilation with the GNU compiler suite (`Makefile.gnu`),
   for compilation with the Intel compiler suite (`Makefile.intel`) and for
   compilation with the XL compiler suite on Blue Gene (`Makefile.xl`). It is
   recommended to use the same compiler used for compiling LAMMPS. Copy one of
   the `Makefile.*` to `Makefile` and edit if necessary. Then compile with

    ```
    make lammps_factories
    make atomistica
    ```

   The final output of this procedure is the file `libatomistica.a`.

2. Copy the files `pair_atomistica.cpp` and `pair_atomistica.h` from the
   `atomistica/src/lammps/pair_style` directory to LAMMPS' `src` directory.

3. Modify your LAMMPS Makefile to link the main executable with
   `libatomistica.a` created by step 2. Add

     ```
     -L/path/to/atomistica/build_lammps -latomistica
     ```

   to the link flags of LAMMPS.

   `pair_atomistica.h` additionally needs
   to be able to find the files `ptrdict.h` and `potentials_factory_c.h` during
   from the atomistica source directory during compilation. These files are
   located here:

     ```
     atomistica/src/support/ptrdict.h 
     atomistica/build_lammps/potentials_factory_c.h
     ```

   Modify your LAMMPS Makefile to add these paths to the header search path by
   adding

     ```
     -I/path/to/atomistica/src/support -I/path/to/atomistica/build_lammps
     ```

   to the compile flags of LAMMPS.

   Atomistica additionally needs to be linked with some LAPACK implementation,
   the Fortran runtime and the MPI library, the details depend on the compiler
   and MPI version.
   
   The following instructions are for _Intel compilers_ only. To link with LAPACK,
   add

     ```
     -mkl=sequential
     ```

   to the link flags of LAMMPS. It will additionally be necessary to link to a
   Fortran runtime library. Failing to link to a Fortran runtime will result in
   undefined symbols when building LAMMPS. For the Intel compiler suite add
   
     ```
     -lifcore -lifport
     ```
     
   to the link flags of LAMMPS. Finally, it will be necessary to link with the 
   Fortran MPI bindings. Try adding
   
     ```
     -l<f77_mpi_library>
     ```
     
   to the link line of LAMMPS. Common names for `<f77_mpi_library>` are
   `mpi_f77` or `mpi_mpifh`; this depends on the MPI implementation that you
   are using.

   The `atomistica/src/lammps/MAKE` subdirectory of the atomistica source
   packages contains a sample makefile, but it is likely outdated. It is
   recommended to modify your LAMMPS makefile by hand!

4. Compile LAMMPS. 

Atomistica potentials are available via the "atomistica" pair style. A typical
LAMMPS script looks like

  ```
  ...  
  units metal  
  ...  
  pair_style atomistica Tersoff   # Tersoff potential  
  pair_coeff * * Si               # LAMMPS atom id 1 is Silicon  
  ...
  ```

If you need to modify the potential parameters create an Atomistica parameter file.

  ```
  Tersoff {  
    A = "1000, 1000, 1000";  
  };
  ```

In the LAMMPS control file use:

  ```
  pair_style atomistica Tersoff params.dat  
  pair_coeff * * Si
  ```  

where "params.dat" is the above file.


# Troubleshooting

## Undefined symbols error:
If, during LAMMPS compilation, you get an error which looks like this:
```
...
./liblammps_mpi.so: error: undefined reference to 'atomistica_startup'
./liblammps_mpi.so: error: undefined reference to 'particles_set_element'
./liblammps_mpi.so: error: undefined reference to 'neighbors_set_pointers'
./liblammps_mpi.so: error: undefined reference to 'ptrdict_read'
...

```
It is possible that you have linked the libraries in an incorrect order. As of 05/12/2023, for a build using ```make mpi``` with GCC and mode = shared, the following lines lead to correct compilation:

```
LMPLINK= -L. -llammps_$@ -L/path/to/atomistica/build_lammps -latomistica -L/path/to/OpenBLAS/0.3.18-GCC-11.2.0/lib/ -lgomp -lgfortran -lstdc++ -lopenblas -lmpi_mpifh

```
in Makefile, and

```
CCFLAGS =	-g -O3 -std=c++11 -I/path/to/atomistica/src/support -I/path/to/atomistica/build_lammps

```
```
LIB = -L/path/to/atomistica/build_lammps -latomistica -L/path/to/OpenBLAS/0.3.18-GCC-11.2.0/lib/ -lgomp -lgfortran -lstdc++ -lopenblas -lmpi_mpifh
```
in makefile.mpi.

