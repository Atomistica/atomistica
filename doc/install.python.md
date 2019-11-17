The following instruction explain how to compile, use and test the Python
interface to Atomistica.

1. Setup your compiler. Edit setup.cfg

   For GNU Fortran/C use:

     [config_fc]
     fcompiler=gnu95
     f90flags=-cpp -ffree-line-length-none [-fopenmp]

   For Intel Fortran/C (intel64) use:

     [config_fc]
     fcompiler=intelem
     f90flags=-fpp [-openmp]
 
     [config_cc]
     compiler=intel

   There may be error messages complaining about undefined symbols when loading
   Atomistica in step 3 below. It is then necessary to force linking to the
   respective libraries. Additional libraries can be specified in the build_ext
   section:

     [build_ext]
     libraries=ifcore,ifport,iomp5,pthread

   This example is the link line typically required to compile with the Intel
   compiler suite and OpenMP enabled. ifcore contains is the Intel Fortran
   runtime and ifport additional portability functions. iomp5 is the Intel
   OpenMP runtime which requires to additionally link to the posix pthread
   library. The respective libraries for the GNU compiler are gfortran and
   gomp.

   You can get a list of available Fortran and C compilers by executing:

     python setup.py build --help-fcompiler
     python setup.py build --help-compiler

   More information can be found here:

     http://docs.python.org/2/install/index.html
     http://thread.gmane.org/gmane.comp.python.numeric.general/15793


2. Compile the Python extension. Execute

     python setup.py build

   This will build and link all Fortran and C sources required for the Python
   interface. You may need to edit setup.py if the LAPACK libraries are not
   automatically detected.

   Note: You will need to

     python setup.py clean
     - or -
     python setup.py build --force

   if the source has changed. Unfortunately, the current numpy distutils
   won't relink the module even if the library has been recompiled without
   --force present.


3. Test if the Python interface imports with errors. Type

     source <path-to-atomistica>/env.sh

   which sets up the environment. To test if Atomistica can be imported into
   Python execute"

     python -c "import atomistica"

   To resolve undefined symbol errors link to the relevant library.
   (See step 1 above.)


4. Test if the Python interface gives correct results. Type

     cd <path-to-atomistica>/tests
     python run_tests.py

   Each test can also be run directly. This will produce more diagnostic output:

     python bulk_properties.py
     python forces_and_virial.py
     python rebo2_molecules.py
