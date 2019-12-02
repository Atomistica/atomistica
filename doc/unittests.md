To compile and run the Fortran unittests do the following:

1. In the build_unittests directory do

    make -f <Makefile> unittests

   Here <Makefile> is some makefile configured for your environment.
   Example makefiles for the Intel and GNU compilers are provided. This should
   generate the file "unittests" (or "unittests.exe" on Windows).

2. Run

    ./unittests
