#
# *** Paths
# 
# MDCORE
#
SRCDIR   = ../src


#
# *** Compilers and parallelization
#
# Serial / OpenMP execution (GNU)
#
FC        = gfortran
F90C      = gfortran
CC        = gcc
CXX       = g++
LD        = gfortran
#
# OpenMP parallelization or hybrid MPI/OpenMP
#
OMP_FLAGS = 
#OMP_FLAGS = -fopenmp


#
# *** LAPACK and BLAS link options here.
#
# cygwin lapack/blas
#
EXTRA_LIB += -llapack -lblas -lstdc++

#
# *** Other settings that rarely need to be touched
#
LIBTOOL  = ar r
#
# * Optimization
#
# Normal (GNU)
#
OPTFLAGS = -O3 -funroll-loops
#
# Debug (GNU)
#
#OPTFLAGS = -g -O0 -fbounds-check

#
# * Defines
#
#   -DMDCORE_MONOLITHIC        Compile modules that work only in the standalone
#                              version of MDCORE (i.e. HarmonicHook, etc.)
#   -DMDCORE_PYTHON            Compile Python specific stuff
#   -DLAMMPS                   Compute LAMMPS specific stuff
#
#   MDCORE_MONOLITHIC, MDCORE_PYTHON, and LAMMPS are (should be) mutually exclusive
#
#   -DHAVE_NETCDF              Compile with NetCDF output module
#   -DHAVE_FFTW3               Compile PME module using FFTW3
#   -DHAVE_MKL                 LAPACK implementation is the MKL
#                              (switches printing of MKL version information)
#   -DHAVE_IFPORT              Compiler is Intel Fortran and the ifport module
#                              is present (switches writing of a restart file
#                              upon SIGTERM, i.e. if wallclock time is reached)
#   -DBROKEN_ISO_C_BINDING     c_loc implementation in iso_c_binding is 
#                              broken (basically all gfortran versions)
#   -DHAVE_CUDA                CUDA is available on the system. Compile code to
#                              use CUDA GPU hardware.
# 
# * libAtoms defines
#
#   -DGETENV_F2003             Fortran 2003 getenv is present (define if you
#                              get undefined references to _getenv_)
#   -DGETARG_F2003             Fortran 2003 getarg is present (define if you
#                              get undefined references to _getarg_)
#
#
#   -DQUIP_ARCH=\"MDCORE\"     libAtoms/QUIP internal versioning
#   -DSIZEOF_FORTRAN_T=8       for libAtoms/QUIP C interoperability
#
# - Would be nice to have all explained eventually.
#
DEFINES  = \
	-DLAMMPS \
	-DNO_BIND_C_OPTIONAL


#
# *** Compilation and linking flags
#     (settings should be made mainly above, not here)
#
GFFLAGS  = \
	$(DEFINES) \
	$(OPTFLAGS) \
	$(MPI_FLAGS) \
	$(OMP_FLAGS)
#
# GNU
#
FFLAGS   = $(GFFLAGS) -x f77-cpp-input
F90FLAGS = $(GFFLAGS) $(EXTRA_INCLUDE) \
	-ffree-form -ffree-line-length-none -x f95-cpp-input
CFLAGS   = -O0

#
# Use LDFLAGS = -static if you want a static binary
#
LDFLAGS  =
LIBS     = $(EXTRA_LIB)

include $(SRCDIR)/makefile.inc

