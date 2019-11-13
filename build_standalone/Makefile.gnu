#
# This Makefile can be used to compiler mdcore standalone on cygwin with the GNU
# compilers.
#


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
MPI_FLAGS =


#
# Initial EXTRA_INCLUDE and EXTRA_LIB
#
# -cxxlib links to the C++ runtime
#
HAVE_NETCDF = 0
EXTRA_FLAGS =
EXTRA_INCLUDE =
EXTRA_LIB = -lstdc++

#
# OpenMP parallelization
#
OMP_FLAGS =
#OMP_FLAGS = -fopenmp


#
# *** Extra includes and libraries
#
# Check for NetCDF
#
ifneq ("$(shell which nf-config)","") 
HAVE_NETCDF = 1
EXTRA_FLAGS += -DHAVE_NETCDF
EXTRA_INCLUDE += $(shell nf-config --fflags)
EXTRA_LIB += $(shell nf-config --flibs)
else
ifneq ("$(shell which nc-config)","")
HAVE_NETCDF = 1
EXTRA_FLAGS += -DHAVE_NETCDF
EXTRA_INCLUDE += $(shell nc-config --fflags)
EXTRA_LIB += $(shell nc-config --flibs)
endif
endif

#
# FIXME!!! Implement FFTW check
#
HAVE_FFTW3 = 0
#EXTRA_FLAGS += -DHAVE_FFTW3
#EXTRA_INCLUDE += -I/j1a/pas/applications/fftw-3.3/include
#EXTRA_LIB += -L/j1a/pas/applications/fftw-3.3/lib -lfftw3


#
# *** LAPACK and BLAS link options here.
#
# lapack/blas
#
EXTRA_LIB += -llapack -lblas

#
## DFT-D3 library
#
#

HAVE_DFTD3 = 0

ifneq ($(HAVE_DFTD3),0)
EXTRA_FLAGS += -DHAVE_DFTD3
DFTD3_PATH = # Please specify the path to dftd3-lib
EXTRA_LIB += -L$(DFTD3_PATH) -ldftd3
EXTRA_INCLUDE += -I$(DFTD3_PATH)
endif


#
# *** Other settings that rarely need to be touched
#
LIBTOOL  = ar r
#
# * Optimization
#
# Normal (GNU)
#
OPTFLAGS = -g -O3 -funroll-loops -fbacktrace
#
# Debug (GNU)
#
#OPTFLAGS = -g -O0 -fbacktrace -fbounds-check

#
# * Defines
#
#   -DMDCORE_PYTHON            Compile Python specific stuff
#
#   -DHAVE_NETCDF              Compile with NetCDF output module
#   -DHAVE_FFTW3               Compile PME module using FFTW3
#   -DHAVE_MKL                 LAPACK implementation is the MKL
#                              (switches printing of MKL version information)
#   -DHAVE_IFPORT              Compiler is Intel Fortran and the ifport module
#                              is present (switches writing of a restart file
#                              upon SIGTERM, i.e. if wallclock time is reached)
#   -DHAVE_CUDA                CUDA is available on the system. Compile code to
#                              use CUDA GPU hardware.
#   -DHAVE_DFTD3               DFT-D3 library
#
# 
DEFINES  = \
	-DNO_BIND_C_OPTIONAL

#
# *** Compilation and linking flags
#     (settings should be made mainly above, not here)
#
GFFLAGS  = \
	$(DEFINES) \
	$(OPTFLAGS) \
	$(OMP_FLAGS) \
	$(EXTRA_FLAGS)

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

