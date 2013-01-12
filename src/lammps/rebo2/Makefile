.SUFFIXES: .c .h .f .f90

###

MDCORE_SRC = ../../..

###

#
# ifort options
#
FC        = ifort
F90C      = ifort
DEFINES   = -DSILENT -DQUIP_ARCH=\"LAMMPS\" -DGETENV_F2003 -DLAMMPS
OPTFLAGS  = -m64 -fpconstant -xHost -O3 -ip -funroll-loops -unroll-aggressive -no-prec-div -no-prec-sqrt 
CFLAGS    = -g -O0 -fPIC
F90FLAGS  = -g -fpp -warn unused -fPIC -traceback $(OPTFLAGS) $(DEFINES)

###

CSRC = \
	cutil.c

F90SRC = \
	error.f90 \
	System.f90 \
	Units.f90 \
	libAtoms_utils_no_module.f90 \
	libAtoms.f90 \
	c_f.f90 \
	tls.f90 \
	io.f90 \
	logging.f90 \
	simple_spline.f90 \
	table2d.f90 \
	table3d.f90 \
	rebo2_default_tables.f90 \
	rebo2.f90 \
	rebo2_scr.f90

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $^ -o $@

all: rebo2

rebo2: $(F90SRC:.f90=.o) $(CSRC:.c=.o)
	ar cr librebo2.a $(F90SRC:.f90=.o) $(CSRC:.c=.o)
	echo "done"

copy:
	cp $(MDCORE_SRC)/src/libAtoms/cutil.c .
	cp $(MDCORE_SRC)/src/libAtoms/error.inc .
	cp $(MDCORE_SRC)/src/libAtoms/libatoms.h .
	cp $(MDCORE_SRC)/src/filter.inc .
	cp $(MDCORE_SRC)/src/macros.inc .
	cp $(MDCORE_SRC)/src/spline.inc .
	cp $(MDCORE_SRC)/src/libAtoms/error.f90 error.f90
	cp $(MDCORE_SRC)/src/libAtoms/System.f90 System.f90
	cp $(MDCORE_SRC)/src/libAtoms/Units.f90 Units.f90
	cp $(MDCORE_SRC)/src/support/c_f.f90 .
	cp $(MDCORE_SRC)/src/support/io.f90 .
	cp $(MDCORE_SRC)/src/support/logging.f90 .
	cp $(MDCORE_SRC)/src/support/simple_spline.f90 .
	cp $(MDCORE_SRC)/src/support/table2d.f90 .
	cp $(MDCORE_SRC)/src/support/table3d.f90 .
	cp $(MDCORE_SRC)/src/support/tls.f90 .
	cp $(MDCORE_SRC)/src/potentials/bop/rebo2/bop_kernel_rebo2.f90 .
	cp $(MDCORE_SRC)/src/potentials/bop/rebo2/rebo2_db.f90 .
	cp $(MDCORE_SRC)/src/potentials/bop/rebo2/rebo2_default_tables.f90 .
	cp $(MDCORE_SRC)/src/potentials/bop/rebo2/rebo2_func.f90 .
	cp $(MDCORE_SRC)/src/potentials/bop/rebo2/rebo2_type.f90 .

clean:
	rm -f *.o *.mod

distclean:
	rm cutil.c
	rm error.inc
	rm filter.inc
	rm macros.inc
	rm spline.inc
	rm error.f90
	rm System.f90
	rm Units.f90
	rm io.f90
	rm logging.f90
	rm simple_spline.f90
	rm table2d.f90
	rm table3d.f90
	rm tls.f90
	rm bop_kernel_rebo2.f90
	rm rebo2_db.f90
	rm rebo2_default_tables.f90
	rm rebo2_func.f90
	rm rebo2_type.f90
