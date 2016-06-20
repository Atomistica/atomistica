#include "macros.inc"

program test_communicator
  use supplib

  use particles
  use neighbors

  use native_io

  use atomistica

  implicit none

  type(MPI_Context) :: mpi

  type(particles_t) :: p 

  integer :: error

  ! ---

  call initialise(mpi)
  call atomistica_startup

  call units_init(eV_A)

  call init(p)

  call read_atoms(p, "aC.dat", error=error)
  HANDLE_ERROR(error)

  call del(p)

  call atomistica_shutdown
  call finalise(mpi)

endprogram test_communicator

