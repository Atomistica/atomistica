!! ======================================================================
!! MDCORE - Interatomic potential library
!! https://github.com/pastewka/mdcore
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================
!>
!! The main mdcore module
!!
!! The main mdcore module
!<
module mdcore
#ifdef HAVE_IFPORT
  use ifport
#endif

  use libAtoms_module

  use logging
  use timer
  use tls

  use versioninfo

!!$  use particles
!!$  use neighbors
!!$  use dynamics
!!$  use molecules
!!$
!!$  use bt_fast_multipole
!!$  use charge_overlap
!!$  use cutoff_coulomb
!!$  use direct_coulomb
!!$  use ewald
!!$  use voltage
!!$
!!$  use albe_bop
!!$  use albe_bop_scr
!!$  use born_mayer
!!$  use exp_bop
!!$  use gupta
!!$  use o2_bop
!!$  use o2_bop_scr
!!$  use rebo2
!!$  use rebo2_scr
!!$  use sliding_p
!!$  use tabulated_eam
!!$  use tersoff
!!$  use variable_charge
!!$  use zbl
!!$  use bmhsw
!!$
!!$  use dense_hamiltonian_type
!!$  use dense_hamiltonian
!!$  use dense_solver
!!$
!!$  use dense_notb

  character(3), parameter, private :: month(12) = &
       (/ "Jan", "Feb", "Mar", "Apr", &
       "May", "Jun", "Jul", "Aug", &
       "Sep", "Oct", "Nov", "Dec" /)

  integer, private :: start_time_and_date(8)
  
contains

  !>
  !! Initialize MDCORE
  !!
  !! Open log file, print a greetings message and start the timer.
  !<
  subroutine mdcore_startup(verbosity) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    integer(c_int), value            :: verbosity

    ! ---

    integer                          :: i
    integer                          :: now(8)

#ifdef HAVE_IFPORT
    character(MAX_HOSTNAM_LENGTH+1)  :: runhost
#endif

#ifdef HAVE_MKL
    character(200)                   :: mklversion
#endif

    ! ---

    call system_initialise(verbosity=verbosity)

#ifdef HAVE_IFPORT
    i = hostnam(runhost)
#endif
    call date_and_time(values=now)

    call logging_start("md.log")

#ifdef LAMMPS
    call prscrlog("Welcome to - LAMMPS+MDCORE+libAtoms -")
#else
    call prscrlog("Welcome to - MDCORE+libAtoms -")
#endif
    call prscrlog
    call prscrlog("   MDCORE revision:      " // trim(mdcore_revision))
    call prscrlog("   MDCORE URL:           " // trim(mdcore_url))
    call prscrlog("   libAtoms revision:    " // trim(libAtoms_revision))
    call prscrlog("   libAtoms URL:         " // trim(libAtoms_URL))
#ifdef HAVE_MKL
    call mklgetversionstring(mklversion)
    call prscrlog("   MKL version:          " // trim(mklversion))
#endif
    call prscrlog
    call prscrlog("   architecture:         " // trim(arch))
    call prscrlog
    call prscrlog("   build host:           " // trim(buildhost))
    call prscrlog("   build date:           " // trim(builddate))
    call prscrlog("   compile options:      " // trim(compileroptions))
    call prscrlog("   compiler version:     " // trim(compilerversion))
    call prscrlog
#ifdef HAVE_IFPORT
    call prscrlog("   run host:             " // trim(runhost))
#endif
    call prscrlog("   run date:             " // month(now(2)) // " " // now(3) // " " // now(1) // " " // now(5) // ":" // now(6) // ":" // now(7))
    call prscrlog

#if defined(_OPENMP) && defined(_MPI)
    call prscrlog("   Hybrid MPI+OpenMP:    " // mpi_n_procs() // " processes, " // omp_get_max_threads() // " threads each")
    call prscrlog
#elif defined(_OPENMP)
    call prscrlog("   Using OpenMP:         " // omp_get_max_threads() // " threads")
    call prscrlog
#elif defined(_MPI)
    call prscrlog("   Using MPI:            " // mpi_n_procs() // " processes")
    call prscrlog
#endif

    call timer_start("MDCORE")

    start_time_and_date  = now

  endsubroutine mdcore_startup


  !>
  !! Finalize MDCORE 
  !!
  !! Print timing information, close log file and write an empty "DONE" file.
  !<
  subroutine mdcore_shutdown() bind(C)
    implicit none

    type(InOutput) :: done_file
    integer :: now(8)

    call timer_stop("MDCORE")

    call date_and_time(values=now)

    call prscrlog
    call prscrlog("   simulation started:  " // month(start_time_and_date(2)) // " " // start_time_and_date(3) // " " // start_time_and_date(1) // " " // start_time_and_date(5) // ":" // start_time_and_date(6) // ":" // start_time_and_date(7))
    call prscrlog("   simulation ended:    " // month(now(2)) // " " // now(3) // " " // now(1) // " " // now(5) // ":" // now(6) // ":" // now(7))
    call prscrlog

    call timer_print_to_log

    call logging_stop

    !$omp parallel
    call tls_del
    !$omp end parallel

#ifdef _MPI
    if (mpi_id() == ROOT) then
#endif
    
    !
    ! Create the (empty) file DONE to let everyone know we finished properly
    !

    call initialise(done_file, "DONE", OUTPUT)
    call finalise(done_file)

#ifdef _MPI
    endif
#endif

#ifndef LAMMPS
    ! Avoid calling MPI_Finalize twice
   call system_finalise
#endif

  endsubroutine mdcore_shutdown

endmodule
