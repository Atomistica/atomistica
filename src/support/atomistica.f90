!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 2 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!! ======================================================================
!>
!! The main atomistica module
!!
!! The main atomistica module
!<

module atomistica
#ifdef _OPENMP
  use omp_lib
#endif

#ifdef HAVE_IFPORT
  use ifport
#endif

  use supplib

  use versioninfo

  character(3), parameter, private :: month(12) = &
       (/ "Jan", "Feb", "Mar", "Apr", &
       "May", "Jun", "Jul", "Aug", &
       "Sep", "Oct", "Nov", "Dec" /)

  integer, private :: start_time_and_date(8)
  
#ifdef LAMMPS
#define DEFINE_GIT_IDENT
#include "../src/lammps/pair_style/pair_atomistica.cpp"
  interface get_atomistica_pair_style_git_ident
     subroutine get_atomistica_pair_style_git_ident(ident) bind(C)
        use, intrinsic :: iso_c_binding
        character(C_CHAR) :: ident(*)
     endsubroutine get_atomistica_pair_style_git_ident
  endinterface

#endif

contains

  !>
  !! Initialize Atomistica
  !!
  !! Open log file, print a greetings message and start the timer.
  !<
  subroutine atomistica_startup() bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    ! ---

    integer                          :: i
    integer                          :: now(8)

#ifdef HAVE_IFPORT
    character(MAX_HOSTNAM_LENGTH+1)  :: runhost
#endif

#ifdef HAVE_MKL
    character(200)                   :: mklversion
#endif

#ifdef LAMMPS
    character(C_CHAR), target        :: pair_style_git_ident(1024)
#endif

    ! ---

#ifdef HAVE_IFPORT
    i = hostnam(runhost)
#endif
    call date_and_time(values=now)

    call logging_start("atomistica.log")

#ifdef LAMMPS
    call prscrlog("Welcome to - LAMMPS+Atomistica -")
#else
    call prscrlog("Welcome to - Atomistica -")
#endif
    call prscrlog
    call prscrlog("   Atomistica revision:  " // trim(atomistica_revision))
    call prscrlog("   Atomistica rev date:  " // trim(atomistica_date))
    call prscrlog("   Atomistica URL:       " // trim(atomistica_url))
    call prscrlog
    call prscrlog("   architecture:         " // trim(arch))
    call prscrlog
    call prscrlog("   build host:           " // trim(buildhost))
    call prscrlog("   build date:           " // trim(builddate))
    call prscrlog("   compiler version:     " // trim(compilerversion))
    call prscrlog("   compile options:      " // trim(compileroptions))
#ifdef HAVE_MKL
    call mkl_get_version_string(mklversion)
    call prscrlog("   MKL version:          " // trim(mklversion))
#endif
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

#ifdef LAMMPS
    ! Check if Atomistica library and LAMMPS' Atomistica pair_style versions are identical.
    call get_atomistica_pair_style_git_ident(pair_style_git_ident)
    if (a2s(c_f_string(c_loc(pair_style_git_ident(1)))) /= ATOMISTICA_PAIR_STYLE_GIT_IDENT) then
       stop "Fatal: GIT blobs of pair_atomistica.cpp used when compiling the Atomistica library and LAMMPS do not agree. Please copy the current pair_atomistica.cpp and pair_atomistica.h to your LAMMPS src directory."
    endif
#endif

!    call rng_init(now(7)+1)

    call timer_start("MDCORE")

    start_time_and_date  = now

  endsubroutine atomistica_startup


  !>
  !! Finalize MDCORE 
  !!
  !! Print timing information, close log file.
  !<
  subroutine atomistica_shutdown() bind(C)
    implicit none

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

  endsubroutine atomistica_shutdown

endmodule
