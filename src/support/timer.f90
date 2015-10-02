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
!! Nested timing routines
!<

#include "macros.inc"

module timer
  use c_f
  use error_module
  use system_module
  use mpi_context_module
  use logging

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  save

  private

  integer, parameter  :: screen = 6
  integer, parameter  :: TIMER_MAX = 200
  integer, parameter  :: TIMER_STR_LEN = 50

  integer             :: nr = 0
  integer             :: current = 0

  character(len=TIMER_STR_LEN)  :: names(TIMER_MAX)
  real(DP)                      :: t(TIMER_MAX, 2)     = 0.0_DP
  real(DP)                      :: times(TIMER_MAX)    = 0.0_DP
  integer                       :: calls(TIMER_MAX)    = 0
  integer                       :: parents(TIMER_MAX)  = 0

#ifdef _OPENMP
  real(DP)                      :: wt(TIMER_MAX)      = 0.0_DP
  real(DP)                      :: wtimes(TIMER_MAX)  = 0.0_DP
#endif

  public  :: timer_start, timer_stop, timer_print, timer_print_to_log

contains

#undef MKL

  !>
  !! Start a timer with a given name
  !<
  subroutine timer_start(name, error)
#ifdef MKL
    use ifport
#endif

    implicit none

    character(*),      intent(in)   :: name
    integer, optional, intent(out)  :: error

    ! ---

    logical  :: found
    integer  :: It(8), p
    real(DP) :: t_now
#ifdef MKL
    real(4)  :: h(2)
#endif

    ! ---

    INIT_ERROR(error)

#ifdef DEBUG
    call print("timer_start: " // name, PRINT_ALWAYS)
#endif

#ifdef MKL
    t_now = etime(h(:))
#else
    call date_and_time(values = It)
    t_now = It(5)*3600.0_DP+It(6)*60.0_DP+It(7)+It(8)/1000.0_DP
#endif

    found = .false.
    do p = 1, nr
       if (trim(name) == trim(names(p)) .and. current == parents(p)) then
          found = .true.
          exit
       end if
    end do

    if (.not. found) then
       nr = nr + 1

       if (nr > TIMER_MAX) then
          RAISE_ERROR("Number of timers exceeded with time " // name // ".", error)
       endif

       t(nr,1)      = t_now            ! start time...
       t(nr,2)      = real( It(3) )    !       ...and day
       times(nr)    = 0.0_DP

#ifdef _OPENMP
       wt(nr)       = omp_get_wtime()
       wtimes(nr)   = 0.0_DP
#endif

       names(nr)    = name
       parents(nr)  = current
       current      = nr

       p = nr
    else
       t(p,1)   = t_now
       t(p,2)   = real(It(3))

#ifdef _OPENMP
       wt(p)    = omp_get_wtime()
#endif

       current  = p
    endif
    calls(p) = calls(p) + 1
    
  endsubroutine timer_start
   
  
  !>
  !! Stop a timer with a given name
  !<
  subroutine timer_stop(name, error)
#ifdef MKL
    use ifport
#endif

    implicit none

    character(*),      intent(in)   :: name
    integer, optional, intent(out)  :: error
    
    ! ---
 
    integer  :: It(8), p
    real(DP) :: t_now, days
#ifdef MKL
    real(4)  :: h(2)
#endif

    ! ---

    INIT_ERROR(error)

#ifdef DEBUG
    call print("timer_stop: " // name, PRINT_ALWAYS)
#endif

#ifdef MKL 
    t_now = etime(h(:))
#else
    call date_and_time(values = It)
    t_now = It(5)*3600.0_DP+It(6)*60.0_DP+It(7)+It(8)/1000.0_DP
#endif

    p = current
    if (trim(name) /= trim(names(p))) then
       times(p) = 0
       RAISE_ERROR("Warning: Timer '" // name // "' not current. Current timer: '" // trim(names(p)) // "'", error)
    else
#ifdef MKL
       times(p)   = times(p) + ( t_now - t(p,1) )
#else
       days       = It(3) - t(p,2)
       times(p)   = times(p) + 86400d0*days + ( t_now - t(p,1) )
#endif

#ifdef _OPENMP
       wtimes(p)  = wtimes(p) + ( omp_get_wtime() - wt(p) )
#endif

       current    = parents(p)
    endif
    
  endsubroutine timer_stop


  !>
  !! Print timings to screen
  !<
  recursive subroutine timer_print_for_parent(un, p, shift)
    implicit none

    integer, intent(in)   :: un
    integer, intent(in)   :: p
    integer, intent(in)   :: shift

    ! ---

    integer, parameter  :: NAME_STR_LEN = 60

    ! ---

    character(NAME_STR_LEN)  :: name

    integer                  :: i
    logical                  :: first
    real(DP)                 :: cum_time

    ! ---

    first     = .true.
    cum_time  = 0.0_DP
    do i = 1, nr

       if (parents(i) == p) then

          if (shift == 0) then
             name = " " // names(i)
          else
             if (first) then
                name = repeat("    ", shift-1) // "   - " // names(i)
             else
                name = repeat("    ", shift-1) // "   - " // names(i)
             endif
          endif

#ifdef _OPENMP
          write (un, '(A60,2X,F12.3,2X,I8)') &
               name, wtimes(i), calls(i)

          cum_time  = cum_time + wtimes(i)
#else
          write (un, '(A60,2X,F12.3,2X,I8)') &
               name, times(i), calls(i)

          cum_time  = cum_time + times(i)
#endif

          call timer_print_for_parent(un, i, shift+1)

          first = .false.
       endif

    enddo

    if (p > 0 .and. .not. first) then
       name = repeat("    ", shift-1) // "   - " // "** remainder **"
#ifdef _OPENMP
       write (un, '(A60,2X,F12.3)') &
            name, wtimes(p) - cum_time
#else
       write (un, '(A60,2X,F12.3)') &
            name, times(p) - cum_time
#endif
    endif

  endsubroutine timer_print_for_parent


  !>
  !! Print timings to screen
  !<
  subroutine timer_print(un)
    implicit none

    integer, intent(in), optional  :: un

    ! ---

    integer  :: l_un

    ! ---

    l_un  = -1
    if (present(un)) then
       l_un  = un
    endif

    write (l_un, '(A)')  "====> TIMINGS <===="    
    write (l_un, '(60X,2X,A12,2X,A8)')  "time[s]", "calls"
    write (l_un, '(60X,2X,A12,2X,A8)')  "-------", "-----"

    call timer_print_for_parent(l_un, 0, 0)
    
  endsubroutine timer_print


  !>
  !! Print timings to log file
  !<
  subroutine timer_print_to_log() bind(C)
    implicit none

    if (mpi_id() == ROOT) then
      call timer_print(ilog)
    endif

  endsubroutine timer_print_to_log


  !>
  !! Start a timer with a given name. Accepts a zero terminated string for 
  !! the name
  !<
  subroutine c_timer_start(name, error) bind(C, name="timer_start")
    use, intrinsic :: iso_c_binding

    type(C_PTR),              value         :: name
#ifdef NO_BIND_C_OPTIONAL
    type(C_PTR),              value         :: error
#else
    integer(C_INT), optional, intent(inout) :: error
#endif

    ! ---

#ifdef NO_BIND_C_OPTIONAL
    integer(C_INT), pointer :: error_fptr
    call c_f_pointer(error, error_fptr)
    call timer_start(a2s(c_f_string(name)), error_fptr)
#else
    call timer_start(a2s(c_f_string(name)), error)
#endif

  endsubroutine c_timer_start


  !>
  !! Stop a timer with a given name. Accepts a zero terminated string for the
  !! name.
  !<
  subroutine c_timer_stop(name, error) bind(C, name="timer_stop")
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),              value         :: name
#ifdef NO_BIND_C_OPTIONAL
    type(C_PTR),              value         :: error
#else
    integer(C_INT), optional, intent(inout) :: error
#endif

    ! ---

#ifdef NO_BIND_C_OPTIONAL
    integer(C_INT), pointer :: error_fptr
    call c_f_pointer(error, error_fptr)
    call timer_stop(a2s(c_f_string(name)), error_fptr)
#else
    call timer_stop(a2s(c_f_string(name)), error)
#endif

  endsubroutine c_timer_stop

endmodule timer
