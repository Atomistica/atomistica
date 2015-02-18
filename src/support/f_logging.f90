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
!! Global logging capabilities
!<
module logging
  use, intrinsic :: iso_c_binding

  use system_module
  use c_f
  use io
  use mpi_context_module

  implicit none

  private

  integer, parameter  :: BYTES_PER_MB  = 1024*1024

  integer             :: ilog = -1
  real(DP)            :: total_memory = 0.0_DP


  interface log_memory_estimate
     module procedure log_memory_estimate_integer4
     module procedure log_memory_estimate_integer4_2
     module procedure log_memory_estimate_integer4_3
     module procedure log_memory_estimate_integer4_4
     module procedure log_memory_estimate_integer8
     module procedure log_memory_estimate_integer8_2
     module procedure log_memory_estimate_integer8_3
     module procedure log_memory_estimate_integer8_4
     module procedure log_memory_estimate_logical
     module procedure log_memory_estimate_logical2
     module procedure log_memory_estimate_logical3
     module procedure log_memory_estimate_logical4
     module procedure log_memory_estimate_real
     module procedure log_memory_estimate_real2
     module procedure log_memory_estimate_real3
     module procedure log_memory_estimate_real4
     module procedure log_memory_estimate_complex
     module procedure log_memory_estimate_complex2
     module procedure log_memory_estimate_complex3
     module procedure log_memory_estimate_complex4
  endinterface

  public :: logging_start, logging_stop, prscrlog, prlog
  public :: log_memory_start, log_memory_stop, log_memory_estimate
  public :: log_memory_general, ilog

contains

  !>
  !! Open log file
  !<
  subroutine logging_start(fn)
    implicit none

    character(*), intent(in)  :: fn

    ! ---

    call logging_stop
    if (mpi_id() == ROOT) then
       ilog = fopen(fn, mode=F_WRITE)
    endif

  endsubroutine logging_start


  !>
  !! Close log file
  !<
  subroutine logging_stop() bind(C)
    implicit none

    if (ilog /= -1) then
      call fclose(ilog)
    endif
    ilog  = -1

  endsubroutine logging_stop


  !>
  !! Record a log message to screen and file
  !<
  subroutine prscrlog(msg)
    implicit none

    character(*), intent(in), optional  :: msg

    ! ---

    if (present(msg)) then
       if (ilog /= -1) then
#if !defined(MDCORE_PYTHON) && !defined(LAMMPS)
       ! Do not print to screen if we're using the Python or LAMMPS module
         write (*, '(A)')  msg
#endif
         write (ilog, '(A)')  msg
       endif
    else
       if (ilog /= -1) then
#if !defined(MDCORE_PYTHON) && !defined(LAMMPS)
         write (*, *)
#endif
         write (ilog, *)
       endif
    endif

  endsubroutine prscrlog


  !>
  !! Record a log message to file only
  !<
  subroutine prlog(msg)
    implicit none

    character(*), intent(in), optional  :: msg

    ! ---

    if (ilog /= -1) then
      if (present(msg)) then
         if (msg(1:1) == "-") then
            write (ilog, '(A)')  trim(msg)
         else
            write (ilog, '(5X,A)')  trim(adjustl(msg))
         endif
      else
         write (ilog, *)
      endif
    endif

  endsubroutine prlog


  !>
  !! Start logging of memory estimates
  !<
  subroutine log_memory_start(name)
    implicit none

    character(*), intent(in)  :: name

    ! ---

    total_memory  = 0.0_DP
    
  endsubroutine log_memory_start


  !>
  !! Stop logging of memory estimates
  !<
  subroutine log_memory_stop(name)
    implicit none

    character(*), intent(in)  :: name

    ! ---

    call prlog("Memory estimate: " // total_memory // " MB")
    
  endsubroutine log_memory_stop


  !>
  !! Print a memory usage estimate to the log file
  !<
  subroutine log_memory_general(bytes, str)
    implicit none

    integer, intent(in)                 :: bytes
    character(*), intent(in), optional  :: str

    ! ---

    real(DP)  :: m

    ! ---

    m  = real(bytes, DP)/BYTES_PER_MB

!    if (m > 1.0) then
!       write (ilog, '(5X,A,A,A,F10.3,A)') &
!            "Memory usage of array ", trim(str), ": ", m, " MB"
!    endif

    total_memory  = total_memory + m

  endsubroutine log_memory_general


#define LOG_MEMORY1(data_type, name, elsize)  \
  subroutine name(arr, str) ; implicit none ; data_type, intent(in)  :: arr(:) ; character(*), intent(in), optional  :: str ; call log_memory_general(size(arr)*elsize, str) ; endsubroutine name

#define LOG_MEMORY2(data_type, name, elsize)  \
  subroutine name(arr, str) ; implicit none ; data_type, intent(in)  :: arr(:, :) ; character(*), intent(in), optional  :: str ; call log_memory_general(size(arr)*elsize, str) ; endsubroutine name

#define LOG_MEMORY3(data_type, name, elsize)  \
  subroutine name(arr, str) ; implicit none ; data_type, intent(in)  :: arr(:, :, :) ; character(*), intent(in), optional  :: str ; call log_memory_general(size(arr)*elsize, str) ; endsubroutine name

#define LOG_MEMORY4(data_type, name, elsize)  \
  subroutine name(arr, str) ; implicit none ; data_type, intent(in)  :: arr(:, :, :, :) ; character(*), intent(in), optional  :: str ; call log_memory_general(size(arr)*elsize, str) ; endsubroutine name


  LOG_MEMORY1(integer(4), log_memory_estimate_integer4, 8)
  LOG_MEMORY2(integer(4), log_memory_estimate_integer4_2, 8)
  LOG_MEMORY3(integer(4), log_memory_estimate_integer4_3, 8)
  LOG_MEMORY4(integer(4), log_memory_estimate_integer4_4, 8)

  LOG_MEMORY1(integer(8), log_memory_estimate_integer8, 8)
  LOG_MEMORY2(integer(8), log_memory_estimate_integer8_2, 8)
  LOG_MEMORY3(integer(8), log_memory_estimate_integer8_3, 8)
  LOG_MEMORY4(integer(8), log_memory_estimate_integer8_4, 8)

  LOG_MEMORY1(logical, log_memory_estimate_logical, 8)
  LOG_MEMORY2(logical, log_memory_estimate_logical2, 8)
  LOG_MEMORY3(logical, log_memory_estimate_logical3, 8)
  LOG_MEMORY4(logical, log_memory_estimate_logical4, 8)

  LOG_MEMORY1(real(DP), log_memory_estimate_real, DP)
  LOG_MEMORY2(real(DP), log_memory_estimate_real2, DP)
  LOG_MEMORY3(real(DP), log_memory_estimate_real3, DP)
  LOG_MEMORY4(real(DP), log_memory_estimate_real4, DP)

  LOG_MEMORY1(complex(DP), log_memory_estimate_complex, 2*DP)
  LOG_MEMORY2(complex(DP), log_memory_estimate_complex2, 2*DP)
  LOG_MEMORY3(complex(DP), log_memory_estimate_complex3, 2*DP)
  LOG_MEMORY4(complex(DP), log_memory_estimate_complex4, 2*DP)


  !>
  !! Record a log message to screen and file
  !<
  subroutine c_prscrlog(msg) bind(C)
    use, intrinsic :: iso_c_binding
    type(C_PTR), value :: msg
    call prscrlog(a2s(c_f_string(msg)))
  endsubroutine c_prscrlog


  !>
  !! Record a log message to file only
  !<
  subroutine c_prlog(msg) bind(C)
    use, intrinsic :: iso_c_binding
    type(C_PTR), value :: msg
    call prlog(a2s(c_f_string(msg)))
  endsubroutine c_prlog

endmodule
