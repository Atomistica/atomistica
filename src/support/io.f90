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
!! Enhanced io feature, i.e., management of unit number
!<

#include "macros.inc"

module io
  use error_module
  use system_module

  implicit none

  private

  integer             :: io_i
  integer, parameter  :: n_possible_units = 1000
  integer, parameter  :: possible_units(n_possible_units) = (/ (io_i, io_i = 100, 100+n_possible_units-1) /)

  integer             :: open_units(n_possible_units) = -1

  public :: F_READ, F_WRITE
  integer, parameter  :: F_READ   = 1000
  integer, parameter  :: F_WRITE  = 1001

  public :: dump
  interface dump
     module procedure dump_real2d, dump_complex2d
  endinterface

  public :: read_ascii
  interface read_ascii
    module procedure reada_real_dim1, reada_int_dim1
  end interface read_ascii

  public :: get_unit, fopen, fclose, read_line

contains

  function get_unit()
    implicit none

    integer  :: get_unit

    ! ---

    integer  :: un, i

    ! ---

    un = -1
    do i = 1, n_possible_units
       if (.not. any(open_units == possible_units(i))) then
          un = possible_units(i)
       endif
    enddo

    if (un == -1) then
       stop "No unassigned unit found."
    endif

    i = 1
    do while (open_units(i) /= -1 .and. i < n_possible_units)
       i = i+1
    enddo
    open_units(i) = un

    get_unit = un

  endfunction get_unit


  !>
  !! Open a file and return the unit number
  !<
  function fopen(fn, mode, error)
    implicit none

    character*(*), intent(in)         :: fn
    integer, intent(in), optional     :: mode
    integer, intent(inout), optional  :: error

    integer                           :: fopen

    ! ---

    integer  :: un, i

    ! ---

    un = get_unit()

    if (present(mode)) then
       if (mode == F_WRITE) then
          open(un, file=fn, iostat=i, action="write")
       else if (mode == F_READ) then
          open(un, file=fn, iostat=i, action="read")
       else
          RAISE_ERROR("Internal error: Wrong mode provided.", error)
       endif
    else
       open(un, file=fn, iostat=i, action="read")
    endif

    if (i == 0) then
       fopen = un
    else
       RAISE_ERROR("Error opening file '" // trim(fn) // "'.", error)
    endif

  endfunction fopen


  !>
  !! Close a file and mark that unit as unused
  !<
  subroutine fclose(un, error)
    implicit none

    integer, intent(in)               :: un
    integer, intent(inout), optional  :: error

    ! ---

    integer  :: i
    logical  :: closed

    ! ---

    closed = .false.
    do i = 1, n_possible_units
       if (open_units(i) == un) then
          open_units(i) = -1
          closed = .true.
       endif
    enddo

    if (.not. closed) then
       RAISE_ERROR("Wrong unit number.", error)
    endif

    close(un)

  endsubroutine fclose


  !>
  !! Dump a real matrix
  !!
  !! Write matrix \param r to file \param fn
  !<
  subroutine dump_real2d(r, fn)
    implicit none

    real(DP), intent(in)      :: r(:, :)   !< matrix
    character(*), intent(in)  :: fn        !< file name

    ! ---

    character(13)  :: fmt
    integer        :: un, i

    ! ---

    write (fmt, '(A,I4.4,A)')  "(", size(r, 2), "ES20.10)"

    un = fopen(fn, F_WRITE)
    write (un, fmt)  (/ ( r(i, :), i = lbound(r, 1), ubound(r, 1) ) /)
    call fclose(un)

  endsubroutine dump_real2d


  !>
  !! Dump a complex matrix
  !!
  !! Write matrix \param r to file \param fn
  !<
  subroutine dump_complex2d(c, fn)
    implicit none

    complex(DP), intent(in)   :: c(:, :)   !< matrix
    character(*), intent(in)  :: fn        !< file name

    ! ---

    character(13)  :: fmt
    integer        :: un, i

    ! ---

    write (fmt, '(A,I4.4,A)')  "(", 2*size(c, 2), "ES20.10)"

    un = fopen(fn, F_WRITE)
    write (un, fmt)  (/ ( c(i, :), i = lbound(c, 1), ubound(c, 1) ) /)
    call fclose(un)

  endsubroutine dump_complex2d


  !>
  !! Read a line from the file
  !!
  !! Read a line from file \param f
  !<
  function read_line(f)
    implicit none

    integer, intent(in) :: f !< file unit
    character(1024)     :: read_line

    ! ---

    read (f, '(A)')  read_line

  endfunction read_line


  !>
  !! Read scalar and array data from ascii files. These
  !! interfaces are not yet heavily overloaded to cater for all intrinsic and
  !! most derived types.
  !<
  subroutine reada_real_dim1(un, da, status, error)
    integer,           intent(in)  :: un
    real(DP),          intent(out) :: da(:)
    integer, optional, intent(out) :: status
    integer, optional, intent(out) :: error

    ! ---

    integer :: my_status

    ! ---

    INIT_ERROR(error)

    if (present(status)) then
       read (un,fmt=*,iostat=status) da
    else
       read (un,fmt=*,iostat=my_status) da
       if (my_status < 0) then
          RAISE_ERROR("End of file.", error)
       endif
       if (my_status > 0) then
          RAISE_ERROR("Error reading.", error)
       endif
    endif
  endsubroutine reada_real_dim1


  !>
  !! Read scalar and array data from ascii files. These
  !! interfaces are not yet heavily overloaded to cater for all intrinsic and
  !! most derived types.
  !<
  subroutine reada_int_dim1(un, ia, status, error)
    integer,           intent(in)  :: un
    integer,           intent(out) :: ia(:)
    integer, optional, intent(out) :: status
    integer, optional, intent(out) :: error

    ! ---

    integer :: my_status

    ! ---

    INIT_ERROR(error)

    if (present(status)) then
       read (un,fmt=*,iostat=status) ia
    else
       read (un,fmt=*,iostat=my_status) ia
       if (my_status < 0) then
          RAISE_ERROR("End of file.", error)
       endif
       if (my_status > 0) then
          RAISE_ERROR("Error reading.", error)
       endif
    endif
  endsubroutine reada_int_dim1

endmodule io
