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

#include "macros.inc"
#include "filter.inc"

!>
!! Filter a type of atom
!!
!! Filter a type of atom
!<
module filter
  use supplib

  use logging

  use data
  use particles

  implicit none

  private

  integer, parameter  :: MAX_EL_STR = 80

  public :: MAX_EL_STR
  public :: filter_from_string, filter_count, filter_sum, filter_average
  public :: filter_mask, filter_pack, filter_unpack, filter_prlog

contains

  !>
  !! Convert a string to a filter
  !!
  !! Convert a string to a filter
  !<
  function filter_from_string(s, p, ierror)
    implicit none

    character*(*), intent(in)         :: s
    type(particles_t), intent(in)     :: p
    integer, intent(inout), optional  :: ierror
    integer                           :: filter_from_string

    ! --

    integer          :: f, i, Z
    character(2)     :: sym
    character(200)   :: t

    ! --

    if (p%nel > 16) then
       RAISE_ERROR("More than 16 elements. Please upgrade the filter module to deal with this.", ierror)
    endif

    f = 0

    if (trim(s) == "*") then

       do i = 1, p%nel
          f = f + ishft(1, i)
       enddo

    else

       t = s
       i = scan(t, ',')
       do while (i /= 0)
          sym  = t(1:i-1)
          t    = t(i+1:)
          i    = scan(t, ',')

          Z = atomic_number(sym)
          if (Z > 0) then
             if (p%Z2el(Z) > 0) then
                f = f + ishft(1, p%Z2el(Z))
!             else
!                RAISE_ERROR("Element '" // trim(sym) // "' known but not present in this simulation.", ierror)
             endif
          else
             RAISE_ERROR("Unknown element '" // trim(sym) // "'.", ierror)
          endif
       enddo
       sym = t

       Z = atomic_number(sym)
       if (Z > 0) then
          if (p%Z2el(Z) > 0) then
             f = f + ishft(1, p%Z2el(Z))
!          else
!             RAISE_ERROR("Element '" // trim(sym) // "' known but not present in this simulation.", ierror)
          endif
       else
          RAISE_ERROR("Unknown element '" // trim(sym) // "'.", ierror)
       endif

    endif

    filter_from_string = f

  endfunction filter_from_string


  !>
  !! Count how many atoms that match this filter we have locally
  !!
  !! Count how many atoms that match this filter we have locally
  !<
  function filter_count(f, p)
    implicit none

    integer, intent(in)            :: f
    type(particles_t), intent(in)  :: p
    integer                        :: filter_count

    ! ---

    integer  :: i, n

    ! ---

    n = 0
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          n = n + 1
       endif
    enddo

    filter_count = n

  endfunction filter_count


  !>
  !! Sum field x
  !!
  !! Sum field x
  !<
  function filter_sum(f, p, x)
    implicit none

    integer, intent(in)            :: f
    type(particles_t), intent(in)  :: p
    real(DP), intent(in)           :: x(p%natloc)
    integer                        :: filter_sum

    ! ---

    integer  :: i

    ! ---

    filter_sum = 0.0_DP
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          filter_sum = filter_sum + x(i)
       endif
    enddo

  endfunction filter_sum


  !>
  !! Sum field x
  !!
  !! Sum field x
  !<
  function filter_average(f, p, x, mpi)
    implicit none

    integer,                     intent(in)  :: f
    type(particles_t),           intent(in)  :: p
    real(DP),                    intent(in)  :: x(p%natloc)
    type(MPI_context), optional, intent(in)  :: mpi
    integer                                  :: filter_average

    ! ---

    integer  :: i, n

    ! ---

    filter_average = 0.0_DP
    n = 0
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          filter_average = filter_average + x(i)
          n = n + 1
       endif
    enddo

#ifdef _MP
    if (present(mpi)) then
       call sum_in_place(mpi, filter_average)
       call sum_in_place(mpi, n)
    endif
#endif

    filter_average = filter_average/n

  endfunction filter_average


  !>
  !! Sum field x
  !!
  !! Sum field x
  !<
  function filter_mask(f, p)
    implicit none

    integer, intent(in)            :: f
    type(particles_t), intent(in)  :: p
    logical                        :: filter_mask(1:p%maxnatloc)

    ! ---

    integer  :: i

    ! ---

    filter_mask = .false.
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          filter_mask(i) = .true.
       endif
    enddo

  endfunction filter_mask


  !>
  !! Pack property into an array
  !!
  !! Pack property into an array
  !<
  subroutine filter_pack(f, p, unpacked, packed)
    implicit none

    integer, intent(in)            :: f
    type(particles_t), intent(in)  :: p
    real(DP), intent(in)           :: unpacked(p%maxnatloc)
    real(DP), intent(inout)        :: packed(*)

    ! ---

    integer  :: i, n

    ! ---

    n = 0
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          n = n + 1
          packed(n)  = unpacked(i)
       endif
    enddo

  endsubroutine filter_pack


  !>
  !! Pack property into an array
  !!
  !! Pack property into an array
  !<
  subroutine filter_unpack(f, p, packed, unpacked)
    implicit none

    integer, intent(in)            :: f
    type(particles_t), intent(in)  :: p
    real(DP), intent(in)           :: packed(*)
    real(DP), intent(inout)        :: unpacked(p%maxnatloc)

    ! ---

    integer  :: i, n

    ! ---

    n = 0
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          n = n + 1
          unpacked(i)  = packed(n)
       endif
    enddo

  endsubroutine filter_unpack


  !>
  !! Dump filter information to log file
  !!
  !! Dump filter information to log file
  !<
  subroutine filter_prlog(f, p, indent)
    implicit none

    integer,           intent(in) :: f
    type(particles_t), intent(in) :: p
    integer, optional, intent(in) :: indent

    ! ---

    integer         :: i
    character(1024) :: s

    ! ---

    s = " "
    do i = 1, p%nel
       if (IS_EL2(f, i)) then
          s = trim(s)//trim(ElementName(p%el2Z(i)))//","
       endif
    enddo

    if (len_trim(s) > 0) then
       s = s(1:len_trim(s)-1)
    endif

    s = trim(s)//" ("
    do i = 1, p%nel
       if (IS_EL2(f, i)) then
          s = trim(s)//"X"
       else
          s = trim(s)//"_"
       endif
    enddo
    s = trim(s)//")"

    if (present(indent)) then
       do i = 1, indent
          s = " " // s
       enddo
    endif

    call prlog(s)

  endsubroutine filter_prlog

endmodule filter
