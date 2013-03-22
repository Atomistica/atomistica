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

#include "macros.inc"
#include "filter.inc"

!>
!! Filter a type of atom
!!
!! Filter a type of atom
!<
module filter
  use supplib

  use particles

  implicit none

  private

  integer, parameter  :: MAX_EL_STR = 80

  public :: MAX_EL_STR
  public :: filter_from_string, filter_count, filter_pack, filter_unpack, filter_prlog

contains

  !>
  !! Convert a string to a filter
  !!
  !! Convert a string to a filter
  !<
  function filter_from_string(s, p, ierror)
    implicit none

    character*(*),     intent(in)    :: s
    type(particles_t), intent(in)    :: p
    integer, optional, intent(inout) :: ierror
    integer                          :: filter_from_string

    ! --

    integer         :: f, i, j, Z
    character(2)    :: sym
    character(200)  :: t

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

          Z    = atomic_number(sym)
          if (Z > 0) then
             do j = 1, p%nel
                if (p%el2Z(j) == Z) then
                   f = f + ishft(1, j)
                endif
             enddo
          else
             RAISE_ERROR("Unknown element '" // trim(sym) // "'.", ierror)
          endif
       enddo
       sym = t

       Z = atomic_number(sym)
       if (Z > 0) then
          do j = 1, p%nel
             if (p%el2Z(j) == Z) then
                f = f + ishft(1, j)
             endif
          enddo
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

integer, intent(in) :: f
    type(particles_t), intent(in) :: p
    integer, optional, intent(in) :: indent

    ! ---

    integer :: i
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
