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
  use libAtoms_module

  use particles

  implicit none

  private

  integer, parameter  :: MAX_EL_STR = 80

  public :: MAX_EL_STR
  public :: filter_from_string

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

endmodule filter
