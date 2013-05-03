!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================
!>
!! Misc stuff
!<

#include "macros.inc"

module misc
  use error_module
  use system_module

  use logging

  implicit none

  private

  public :: swap
  interface swap
     module procedure swap_real, swap_int, swap_logical, swap_sym
  endinterface

  public :: equal, uppercase, xyz2index, index2xyz

contains

  !>
  !! Convert x,y,z coordinates into a global (scalar,compact) index
  !<
  pure function xyz2index(x, y, z, n)
    implicit none

    integer, intent(in)  :: x, y, z
    integer, intent(in)  :: n(3)
    integer              :: xyz2index

    ! ---

    xyz2index = (x-1) + n(1)*((y-1) + n(2)*(z-1)) + 1

  endfunction xyz2index


  !>
  !! Extract x,y,z coordinates from a global (scalar,compact) index
  !<
  pure subroutine index2xyz(i, n, x, y, z)
    implicit none

    integer, intent(in)   :: i
    integer, intent(in)   :: n(3)
    integer, intent(out)  :: x, y, z

    ! ---

    x = mod( i - 1, n(1) )
    y = mod( (i - 1 - x)/n(1), n(2) )
    z = ( i - 1 - x - n(1)*y )/(n(1)*n(2))

    x = x+1
    y = y+1
    z = z+1
    
  endsubroutine index2xyz


  !>
  !! Convert string to upper case
  !<
  pure function uppercase(s)
    implicit none

    character(*), intent(in)  :: s

    character(len(s))         :: uppercase

    ! ---

    integer  :: i

    ! ---

    do i = 1, len(s)
       if (ichar(s(i:i)) >= ichar('a') .and. ichar(s(i:i)) <= ichar('z')) then
          uppercase(i:i)  = char(ichar(s(i:i)) + ichar('A')-ichar('a'))
       else
          uppercase(i:i)  = s(i:i)
       endif
    enddo

  endfunction uppercase


  !>
  !! Check if two stings are equal, independent of case
  !<
  pure function equal(s1, s2)
    implicit none

    character(*), intent(in)  :: s1
    character(*), intent(in)  :: s2

    logical                   :: equal

    ! ---

    equal  = uppercase(trim(s1)) == uppercase(trim(s2))

  endfunction equal


  !>
  !! Swap a real
  !<
  elemental subroutine swap_real(val1, val2)
    implicit none

    real(DP), intent(inout)  :: val1, val2
    real(DP)                 :: tmp

    tmp = val1
    val1 = val2
    val2 = tmp

  endsubroutine swap_real


  !>
  !! Swap an integer
  !<
  elemental subroutine swap_int(val1, val2)
    implicit none

    integer, intent(inout)  :: val1, val2
    integer                 :: tmp

    tmp = val1
    val1 = val2
    val2 = tmp

  endsubroutine swap_int


  !>
  !! Swap a logical
  !<
  elemental subroutine swap_logical(val1, val2)
    implicit none

    logical, intent(inout)  :: val1, val2
    logical                 :: tmp

    tmp = val1
    val1 = val2
    val2 = tmp

  endsubroutine swap_logical


  !>
  !! Swap element symbols
  !<
  elemental subroutine swap_sym(val1, val2)
    implicit none

    character*4, intent(inout)  :: val1, val2
    character*4                 :: tmp

    tmp = val1
    val1 = val2
    val2 = tmp

  endsubroutine swap_sym

endmodule misc
