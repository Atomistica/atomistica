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

  public :: resize
  interface resize
     module procedure resize_int
     module procedure resize_real, resize_real2
     module procedure resize_complex
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


  !>
  !! Resize integer array while keeping contents
  !!
  !! Resize array v to length n, keeping the contents. (If n is smaller than
  !! the current size, obviously the elements exceeding the new size are
  !! lost.)
  !!
  !! If the new size is zero, the array will be deallocated. Likewise, if the
  !! old size is zero, the array will be allocated.
  !<
  subroutine resize_int(v, n)
    implicit none

    integer, dimension(:), allocatable, intent(inout) :: v  !< Array to be resized
    integer                                           :: n  !< New length

    ! ---

    integer                            :: i     ! loops
    integer, dimension(:), allocatable :: newv  ! new array
    integer                            :: on    ! max i to copy

    ! ---

    if(n==0) then
       if(allocated(v)) then
          deallocate(v)
       end if
    else
       if (allocated(v) .and. size(v) == n)  return

       allocate(newv(n))
       newv = 0

       if(allocated(v)) then
          on = min(size(v), n)
          do i = 1, on
             newv(i) = v(i)
          enddo

          deallocate(v)
       endif

       allocate(v(n))
       v = newv
       deallocate(newv)
    end if

  end subroutine resize_int


  !>
  !! Resize real array while keeping contents
  !!
  !! Resize array v to length n, keeping the contents. (If n is smaller than
  !! the current size, obviously the elements exceeding the new size are
  !! lost.)
  !!
  !! If the new size is zero, the array will be deallocated. Likewise, if the
  !! old size is zero, the array will be allocated.
  !<
  subroutine resize_real(v, n)
    implicit none

    real(DP), dimension(:), allocatable, intent(inout) :: v  !< Array to be resized
    integer                                            :: n  !< New length

    ! ---

    integer                             :: i     ! loops
    real(DP), dimension(:), allocatable :: newv  ! new array
    integer                             :: on    ! max i to copy

    ! ---

    if(n==0) then
       if(allocated(v)) then
          deallocate(v)
       end if
    else
       if (allocated(v) .and. size(v) == n)  return

       allocate(newv(n))
       newv = 0

       if(allocated(v)) then
          on = min(size(v), n)
          do i = 1, on
             newv(i) = v(i)
          enddo

          deallocate(v)
       endif

       allocate(v(n))
       v = newv
       deallocate(newv)
    end if

  end subroutine resize_real


  !>
  !! Resize real array while keeping contents
  !!
  !! Resize 2d array v to length n1 x n2, keeping the contents.
  !! (If n is smaller than
  !! the current size, obviously the elements exceeding the new size are
  !! lost.)
  !!
  !! If the new size is zero, the array will be deallocated. Likewise, if the
  !! old size is zero, the array will be allocated.
  !<
  subroutine resize_real2(v, n1, n2)
    implicit none

    real(DP), dimension(:, :), allocatable, intent(inout) :: v   !< Array to be resized
    integer                                               :: n1  !< New length
    integer                                               :: n2  !< New length

    ! ---

    integer                                :: i, j        ! loops
    real(DP), dimension(:, :), allocatable :: newv        ! new array
    integer                                :: on1, on2    ! max i to copy

    ! ---

    if(n1 == 0 .or. n2 == 0) then
       if(allocated(v)) then
          deallocate(v)
       end if
    else
       if (allocated(v) .and. size(v,1) == n1 .and. size(v,2) == n2)  return

       allocate(newv(n1, n2))
       newv = 0

       if(allocated(v)) then
          on1 = min(size(v, 1), n1)
          on2 = min(size(v, 2), n2)
          do i = 1, on2
             do j = 1, on1
                newv(j, i) = v(j, i)
             enddo
          enddo

          deallocate(v)
       endif

       allocate(v(n1, n2))
       v = newv
       deallocate(newv)
    end if

  end subroutine resize_real2


  !>
  !! Resize real array while keeping contents
  !!
  !! Resize array v to length n, keeping the contents. (If n is smaller than
  !! the current size, obviously the elements exceeding the new size are
  !! lost.)
  !!
  !! If the new size is zero, the array will be deallocated. Likewise, if the
  !! old size is zero, the array will be allocated.
  !<
  subroutine resize_complex(v, n)
    implicit none

    complex(DP), dimension(:), allocatable, intent(inout) :: v  !< Array to be resized
    integer                                               :: n  !< New length

    ! ---

    integer                                :: i     ! loops
    complex(DP), dimension(:), allocatable :: newv  ! new array
    integer                                :: on    ! max i to copy

    ! ---

    if(n==0) then
       if(allocated(v)) then
          deallocate(v)
       end if
    else
       if (allocated(v) .and. size(v) == n)  return

       allocate(newv(n))
       newv = 0

       if(allocated(v)) then
          on = min(size(v), n)
          do i = 1, on
             newv(i) = v(i)
          enddo

          deallocate(v)
       endif

       allocate(v(n))
       v = newv
       deallocate(newv)
    end if

  end subroutine resize_complex

endmodule misc
