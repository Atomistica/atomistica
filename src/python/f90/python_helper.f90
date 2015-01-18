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

!>
!! Python helper module
!!
!! Provides access to Fortran type objects from C
!<
module python_helper
  use supplib

  use particles
  use neighbors

  implicit none

  integer, parameter  :: MAX_STR_LEN = 1024

contains

  !>
  !! Convert zero terminated string to Fortran string
  !<
  function z2s(z) result(s)
    use, intrinsic :: iso_c_binding

    implicit none

    character(kind=c_char, len=1), intent(in)  :: z(*)
    character(MAX_STR_LEN)                     :: s

    ! ---

    integer  :: i

    ! ---

    s = ""

    i = 1
    do while (z(i) /= C_NULL_CHAR)
       s(i:i) = z(i)
       i = i+1
    enddo

  endfunction z2s


  !
  ! Particles stuff
  !


  !>
  !! Check if a field with this name already exists
  !<
  function f_data_exists(this_cptr, name, data_type) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none
 
    type(c_ptr),                   value        :: this_cptr
    character(kind=c_char, len=1)               :: name(*)
    integer(c_int),                intent(out)  :: data_type

    logical(c_bool)                             :: f_data_exists

    ! ---

    type(data_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    f_data_exists = exists(this, z2s(name), data_type)

  endfunction f_data_exists


  !>
  !! Return length of array stored in a data object
  !<
  function data_get_len(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    integer(c_int)      :: data_get_len

    ! ---

    type(data_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    data_get_len  = this%len

  endfunction data_get_len


  !>
  !! Return a pointer to the field data
  !<
  subroutine real_ptr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    real(c_double), pointer  :: ptr(:)

    ! ---

    call c_f_pointer(this_cptr, this)
    call ptr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1))

  endsubroutine real_ptr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine integer_ptr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    integer(c_int), pointer  :: ptr(:)

    ! ---

    call c_f_pointer(this_cptr, this)
    call ptr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1))

  endsubroutine integer_ptr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine realx_ptr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    real(c_double), pointer  :: ptr(:, :)

    ! ---

    call c_f_pointer(this_cptr, this)
    call ptr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1, 1))

  endsubroutine realx_ptr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine realxxx_ptr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    real(c_double), pointer  :: ptr(:, :, :)

    ! ---

    call c_f_pointer(this_cptr, this)
    call ptr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1, 1, 1))

  endsubroutine realxxx_ptr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine real_attr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    real(c_double), pointer  :: ptr

    ! ---

    call c_f_pointer(this_cptr, this)
    call attr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr)

  endsubroutine real_attr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine real3_attr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    real(c_double), pointer  :: ptr(:)

    ! ---

    call c_f_pointer(this_cptr, this)
    call attr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1))

  endsubroutine real3_attr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine real3x3_attr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    real(c_double), pointer  :: ptr(:, :)

    ! ---

    call c_f_pointer(this_cptr, this)
    call attr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1, 1))

  endsubroutine real3x3_attr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine integer_attr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    integer(c_int), pointer  :: ptr

    ! ---

    call c_f_pointer(this_cptr, this)
    call attr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr)

  endsubroutine integer_attr_by_name


  !>
  !! Return a pointer to the field data
  !<
  subroutine integer3_attr_by_name(this_cptr, name, ptr_cptr, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),                  value          :: this_cptr
    character(kind=c_char, len=1)                :: name(*)
    type(c_ptr),                  intent(out)    :: ptr_cptr
    integer(c_int),               intent(inout)  :: ierror

    ! ---

    type(data_t),   pointer  :: this
    integer(c_int), pointer  :: ptr(:)

    ! ---

    call c_f_pointer(this_cptr, this)
    call attr_by_name(this, z2s(name), ptr, ierror)
    ptr_cptr = c_loc(ptr(1))

  endsubroutine integer3_attr_by_name


  !>
  !! Open log file
  !<
  subroutine f_logging_start(fn) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    character(kind=c_char, len=1)  :: fn(*)

    ! ---

    call logging_start(z2s(fn))

  endsubroutine f_logging_start


  !>
  !! Return error string
  !<
  subroutine get_full_error_string(str) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    character(kind=c_char, len=1), intent(inout)  :: str(*)

    ! ---

    integer          :: i, l
    character(1000)  :: errstr

    ! ---
    
    errstr = get_error_string_and_clear()
    l = len_trim(errstr)
    do i = 1, l
       str(i) = errstr(i:i)
    enddo
    str(l+1) = C_NULL_CHAR

  endsubroutine get_full_error_string

endmodule python_helper
