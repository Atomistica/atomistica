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

!>
!! Python helper module
!!
!! Provides access to Fortran type objects from C
!<
module python_helper
  use libAtoms_module

  use logging

  use data
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
  function data_exists(this_cptr, name, data_type) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none
 
    type(c_ptr),                   value        :: this_cptr
    character(kind=c_char, len=1)               :: name(*)
    integer(c_int),                intent(out)  :: data_type

    logical(c_bool)                             :: data_exists

    ! ---

    type(data_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    data_exists = exists(this, z2s(name), data_type)

  endfunction data_exists


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
