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
!! Dictionary of pointers
!<

module ptrdict
  use, intrinsic :: iso_c_binding

  implicit none

  !
  ! Interface to the C-routines
  !

  interface
     function ptrdict_register_section(this, name, description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value  :: this
       character(kind=c_char, len=1)         :: name(*)
       character(kind=c_char, len=1)         :: description(*)

       type(c_ptr)                           :: ptrdict_register_section
     endfunction ptrdict_register_section


     function ptrdict_register_module(this, enabled, name, description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: enabled
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)

       type(c_ptr)                            :: ptrdict_register_module
     endfunction ptrdict_register_module


     subroutine ptrdict_register_integer_property(this, ptr, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value      :: this
       type(c_ptr),                   value      :: ptr
       character(kind=c_char, len=1)             :: name(*)
       character(kind=c_char, len=1)             :: description(*)
     endsubroutine ptrdict_register_integer_property


     subroutine ptrdict_register_real_property(this, ptr, name, &
       description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value      :: this
       type(c_ptr),                   value      :: ptr
       character(kind=c_char, len=1)             :: name(*)
       character(kind=c_char, len=1)             :: description(*)
     endsubroutine ptrdict_register_real_property


     subroutine ptrdict_register_boolean_property(this, ptr, name, &
       description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value      :: this
       type(c_ptr),                   value      :: ptr
       character(kind=c_char, len=1)             :: name(*)
       character(kind=c_char, len=1)             :: description(*)
     endsubroutine ptrdict_register_boolean_property


     subroutine ptrdict_register_string_property(this, ptr, maxlen, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value  :: this
       type(c_ptr),                   value  :: ptr
       integer(c_int),                value  :: maxlen
       character(kind=c_char, len=1)         :: name(*)
       character(kind=c_char, len=1)         :: description(*)
     endsubroutine ptrdict_register_string_property


     subroutine ptrdict_register_point_property(this, ptr, name, &
       description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value  :: this
       type(c_ptr),                   value  :: ptr
       character(kind=c_char, len=1)         :: name(*)
       character(kind=c_char, len=1)         :: description(*)
     endsubroutine ptrdict_register_point_property

     subroutine ptrdict_register_intpoint_property(this, ptr, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value  :: this
       type(c_ptr),                   value  :: ptr
       character(kind=c_char, len=1)         :: name(*)
       character(kind=c_char, len=1)         :: description(*)
     endsubroutine ptrdict_register_intpoint_property


     subroutine ptrdict_register_enum_property(this, ptr, nchoices, lenchoice, &
          choices, name, description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: nchoices
       integer(c_int),                value   :: lenchoice
       character(kind=c_char, len=1)          :: choices(lenchoice, nchoices)
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_enum_property


     subroutine ptrdict_register_list_property(this, ptr, maxlen, len, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: maxlen
       type(c_ptr),                   value   :: len
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_list_property


     subroutine ptrdict_register_integer_list_property(this, ptr, maxlen, len, &
       name, description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: maxlen
       type(c_ptr),                   value   :: len
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_integer_list_property
     

     subroutine ptrdict_register_string_list_property(this, ptr, strlen, &
          maxlen, len, name, description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       integer(c_int),                value   :: strlen
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: maxlen
       type(c_ptr),                   value   :: len
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_string_list_property


     subroutine ptrdict_register_array1d_property(this, ptr, nx, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: nx
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_array1d_property


     subroutine ptrdict_register_array2d_property(this, ptr, nx, ny, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: nx
       integer(c_int),                value   :: ny
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_array2d_property


     subroutine ptrdict_register_array3d_property(this, ptr, nx, ny, nz, name, &
          description) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value   :: this
       type(c_ptr),                   value   :: ptr
       integer(c_int),                value   :: nx
       integer(c_int),                value   :: ny
       integer(c_int),                value   :: nz
       character(kind=c_char, len=1)          :: name(*)
       character(kind=c_char, len=1)          :: description(*)
     endsubroutine ptrdict_register_array3d_property


     subroutine ptrdict_cleanup(root) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr), value  :: root
     endsubroutine ptrdict_cleanup

     
     subroutine ptrdict_read(root, fn) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value  :: root
       character(kind=c_char, len=1)         :: fn(*)
     endsubroutine ptrdict_read

     
     subroutine ptrdict_write(root, fn) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       type(c_ptr),                   value  :: root
       character(kind=c_char, len=1)         :: fn(*)
     endsubroutine ptrdict_write
  endinterface

  !
  ! Fortran-90 wrapper
  !

  type ptrdict_t
     type(c_ptr) :: ptrdict
  endtype ptrdict_t

  interface register_section
     module procedure fptrdict_register_section
  endinterface

  interface register
     module procedure fptrdict_register_integer_property
     module procedure fptrdict_register_real_property
     module procedure fptrdict_register_boolean_property
     module procedure fptrdict_register_string_property
     module procedure fptrdict_register_array1d_property
     module procedure fptrdict_register_array2d_property
     module procedure fptrdict_register_array3d_property
  endinterface

  interface cleanup
     module procedure fptrdict_cleanup
  endinterface

  interface read
     module procedure fptrdict_read
  endinterface

  interface write
     module procedure fptrdict_write
  endinterface

contains

  function fptrdict_register_section(this, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    type(ptrdict_t)                    :: fptrdict_register_section

    ! ---

    if (present(description)) then
       fptrdict_register_section%ptrdict = &
          ptrdict_register_section(this%ptrdict, CSTR(name), CSTR(description))
    else
       fptrdict_register_section%ptrdict = &
          ptrdict_register_section(this%ptrdict, CSTR(name), CSTR("N/A"))
    endif

  endfunction fptrdict_register_section


  subroutine fptrdict_register_integer_property(this, ptr, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    integer(C_INT),         target     :: ptr
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_integer_property(this%ptrdict, c_loc(ptr), &
                                              CSTR(name), CSTR(description))
    else
       call ptrdict_register_integer_property(this%ptrdict, c_loc(ptr), &
                                              CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_integer_property


  subroutine fptrdict_register_real_property(this, ptr, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_real_property(this%ptrdict, c_loc(ptr), &
                                           CSTR(name), CSTR(description))
    else
       call ptrdict_register_real_property(this%ptrdict, c_loc(ptr), &
                                           CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_real_property


  subroutine fptrdict_register_boolean_property(this, ptr, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    logical(C_BOOL),        target     :: ptr
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_boolean_property(this%ptrdict, c_loc(ptr), &
                                              CSTR(name), CSTR(description))
    else
       call ptrdict_register_boolean_property(this%ptrdict, c_loc(ptr), &
                                              CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_boolean_property


  subroutine fptrdict_register_string_property(this, ptr, maxlen, name, &
                                               description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    character(*),           target     :: ptr
    integer(c_int),         intent(in) :: maxlen
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_string_property(this%ptrdict, c_loc(ptr(1:1)), maxlen, &
                                             CSTR(name), CSTR(description))
    else
       call ptrdict_register_string_property(this%ptrdict, c_loc(ptr(1:1)), maxlen, &
                                             CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_string_property


  subroutine fptrdict_register_array1d_property_s(this, ptr, nx, name, &
                                                  description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr(nx)
    integer,                intent(in) :: nx 
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_array1d_property(this%ptrdict, c_loc(ptr), nx, &
                                              CSTR(name), CSTR(description))
    else
       call ptrdict_register_array1d_property(this%ptrdict, c_loc(ptr), nx, &
                                              CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_array1d_property_s


  subroutine fptrdict_register_array1d_property(this, ptr, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr(:)
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    call fptrdict_register_array1d_property_s(this, ptr, size(ptr, 1), name, &
                                              description)

  endsubroutine fptrdict_register_array1d_property


  subroutine fptrdict_register_array2d_property_s(this, ptr, nx, ny, name, &
                                                  description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr(nx, ny)
    integer,                intent(in) :: nx, ny
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_array2d_property(this%ptrdict, c_loc(ptr), &
                                              nx, ny, &
                                              CSTR(name), CSTR(description))
    else
       call ptrdict_register_array2d_property(this%ptrdict, c_loc(ptr), &
                                              nx, ny, &
                                              CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_array2d_property_s


  subroutine fptrdict_register_array2d_property(this, ptr, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr(:, :)
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    call fptrdict_register_array2d_property_s(this, ptr, &
                                              size(ptr, 1), size(ptr, 2), &
                                              CSTR(name), CSTR(description))

  endsubroutine fptrdict_register_array2d_property


  subroutine fptrdict_register_array3d_property_s(this, ptr, nx, ny, nz, &
                                                  name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr(nx, ny, nz)
    integer,                intent(in) :: nx, ny, nz
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    if (present(description)) then
       call ptrdict_register_array3d_property(this%ptrdict, c_loc(ptr), &
                                              size(ptr, 1), size(ptr, 2), &
                                              size(ptr, 3), &
                                              CSTR(name), CSTR(description))
    else
       call ptrdict_register_array3d_property(this%ptrdict, c_loc(ptr), &
                                              size(ptr, 1), size(ptr, 2), &
                                              size(ptr, 3), &
                                              CSTR(name), CSTR("N/A"))
    endif

  endsubroutine fptrdict_register_array3d_property_s


  subroutine fptrdict_register_array3d_property(this, ptr, name, description)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t),        intent(in) :: this
    real(C_DOUBLE),         target     :: ptr(:, :, :)
    character(*),           intent(in) :: name
    character(*), optional, intent(in) :: description

    ! ---

    call fptrdict_register_array3d_property_s(this, ptr, &
                                              size(ptr, 1), size(ptr, 2), &
                                              size(ptr, 3), &
                                              CSTR(name), CSTR(description))

  endsubroutine fptrdict_register_array3d_property


  subroutine fptrdict_cleanup(root)
    implicit none

    type(ptrdict_t), intent(in) :: root

    ! ---

    call ptrdict_cleanup(root%ptrdict)

  endsubroutine fptrdict_cleanup

     
  subroutine fptrdict_read(root, fn)
    implicit none

    type(ptrdict_t), intent(in) :: root
    character(*),    intent(in) :: fn

    ! ---

    call ptrdict_read(root%ptrdict, fn)

  endsubroutine fptrdict_read

     
  subroutine fptrdict_write(root, fn)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ptrdict_t), intent(in) :: root
    character(*),    intent(in) :: fn

    ! ---

    call ptrdict_write(root%ptrdict, fn)

  endsubroutine fptrdict_write
  
endmodule ptrdict
