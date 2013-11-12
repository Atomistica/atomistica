!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
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

endmodule ptrdict
