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

!>
!! Provide tools for C/Fortran interoperability
!<
module c_f
  use, intrinsic :: iso_c_binding

  private

  character(kind=C_CHAR), save, target :: dummy_string(6) = "(null)"
  character(kind=C_CHAR), save, target :: one_string(2) = "?"

  public :: c_f_string

contains

  !>
  !! Convert a null-terminated C string into a Fortran character array pointer
  !<
  function c_f_string(cptr) result(fptr)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),            intent(in) :: cptr ! The C address
    character(kind=C_CHAR), pointer    :: fptr(:)

    ! ---
      
    interface ! strlen is a standard C function from <string.h>
       ! int strlen(char *string)
       function strlen(string) result(len) bind(C, name="strlen")
         use, intrinsic :: iso_c_binding
         type(C_PTR), value :: string ! A C pointer
       endfunction strlen
    endinterface

    ! ---

    integer :: i
    character(kind=C_CHAR), pointer :: tmp(:)

    ! ---
      
    if (c_associated(cptr)) then
       i = strlen(cptr)
       if (i > 1) then
          call c_f_pointer(fptr=fptr, cptr=cptr, shape=[i])
       else
          ! Somehow this cannot handle len 1 strings
          call c_f_pointer(fptr=tmp, cptr=cptr, shape=[2])
          one_string(1) = tmp(1)
          one_string(2) = ' '
          fptr => one_string
       endif
    else
       ! To avoid segfaults, associate FPTR with a dummy target:
       fptr => dummy_string
    endif
            
  endfunction c_f_string

endmodule c_f
