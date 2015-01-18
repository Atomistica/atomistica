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
! @meta
!   shared
!   classtype:no_integration_t classname:NoIntegration interface:integrators
! @endmeta

!>
!! Switch off integration altogether (for testing purposes)
!<

#include "macros.inc"

module no_integration
  use libAtoms_module

  use particles

  implicit none

  private

  public :: no_integration_t
  type no_integration_t

     integer   :: enabled  = 0

  endtype no_integration_t

  public :: register
  interface register
    module procedure no_integration_register
  endinterface

contains

  subroutine no_integration_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(no_integration_t), target, intent(inout)  :: this
    type(c_ptr),                    intent(in)     :: cfg
    type(c_ptr),                    intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("NoIntegration"), &
         CSTR("No integration."))

  endsubroutine no_integration_register

endmodule no_integration
