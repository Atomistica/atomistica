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
!   shared:directory
! @endmeta

!<
!! Coulomb dispatch module.
!!
!! Coulomb dispatch module, the Python wrapper.
!!
!! This module contains a single Coulomb class which manages the individual Coulomb solver.
!! Since Fortran 90 does not support inheritance this is done manually, within this module.
!!
!! Additionally, the coulomb_t class manages conversion between different systems of units.
!!
!! Important: This is also the reference interface for all Coulomb modules.
!!
!! A typical use case would be:
!!
!!   type(particles_t)      :: p
!!   real(DP), allocatable  :: q(:)
!!   type(neighbors_t)      :: nl
!!
!!   type(coulomb_t)        :: coul
!!
!!   allocate(coul%direct_coulomb)
!!   call init(coul%direct_coulomb)   ! DirectCoulomb init takes no parameters
!!
!!   ... some code ...
!!
!!   call del(coul)
!!
!! Note on units:
!!   In eV/A units 1/epsilon_0 = 4 pi Hartree Bohr
!! 
!>

#include "macros.inc"

#include "have.inc"

module coulomb
  use, intrinsic :: iso_c_binding

  use supplib

  use particles
  use neighbors

#include "coulomb.inc"

  implicit none

  private

  public :: C_PTR
  public :: coulomb_set_Hubbard_U, coulomb_potential

  interface
    subroutine py_coulomb_set_Hubbard_U(this_cptr, p_cptr, U, ierror) bind(C)
      use, intrinsic :: iso_c_binding

      implicit none

      type(C_PTR),    value       :: this_cptr
      type(C_PTR),    value       :: p_cptr
      real(C_DOUBLE), intent(in)  :: U(*)
      integer(C_INT), intent(out) :: ierror
    endsubroutine py_coulomb_set_Hubbard_U

    subroutine py_coulomb_potential(this_cptr, p_cptr, nl_cptr, q, phi, &
         ierror) bind(C)
      use, intrinsic :: iso_c_binding

      implicit none

      type(C_PTR),    value         :: this_cptr
      type(C_PTR),    value         :: p_cptr
      type(C_PTR),    value         :: nl_cptr
      real(C_DOUBLE), intent(in)    :: q(*)
      real(C_DOUBLE), intent(inout) :: phi(*)
      integer(C_INT), intent(out)   :: ierror

    endsubroutine py_coulomb_potential
  endinterface

contains

  subroutine coulomb_set_Hubbard_U(this_cptr, p, U, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)  :: this_cptr
    type(particles_t), target      :: p
    real(DP),          intent(in)  :: U(*)
    integer, optional, intent(out) :: ierror

    ! ---
    
    integer :: ierror_loc
    
    ! ---

    if (present(ierror)) then
       INIT_ERROR(ierror)
       call py_coulomb_set_Hubbard_U(this_cptr, c_loc(p), U, ierror)
       PASS_ERROR(ierror)
    else
       ierror_loc = ERROR_NONE
       call py_coulomb_set_Hubbard_U(this_cptr, c_loc(p), U, ierror_loc)
       HANDLE_ERROR(ierror_loc)
    endif

  endsubroutine coulomb_set_Hubbard_U

  subroutine coulomb_potential(this_cptr, p, nl, q, phi, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), target        :: p
    type(neighbors_t), target        :: nl
    real(DP),          intent(in)    :: q(p%maxnatloc)
    real(DP),          intent(inout) :: phi(p%maxnatloc)
    integer, optional, intent(inout) :: ierror

    ! ---

    integer :: ierror_loc
    
    ! ---

    if (present(ierror)) then
       INIT_ERROR(ierror)
       call py_coulomb_potential(this_cptr, c_loc(p), c_loc(nl), q, phi, &
          ierror)
       PASS_ERROR(ierror)
    else
       ierror_loc = ERROR_NONE
       call py_coulomb_potential(this_cptr, c_loc(p), c_loc(nl), q, phi, &
          ierror_loc)
       HANDLE_ERROR(ierror_loc)
    endif
    
  endsubroutine coulomb_potential

endmodule coulomb
