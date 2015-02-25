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
!! Coulomb dispatch module.
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

  implicit none

  private

  ! Note: coulomb_t is hidden. Everything is passed as type(C_PTR) to hide the
  ! complexity of coulomb_t from the compiler. This speeds up compile times
  ! and avoids nasty compiler crashes. However, this invalidates Fortran
  ! interfaces since the compiler can't match a generic call to datatype.

  public :: C_PTR

  public :: coulomb_alloc, coulomb_free

  public :: coulomb_del, coulomb_bind_to, coulomb_set_Hubbard_U
  public :: coulomb_potential, coulomb_energy_and_forces

contains

  !>
  !! Allocator
  !!
  !! Allocate memory for new coulomb instance
  !<
  subroutine coulomb_alloc(this_cptr)
    implicit none

    type(C_PTR), intent(out) :: this_cptr

  endsubroutine coulomb_alloc


  !>
  !! Free memory
  !!
  !! Free memory occupied by a coulomb instance
  !<
  subroutine coulomb_free(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

  endsubroutine coulomb_free


  !>
  !! Destructor
  !!
  !! Delete the Coulomb dispatch object and all allocated objects driver
  !<
  subroutine coulomb_del(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

  endsubroutine coulomb_del


  !>
  !! Bind to a certain Particles and Neighbors object
  !!
  !! Bind to a certain Particles and Neighbors object
  !<
  subroutine coulomb_bind_to(this_cptr, p, nl, ierror)
    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(inout) :: p
    type(neighbors_t), intent(inout) :: nl
    integer, optional, intent(inout) :: ierror

    ! ---

  endsubroutine coulomb_bind_to


  !>
  !! Set Hubbard-Us for all the elements
  !!
  !! Set Hubbard-Us for all the elements
  !<
  subroutine coulomb_set_Hubbard_U(this_cptr, p, U, ierror)
    implicit none

    type(C_PTR),       intent(in)  :: this_cptr
    type(particles_t), intent(in)  :: p
    real(DP),          intent(in)  :: U(:)
    integer, optional, intent(out) :: ierror

    ! ---

  endsubroutine coulomb_set_Hubbard_U


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !!
  !! Calculate the electrostatic potential of every atom (for variable charge models). Note that \param phi
  !! will be overriden.
  !<
  subroutine coulomb_potential(this_cptr, p, nl, q, phi, ierror)
    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(in)    :: p
    type(neighbors_t), intent(inout) :: nl
    real(DP),          intent(in)    :: q(p%maxnatloc)
    real(DP),          intent(inout) :: phi(p%maxnatloc)
    integer, optional, intent(inout) :: ierror

    ! ---

  endsubroutine coulomb_potential


  !>
  !! Calculate the total energy and all forces
  !!
  !! Returns the total (Coulomb) energy, all forces and optionally the virial contribution.
  !! Note that only the diagonal of the virial is correct right now.
  !!
  !! This assumes that both, positions and charges, of the atoms have changed.
  !<
  subroutine coulomb_energy_and_forces(this_cptr, p, nl, q, epot, f, wpot, &
       ierror)
    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(particles_t),  intent(in)    :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(in)    :: q(p%maxnatloc)
    real(DP),           intent(inout) :: epot
    real(DP), optional, intent(inout) :: f(3, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot(3, 3)
    integer,  optional, intent(inout) :: ierror

    ! ---

  endsubroutine coulomb_energy_and_forces

endmodule coulomb

