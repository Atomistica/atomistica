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

!<
!! Potentials dispatch module.
!!
!! This module contains a single Potentials class which manages the individual
!! interatomic potentials. Since Fortran 90 does not support inheritance this
!! is done manually, within this module.
!!
!! This is also the reference interface for all Potentials.
!!
!! A typical use case would be:
!!
!!   type(particles_t)  :: p
!!   type(neighbors_t)  :: nl
!!
!!   type(potentials_t)  :: pot
!!
!!   allocate(pot%rebo2(1))
!!   call init(pot%rebo2(1))   ! REBO2 init takes no parameters
!!
!!   ... some code ...
!!
!!   call energy_and_forces(pot, p, nl, epot, f, wpot)
!!
!!   ... some code ...
!!
!!   call del(pot)
!!
!>

#include "macros.inc"

#include "have.inc"

module potentials
  use, intrinsic :: iso_c_binding

  use supplib

  use particles
  use neighbors
  use dynamics

  use coulomb

  use {classname}

  implicit none

  private

    public :: potentials_t
    type potentials_t

     !
     ! general buffers, will be filled if allocated
     !

     ! Note: The first two are pointes such that they may be stored in the 
     ! Atoms dynamic data structure
     real(DP), pointer  :: epot_per_at(:) => NULL()
     real(DP), pointer  :: wpot_per_at(:, :, :) => NULL()
     real(DP), pointer  :: epot_per_bond(:) => NULL()
     real(DP), pointer  :: f_per_bond(:, :) => NULL()
     real(DP), pointer  :: wpot_per_bond(:, :, :) => NULL()

     !
     ! embedded atom potentials
     !

     type({classtype}), allocatable :: {classname}(:)

  endtype potentials_t


  ! Note: potentials_t is hidden. Everything is passed as type(C_PTR) to hide
  ! the complexity of potentials_t from the compiler. This speeds up compile
  ! times and avoids nasty compiler crashes. However, this invalidates Fortran
  ! interfaces since the compiler can't match a generic call to datatype.

  public :: potentials_alloc, potentials_free
  public :: potentials_register_data, potentials_init
  public :: potentials_del, potentials_set_Coulomb
  public :: potentials_bind_to, potentials_energy_and_forces

contains

  !>
  !! Allocator
  !!
  !! Allocate memory for new potentials instance
  !<
  subroutine potentials_alloc(this_cptr)
    implicit none

    type(C_PTR), intent(out) :: this_cptr

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    allocate(this)
    this_cptr = c_loc(this)

  endsubroutine potentials_alloc


  !>
  !! Free memory
  !!
  !! Free memory occupied by a potentials instance
  !<
  subroutine potentials_free(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    deallocate(this)

  endsubroutine potentials_free


  !>
  !! Register per-atom fields with a particles object
  !!
  !! Call the register_data of all potentials.
  !<
  subroutine potentials_register_data(this_cptr, p, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(inout) :: p
    integer, optional, intent(out)   :: ierror

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(this_cptr, this)

#define REGISTER_DATA(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call register_data(this%x(i), p, ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    REGISTER_DATA({classname})

#undef REGISTER_DATA

  endsubroutine potentials_register_data

  !>
  !! Constructor
  !!
  !! Call the constructors of all potentials, and
  !! removes the respective lists from memory.
  !!
  !! Note: This is used by the standalone code only.
  !<
  subroutine potentials_init(this_cptr)
    implicit none

    type(C_PTR), intent(in)  :: this_cptr

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    call c_f_pointer(this_cptr, this)

#define INIT(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call init(this%x(i)) ; enddo ; endif

    INIT({classname})

#undef INIT

  endsubroutine potentials_init


  !>
  !! Destructor
  !!
  !! Call the destructors of all potentials, and
  !! removes the respective lists from memory.
  !<
  subroutine potentials_del(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    call c_f_pointer(this_cptr, this)

#define DEL(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call del(this%x(i)) ; enddo ; deallocate(this%x) ; endif

    DEL({classname})

#undef DEL

  endsubroutine potentials_del


  !>
  !! Set the Coulomb solver
  !!
  !! Set the Coulomb solver
  !<
  subroutine potentials_set_Coulomb(this_cptr, coul, ierror)
    use, intrinsic :: iso_c_binding
   
    implicit none

    type(C_PTR),        intent(in)  :: this_cptr
    type(C_PTR),        intent(in)  :: coul
    integer,  optional, intent(out) :: ierror

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(this_cptr, this)

#define SET_COULOMB(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call set_Coulomb(this%x(i), coul, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    SET_COULOMB({classname})
    
#undef SET_COULOMB

  endsubroutine potentials_set_Coulomb



  !>
  !! Bind the potentials to a certain Particles and Neighbors object
  !!
  !! Bind the potentials to a certain Particles and Neighbors object. This will
  !! tell the potential to initialize its internal buffers according to the 
  !! array sizes used by the Particles and Neighbors object. All subsequent
  !! calls to energy and forces *must* be carried out with the same Particles
  !! and Neighbors object. If either Particles or Neighbors object changes,
  !! bind_to will need to be called again.
  !<
  subroutine potentials_bind_to(this_cptr, p, nl, coul, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    type(C_PTR),        intent(in)    :: coul
    integer,  optional, intent(out)   :: ierror

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(this_cptr, this)

#define BIND_TO(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call bind_to(this%x(i), p, nl, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    BIND_TO({classname})

#undef BIND_TO

    if (coulomb_is_enabled(coul)) then
       call coulomb_bind_to(coul, p, nl, ierror=ierror)
       PASS_ERROR(ierror)
    endif

#define BIND_TO_WITH_COUL(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call bind_to_with_coul(this%x(i), p, nl, coul, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    BIND_TO_WITH_COUL({classname})

#undef BIND_TO_WITH_COUL

  endsubroutine potentials_bind_to


  !>
  !! Compute energies and forces
  !!
  !! Calls the energy and forces of all potentials.
  !<
  subroutine potentials_energy_and_forces(this_cptr, dyn, nl, coul, q, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),           intent(in)    :: this_cptr
    type(dynamics_t),      intent(inout) :: dyn
    type(neighbors_t),     intent(inout) :: nl
    type(C_PTR), optional, intent(in)    :: coul
    real(DP),    optional, intent(inout) :: q(dyn%p%maxnatloc)
    integer,     optional, intent(out)   :: ierror

    ! ---

    type(potentials_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(this_cptr, this)

    dyn%epot  = 0.0_DP
    dyn%f     = 0.0_DP
    dyn%wpot  = 0.0_DP

    if (associated(this%epot_per_at))    this%epot_per_at    = 0.0_DP
    if (associated(this%wpot_per_at))    this%wpot_per_at    = 0.0_DP

    ! ---

! In Fortran 2008 we can pass unallocated array as non-existent optional
! arguments. Let's hope all compilers support this.

#define ENERGY_AND_FORCES(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call energy_and_forces(this%x(i), dyn%p, nl, dyn%epot, dyn%f, dyn%wpot, epot_per_at=this%epot_per_at, wpot_per_at=this%wpot_per_at, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    ENERGY_AND_FORCES({classname})

#undef ENERGY_AND_FORCES

#define ENERGY_AND_FORCES_WITH_CHARGES(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call energy_and_forces_with_charges(this%x(i), dyn%p, nl, dyn%epot, dyn%f, dyn%wpot, q=q, epot_per_at=this%epot_per_at, wpot_per_at=this%wpot_per_at, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    ENERGY_AND_FORCES_WITH_CHARGES({classname})

#undef ENERGY_AND_FORCES_WITH_CHARGES

#define ENERGY_AND_FORCES_WITH_CHARGES_AND_COUL(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call energy_and_forces_with_charges_and_coul(this%x(i), dyn%p, nl, q, coul, dyn%epot, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    ENERGY_AND_FORCES_WITH_CHARGES_AND_COUL({classname})

#undef ENERGY_AND_FORCES_WITH_CHARGES_AND_COUL

#define ENERGY_AND_FORCES_WITH_DYN(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call energy_and_forces_with_dyn(this%x(i), dyn, nl, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    ENERGY_AND_FORCES_WITH_DYN({classname})

#undef ENERGY_AND_FORCES_WITH_DYN

    if (present(coul) .and. (present(q))) then
       if (coulomb_is_enabled(coul)) then
          call coulomb_energy_and_forces(coul, dyn%p, nl, q, dyn%epot, dyn%f, dyn%wpot, error=ierror)
          PASS_ERROR(ierror)
       endif
    endif

  endsubroutine potentials_energy_and_forces

endmodule potentials

