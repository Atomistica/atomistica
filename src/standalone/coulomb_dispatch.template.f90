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

  use data
  use particles
  use neighbors

  use {classname}
  
  implicit none

  private

  character(MAX_NAME_STR), parameter :: PHI_STR  = "electrostatic_potential"
  character(MAX_NAME_STR), parameter :: E_STR    = "electric_field"

  integer, parameter  :: PHI_TAG  = F_TO_TRAJ
  integer, parameter  :: E_TAG    = F_TO_TRAJ

  public :: coulomb_t
  type coulomb_t

     integer :: p_pos_rev  = -1     !< Last revision of the particles object, to detect changes
     integer :: p_other_rev  = -1   !< Last revision of the particles object, to detect changes

     !
     ! Dispatch table
     !

     type({classtype}), allocatable :: {classname}

     !
     ! Internal potential and electric field arrays
     !

     real(DP), pointer :: phi(:)               !< Electrostatic potential
     real(DP), pointer :: E(:, :)              !< Electric field

     real(DP) :: epot = 0.0_DP        !< Coulomb interaction energy
     real(DP) :: wpot(3, 3) = 0.0_DP  !< Virial

  endtype coulomb_t

  ! Note: coulomb_t is hidden. Everything is passed as type(C_PTR) to hide the
  ! complexity of coulomb_t from the compiler. This speeds up compile times
  ! and avoids nasty compiler crashes. However, this invalidates Fortran
  ! interfaces since the compiler can't match a generic call to datatype.

  public :: C_PTR

  public :: coulomb_alloc, coulomb_free, coulomb_register_data
  public :: coulomb_init, coulomb_is_enabled
  public :: coulomb_del, coulomb_bind_to, coulomb_set_Hubbard_U
  public :: coulomb_potential, coulomb_potential_and_field
  public :: coulomb_energy_and_forces, coulomb_get_potential_energy

contains

  !>
  !! Allocator
  !!
  !! Allocate memory for new coulomb instance
  !<
  subroutine coulomb_alloc(this_cptr)
    implicit none

    type(C_PTR), intent(out) :: this_cptr

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    allocate(this)
    this_cptr = c_loc(this)

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

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    deallocate(this)

  endsubroutine coulomb_free


  !>
  !! Constructor
  !!
  !! Register the phi and E fields
  !<
  subroutine coulomb_register_data(this_cptr, p)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(inout) :: p

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    call add_real(p%data, PHI_STR, PHI_TAG)
    call add_real3(p%data, E_STR, E_TAG)

  endsubroutine coulomb_register_data


  !>
  !! Constructor
  !!
  !! Call the constructor of all allocated Coulomb objects
  !<
  subroutine coulomb_init(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

#define INIT(x)  if (allocated(this%x)) then ; call init(this%x) ; endif

    INIT({classname})

#undef INIT

  endsubroutine coulomb_init


  !>
  !! Check whether any Coulomb module is enabled
  !!
  !! Check whether any Coulomb module is enabled
  !<
  logical function coulomb_is_enabled(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    coulomb_is_enabled  = any( (/ &
         allocated(this%{classname}), &
         .false. &
         /) )

  endfunction coulomb_is_enabled


  !>
  !! Destructor
  !!
  !! Delete the Coulomb dispatch object and all allocated objects driver
  !<
  subroutine coulomb_del(this_cptr)
    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

#define DEL(x)  if (allocated(this%x)) then ; call del(this%x) ; deallocate(this%x) ; endif

    DEL({classname})

#undef DEL

!    if (allocated(this%phi)) then
!       deallocate(this%phi)
!    endif
!    if (allocated(this%E)) then
!       deallocate(this%E)
!    endif

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

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    call ptr_by_name(p%data, PHI_STR, this%phi, ierror)
    PASS_ERROR(ierror)
    call ptr_by_name(p%data, E_STR, this%E, ierror)
    PASS_ERROR(ierror)

#define BIND_TO(x)  if (allocated(this%x)) then ; call bind_to(this%x, p, nl, ierror) ; PASS_ERROR(ierror) ; endif

    BIND_TO({classname})

#undef BIND_TO

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

    type(coulomb_t), pointer :: this

    ! ---

    INIT_ERROR(ierror)

    call c_f_pointer(this_cptr, this)

#define SET_HUBBARD_U(x)  if (allocated(this%x)) then ; call set_Hubbard_U(this%x, p, U, error=ierror) ; PASS_ERROR(ierror) ; endif
 
    SET_HUBBARD_U({classname})

#undef SET_HUBBARD_U

  endsubroutine coulomb_set_Hubbard_U


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !!
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !<
  subroutine coulomb_potential_and_field(this_cptr, p, nl, q, epot, wpot, phi, &
       E, ierror)
    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(particles_t),  intent(in)    :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(in)    :: q(p%maxnatloc)
    real(DP), optional, intent(inout) :: phi(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot
    real(DP), optional, intent(inout) :: E(3, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot(3, 3)
    integer,  optional, intent(inout) :: ierror

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    if (have_positions_changed(p, this%p_pos_rev) .or. has_other_changed(p, this%p_other_rev)) then

       this%epot              = 0.0_DP
       this%wpot              = 0.0_DP

       ! Note: This needs to be 1:p%nat, NOT 1:p%natloc for the fast-multipole solver
       this%phi(1:p%nat)      = 0.0_DP
       VEC3(this%E, 1:p%nat)  = 0.0_DP

#define POTENTIAL_AND_FIELD(x)  if (allocated(this%x)) then ; call potential_and_field(this%x, p, nl, q, this%phi, this%epot, this%E, this%wpot, ierror) ; PASS_ERROR(ierror) ; endif

       POTENTIAL_AND_FIELD({classname})

#undef POTENTIAL_AND_FIELD

       if (system_of_units == eV_A .or. system_of_units == eV_A_fs) then
          this%phi    = this%phi * Hartree * Bohr
          this%epot   = this%epot * Hartree * Bohr

          this%E      = this%E * Hartree * Bohr
          this%wpot   = this%wpot * Hartree * Bohr
       endif

       !this%epot = 0.5_DP*dot_product(q(1:p%natloc), this%phi(1:p%natloc))
       !write (*, *)  this%epot, 0.5_DP*dot_product(q(1:p%natloc), this%phi(1:p%natloc))

    endif

    if (present(epot)) then
       epot                 = this%epot
    endif
    if (present(wpot)) then
       wpot                 = this%wpot
    endif
    if (present(phi)) then
       phi(1:p%natloc)      = this%phi(1:p%natloc)
    endif
    if (present(E)) then
       VEC3(E, 1:p%natloc)  = VEC3(this%E, 1:p%natloc)
    endif

  endsubroutine coulomb_potential_and_field


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

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    !call potential_and_field(this, p, nl, q, phi=phi, ierror=ierror)
    !PASS_ERROR(ierror)

    phi  = 0.0_DP

#define POTENTIAL(x)  if (allocated(this%x)) then ; call potential(this%x, p, nl, q, phi, ierror) ; PASS_ERROR(ierror) ; endif

    POTENTIAL({classname})

    if (system_of_units == eV_A .or. system_of_units == eV_A_fs) then
       phi  = phi * Hartree * Bohr
    endif

#undef POTENTIAL

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

    type(coulomb_t), pointer :: this

    ! ---

    integer   :: i

    ! ---

    call c_f_pointer(this_cptr, this)
    ! FIXME! Pressure tensor does not work yet, only
    ! hydrostatic pressure works.

    call coulomb_potential_and_field(this_cptr, p, nl, q, ierror=ierror)
    PASS_ERROR(ierror)

    if (present(f)) then
       do i = 1, p%natloc
          VEC3(f, i)  = VEC3(f, i) + q(i)*VEC3(this%E, i)
       enddo
    endif

    epot     = epot + this%epot
    if (present(wpot)) then
       wpot  = wpot + this%wpot
    endif

  endsubroutine coulomb_energy_and_forces


  !>
  !! Return potential energy
  !!
  !! Returns previously calculated potential energy
  !<
  real(DP) function coulomb_get_potential_energy(this_cptr)
    implicit none

    type(C_PTR),        intent(in)    :: this_cptr

    ! ---

    type(coulomb_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    coulomb_get_potential_energy = this%epot

  endfunction coulomb_get_potential_energy

endmodule coulomb

