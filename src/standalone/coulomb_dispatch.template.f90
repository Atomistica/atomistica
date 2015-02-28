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

    INIT_ERROR(ierror)

    call c_f_pointer(this_cptr, this)

    phi = 0.0_DP

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
  subroutine coulomb_energy_and_forces(this_cptr, p, nl, q, epot_out, f_out, &
       wpot_out, error)
    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(particles_t),  intent(in)    :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(in)    :: q(p%nat)
    real(DP),           intent(inout) :: epot_out
    real(DP),           intent(inout) :: f_out(3, p%nat)
    real(DP),           intent(inout) :: wpot_out(3, 3)
    integer,  optional, intent(out)   :: error

    ! ---

    type(coulomb_t), pointer :: this

    real(DP) :: f(3, p%nat)

    ! ---

    INIT_ERROR(error)

    call c_f_pointer(this_cptr, this)

    this%epot = 0.0_DP
    f         = 0.0_DP
    this%wpot = 0.0_DP

#define ENERGY_AND_FORCES(x)  if (allocated(this%x)) then ; call energy_and_forces(this%x, p, nl, q, this%epot, f, this%wpot, error=error) ; PASS_ERROR(error) ; endif

    ENERGY_AND_FORCES({classname})

#undef ENERGY_AND_FORCES
 
    if (system_of_units == eV_A .or. system_of_units == eV_A_fs) then
       this%epot = this%epot * Hartree*Bohr
       this%wpot = this%wpot * Hartree*Bohr

       f_out = f_out + Hartree*Bohr*f
    else
       f_out = f_out + f
    endif

    epot_out = epot_out + this%epot
    wpot_out = wpot_out + this%wpot

  endsubroutine coulomb_energy_and_forces

endmodule coulomb

