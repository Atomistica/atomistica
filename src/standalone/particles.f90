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
!! Particle information
!!
!! Position and cell information are stored in the data structures
!! of this module. Velocity and force information are stored in a different
!! one, see dynamics.f90.
!<
module particles
  use, intrinsic :: iso_c_binding

  use supplib

  implicit none

  private

  !>
  !! Highest element number stored in the periodic table module
  !<
  integer, parameter :: MAX_Z = ubound(ElementName, 1)
  public :: MAX_Z

  public :: TI_ATTR_STR, Z_STR, EL_STR, INDEX_STR, M_STR, R_NON_CYC_STR, R_CONT_STR, G_STR, V_STR, F_STR, Q_STR
  public :: SHEAR_DX_STR, CELL_STR

  character(MAX_NAME_STR), parameter  :: TI_ATTR_STR      = "time"

  character(MAX_NAME_STR), parameter  :: Z_STR            = "Z"
  character(MAX_NAME_STR), parameter  :: Z_ALIAS_STR      = "atom_types"
  character(MAX_NAME_STR), parameter  :: EL_STR           = "internal_element_number"
  character(MAX_NAME_STR), parameter  :: INDEX_STR        = "atom_index"
  character(MAX_NAME_STR), parameter  :: M_STR            = "masses"
  character(MAX_NAME_STR), parameter  :: R_NON_CYC_STR    = "coordinates"                  ! ... are allowed to leave the cell between neighbor list updates
  character(MAX_NAME_STR), parameter  :: R_CONT_STR       = "continuous_coordinates"     ! ... will never be wrapped back to the cell
  character(MAX_NAME_STR), parameter  :: G_STR            = "groups"

  ! Managed by the "dynamics" object
  character(MAX_NAME_STR), parameter  :: V_STR  = "velocities"
  character(MAX_NAME_STR), parameter  :: F_STR  = "forces"

  ! Charges are "managed" by the Coulomb objects
  character(MAX_NAME_STR), parameter  :: Q_STR            = "charges"

  character(MAX_NAME_STR), parameter  :: SHEAR_DX_STR     = "shear_dx"
  character(MAX_NAME_STR), parameter  :: CELL_STR         = "cell"

  public :: F_CONSTANT, F_VERBOSE_ONLY, F_RESTART, F_TO_TRAJ, F_COMMUNICATE, F_COMM_GHOSTS, F_COMM_FORCES, F_TO_ENER, F_ALL, Q_TAG

  integer, parameter  :: F_CONSTANT      = 1         ! Field does not vary over time
  integer, parameter  :: F_VERBOSE_ONLY  = 2         ! Internal use only
  integer, parameter  :: F_RESTART       = 4         ! Necessary for a clean restart
  integer, parameter  :: F_TO_TRAJ       = 8         ! Output to trajectory file
  integer, parameter  :: F_COMMUNICATE   = 16        ! Communicate this field
  integer, parameter  :: F_COMM_GHOSTS   = 32        ! Communicate this field for ghost particles
  integer, parameter  :: F_COMM_FORCES   = 64        ! Communicate this property back to the ghost particle
  integer, parameter  :: F_TO_ENER       = 128       ! Output to ener.out file

  integer, parameter  :: F_ALL  = F_CONSTANT + F_VERBOSE_ONLY + F_RESTART + F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS
  
  integer, parameter  :: Q_TAG  = F_TO_TRAJ + F_RESTART + F_COMMUNICATE + F_COMM_GHOSTS

  
  !
  ! This stores the static information,
  ! i.e. the *positions*
  !
  
  public :: particles_t
  type particles_t

     !
     ! Is this particles-object initialized?
     !

     logical                :: initialized  = .false.

     integer                :: pos_rev          = 0          !> Have the positions been changed?
     integer                :: other_rev        = 0          !> Has anything else been changed?

     !
     ! Simulation box
     !

     logical                :: orthorhombic_cell_is_required = .false.
     logical                :: cell_is_orthorhombic

     real(DP), pointer      :: Abox(:, :)
     real(DP)               :: Bbox(3, 3)

     !
     ! Simulation box (on this processor only)
     !

     real(DP)               :: lower(3)
     real(DP)               :: upper(3)
     
     real(DP)               :: lower_with_border(3)
     real(DP)               :: upper_with_border(3)

     !>
     !! Communication border
     !<
     real(DP)               :: border = 0.0_DP
     
     !
     ! Periodicity
     !

     integer :: pbc(3)
     logical :: locally_pbc(3)

     !
     ! Lees-Edwards boundary conditions
     !

     real(DP), pointer      :: shear_dx(:)
     real(DP)               :: shear_dv(3)

     !
     ! Interaction range
     !

     real(DP)               :: accum_max_dr

     !
     ! Degrees of freedom (=totnat in the unconstrained case)
     !

     integer                :: dof

     !
     ! Particle number information
     !

     integer                :: nat                  ! number of particles in system
                                                    ! (including ghost particles)
     integer                :: natloc               ! number of particles on this processor
                                                    ! (excluding ghost particles)
     integer                :: maxnatloc            ! maximum number of particles on this processor
     integer                :: totnat               ! total number of particles on all processors

     !
     ! All particel data is managed by the *data* field. The other fields are pointers to the
     ! entries of data.
     !

     type(data_t), pointer  :: data

     character(4), pointer  :: sym(:)
     integer, pointer       :: Z(:)                 ! element number
     integer, pointer       :: el(:)                ! element number
     integer, pointer       :: index(:)             ! global index
     real(DP), pointer      :: m(:)                 ! mass

     ! These position are always global, but inside the box and may
     ! be outside the local box (MPI version).
     ! These are identical for all processes.
#ifndef IMPLICIT_R
     real(DP), pointer      :: r(:, :)              ! positions
#endif

     ! These positions are always local and may be outside the global box.
     ! These differ on for different processes.
     real(DP), pointer      :: r_non_cyc(:, :)      ! displacement from last binning

     ! These positions are continouos from the start of the simulation.
     ! No pbc boundary condition are applied.
     real(DP), pointer      :: r_cont(:, :)         ! displacement from beginning of simulation

     integer, pointer       :: g(:)                 !< group (for the user's bookkeeping, should not be touched in the code)

     real(DP), allocatable  :: sort_index(:)

#ifdef _MP
     integer, allocatable   :: from_rank(:)         ! Ghost particles: Where do they come from
#endif

     !
     ! Global to local index transformation.
     ! Needed for the next list.
     ! (Tommi Jaervi, 13.11.2009: Moved next list to molecules module
     !  but probably this transformation is needed for MPI anyway?)
     !

     integer, allocatable   :: global2local(:)

     !
     ! Some statistics, i.e. which elements occur in the simulation
     !

     integer                :: nZ(MAX_Z)
     integer                :: nel           !> number of distinct elements
     integer                :: el2Z(MAX_Z)   !> id - i.e. from 1 to nel
     integer                :: Z2el(MAX_Z)   !> reverse mapping
     
     !
     ! Tag - this is used to attach the umbrella Python instance
     !
     
     type(C_PTR) :: tag

  endtype particles_t

!  integer, allocatable   :: g_index(:)

  !
  ! Some unit conversion stuff
  !

  public :: eVA_to_fs_sq, sqrt_eVA_to_fs_sq

  real(DP), parameter  :: eVA_to_fs_sq = 1000.0_DP/(6.0221415_DP*1.60217646_DP) ! femtoseconds
  real(DP), parameter  :: sqrt_eVA_to_fs_sq = 10.1805056398418_DP

  !
  ! The system of units to be used
  !

  public :: eV_A, eV_A_fs, H_Bohr, n_units, len_unit_str, STR_eVA, STR_HBohr
  public :: unit_strs

  integer, parameter  :: NA            = 0
  integer, parameter  :: eV_A          = 1
  integer, parameter  :: eV_A_fs       = 2
  integer, parameter  :: H_Bohr        = 3
  integer, parameter  :: n_units       = 4
  integer, parameter  :: len_unit_str  = 20

  ! This is need for xlf
  character(len_unit_str), parameter  :: STR_NA     = CSTR("N/A")
  character(len_unit_str), parameter  :: STR_eVA    = CSTR("eV/A")
  character(len_unit_str), parameter  :: STR_eVAfs  = CSTR("eV/A/fs")
  character(len_unit_str), parameter  :: STR_HBohr  = CSTR("H/Bohr")
  character(len_unit_str), parameter  :: unit_strs(n_units) = &
       (/ STR_NA, STR_eVA, STR_eVAfs, STR_HBohr /)

  public :: system_of_units, centered_box, length_to_A, length_to_Bohr
  public :: energy_to_eV, K_to_energy
  public :: time_to_fs, velocity_to_Afs, pressure_to_GPa
  public :: mass_to_g_mol

  integer,       target :: system_of_units = NA
  logical(BOOL), target :: centered_box = .false.
  
  real(DP)            :: length_to_A      = 1.0_DP
  real(DP)            :: length_to_Bohr   = 1.0_DP/Bohr
  real(DP)            :: energy_to_eV     = 1.0_DP
  real(DP)            :: K_to_energy      = Boltzmann_K
  real(DP)            :: time_to_fs       = sqrt_eVA_to_fs_sq
  real(DP)            :: velocity_to_Afs  = 1.0_DP/sqrt_eVA_to_fs_sq
  real(DP)            :: pressure_to_GPa  = elem_charge * 1d21
  real(DP)            :: mass_to_g_mol    = 1.0_DP

  public :: length_str, energy_str, time_str, force_str, pressure_str

  character(10)       :: length_str       = "A"
  character(10)       :: energy_str       = "eV"
  character(10)       :: time_str         = "10fs"
  character(10)       :: force_str        = "eV/A"
  character(10)       :: pressure_str     = "eV/A^3"
  character(10)       :: mass_str         = "g/mol"

  public :: init
  interface init
     module procedure particles_init, particles_init_from_particles
  endinterface

  public :: initialized
  interface initialized
     module procedure particles_initialized
  endinterface

  public :: allocate
  interface allocate
     module procedure particles_allocate
  endinterface

  public :: allocated
  interface allocated
     module procedure particles_allocated
  endinterface

  public :: del
  interface del
     module procedure particles_del
  endinterface

  public :: assign_ptrs
  interface assign_ptrs
     module procedure particles_assign_ptrs
  endinterface

!  interface align
!     module procedure particles_align
!  endinterface

  public :: center
  interface center
     module procedure particles_center
  endinterface

  public :: I_changed_positions
  interface I_changed_positions
     module procedure particles_I_changed_positions
  endinterface

  public :: have_positions_changed
  interface have_positions_changed
     module procedure particles_have_positions_changed
  endinterface

  public :: I_changed_other
  interface I_changed_other
     module procedure particles_I_changed_other
  endinterface

  public :: has_other_changed
  interface has_other_changed
     module procedure particles_has_other_changed
  endinterface

  public :: inbox
  interface inbox
     module procedure particles_inbox
  endinterface

  public :: update_elements
  interface update_elements
     module procedure particles_update_elements
  endinterface

  public :: compute_kinetic_energy_and_virial
  interface compute_kinetic_energy_and_virial
     module procedure particles_compute_kinetic_energy_and_virial
  endinterface

  public :: move
  interface move
     module procedure particles_move
  endinterface

  public :: pnc2pos
  interface pnc2pos
     module procedure particles_pnc2pos
  endinterface

  public :: remove
  interface remove
     module procedure particles_remove
  endinterface

  public :: require_orthorhombic_cell
  interface require_orthorhombic_cell
     module procedure particles_require_orthorhombic_cell
  endinterface

  public :: set_cell
  interface set_cell
     module procedure particles_set_cell, particles_set_cell_orthorhombic
  endinterface

  public :: set_lees_edwards
  interface set_lees_edwards
     module procedure particles_set_lees_edwards
  endinterface

  public :: set_from_Z
  interface set_from_Z
     module procedure particles_set_from_Z
  endinterface

  public :: set_total_nat
  interface set_total_nat
     module procedure particles_set_total_nat
  endinterface

  public :: sort
  interface sort
     module procedure sort_particles
  endinterface

  public :: swap
  interface swap
     module procedure swap_particles
  endinterface

  public :: translate
  interface translate
     module procedure particles_translate
  endinterface

  public :: volume
  interface volume
     module procedure particles_volume
  endinterface

  public :: operator(+)
  interface operator(+)
     module procedure particles_add
  endinterface

  public :: operator(*)
  interface operator(*)
     module procedure particles_mul
  endinterface

  public :: in_bounds
  interface in_bounds
     module procedure cyclic_in_bounds
  endinterface

  public :: in_cell
  interface in_cell
     module procedure cyclic_in_cell, cyclic_in_cell_vec
  endinterface

  public :: in_cellc
  interface in_cellc
     module procedure cyclic_in_cellc, cyclic_in_cellc_vec
  endinterface

  public :: request_border
  interface request_border
     module procedure particles_request_border
  endinterface request_border

  public :: get_true_cell
  interface get_true_cell
    module procedure particles_get_true_cell
  endinterface

  public :: units_init
  public :: particles_dump_info

contains

  !>
  !! Initially set/change cell size
  !!
  !! Initially set/change cell size
  !<
  subroutine particles_set_cell(this, Abox, pbc, scale_atoms, error)
    implicit none

    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: Abox(3, 3)
    logical, optional, intent(in)    :: pbc(3)
    logical, optional, intent(in)    :: scale_atoms
    integer, optional, intent(inout) :: error
    
    ! ---

    real(DP) :: A(3,3), fac(3, 3)
    integer  :: i, in, ipiv(3)

    ! --

!    call info("- particles_set_cell")

    if (present(pbc)) then
       where (pbc)
          this%pbc = 1
       elsewhere
          this%pbc = 0
       endwhere
    endif

    if (.not. all(this%pbc /= 0)) then
       call require_orthorhombic_cell(this, error)
       PASS_ERROR(error)
    endif

    if (present(scale_atoms) .and. scale_atoms) then
       fac = matmul(Abox, this%Bbox)
       !$omp  parallel do default(none) &
       !$omp& shared(this) firstprivate(fac)
       do i = 1, this%natloc
#ifndef IMPLICIT_R
          POS3(this, i) = matmul(fac, POS3(this, i))
#endif
          PNC3(this, i) = matmul(fac, PNC3(this, i))
          PCN3(this, i) = matmul(fac, PCN3(this, i))
       enddo
    endif

!    this%Abox  = ( Abox + transpose(Abox) )/2
    this%Abox  = Abox

    this%Bbox  = 0.0_DP
    do i = 1, 3
       this%Bbox(i, i)  = 1.0_DP
    enddo

    this%cell_is_orthorhombic = &
         abs(dot_product(this%Abox(1, :), this%Abox(2, :))) < 1d-12 .and. &
         abs(dot_product(this%Abox(2, :), this%Abox(3, :))) < 1d-12 .and. &
         abs(dot_product(this%Abox(3, :), this%Abox(1, :))) < 1d-12

!    if (.not. this%cell_is_orthorhombic) then
!       call info("     Cell is not orthorhombic.")
!    endif

!    call info("     " // this%Abox(1, :))
!    call info("     " // this%Abox(2, :))
!    call info("     " // this%Abox(3, :))

    A  = this%Abox
    call dgesv(3, 3, A, 3, ipiv, this%Bbox, 3, in)

    if (in /= 0) then
       RAISE_ERROR("Failed to determine the reciprocal lattice. Cell = " // this%Abox(:, 1) // ", " // this%Abox(:, 2) // ", " // this%Abox(:, 3), error)
    endif

    if (.not. all(this%pbc /= 0)) then
       call require_orthorhombic_cell(this, error)
       PASS_ERROR(error)
    endif

    if (.not. this%cell_is_orthorhombic .and. this%orthorhombic_cell_is_required) then
       RAISE_ERROR("This cell is non-orthorhombic, however, an orthorhombic cell was required.", error)
    endif

    this%lower  = (/ 0.0, 0.0, 0.0 /)
    this%upper  = (/ this%Abox(1, 1), this%Abox(2, 2), this%Abox(3, 3) /)

    this%lower_with_border = this%lower
    this%upper_with_border = this%upper

!    call info

  endsubroutine particles_set_cell


  !>
  !! Initially set/change cell size
  !!
  !! Initially set/change cell size
  !<
  subroutine particles_set_cell_orthorhombic(this, cell, pbc, scale_atoms, error)
    implicit none

    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: cell(3)
    logical, optional, intent(in)    :: pbc(3)
    logical, optional, intent(in)    :: scale_atoms
    integer, optional, intent(inout) :: error

    ! ---

    real(DP)  :: cell3x3(3, 3)

    ! ---

    cell3x3        = 0.0_DP
    cell3x3(1, 1)  = cell(1)
    cell3x3(2, 2)  = cell(2)
    cell3x3(3, 3)  = cell(3)

    call particles_set_cell(this, cell3x3, pbc=pbc, scale_atoms=scale_atoms, error=error)

  endsubroutine particles_set_cell_orthorhombic


  !**********************************************************************
  ! Set Lees-Edwards boundary conditions
  !**********************************************************************
  subroutine particles_set_lees_edwards(this, dx, dv, error)
    implicit none

    type(particles_t), intent(inout)  :: this
    real(DP), intent(in)              :: dx(3)
    real(DP), intent(in), optional    :: dv(3)
    integer, intent(inout), optional  :: error

    ! ---

    integer   :: i

    real(DP)  :: old_dx(3)

    ! ---

    call require_orthorhombic_cell(this, error)
    PASS_ERROR(error)

    old_dx         = this%shear_dx

    this%shear_dx  = dx

    if (present(dv)) then
       this%shear_dv  = dv
    endif

    do i = 1, 2
       do while (this%shear_dx(i) >= this%Abox(i, i)/2)
          this%shear_dx(i)  = this%shear_dx(i) - this%Abox(i, i)
       enddo
       do while (this%shear_dx(i) < -this%Abox(i, i)/2)
          this%shear_dx(i)  = this%shear_dx(i) + this%Abox(i, i)
       enddo
    enddo

    this%accum_max_dr  = this%accum_max_dr + norm( in_bounds(this, this%shear_dx - old_dx) )

    call I_changed_positions(this)
          
  endsubroutine particles_set_lees_edwards


  !**********************************************************************
  ! Python interface: Allocate a particle object
  !**********************************************************************
  subroutine particles_alloc(t)
    implicit none
    type(particles_t), pointer  :: t
    allocate(t)
  endsubroutine particles_alloc


  !**********************************************************************
  ! Python interface: Deallocate a particle object
  !**********************************************************************
  subroutine particles_dealloc(t)
    implicit none
    type(particles_t), pointer  :: t
    deallocate(t)
  endsubroutine particles_dealloc


  !>
  !! Initialize particle information
  !!
  !! Initialize particle information.
  !<
  subroutine particles_init(this)
    implicit none

    type(particles_t), intent(inout)   :: this

    ! ---

    this%initialized                    = .true.
    this%orthorhombic_cell_is_required  = .false.
    this%cell_is_orthorhombic           = .true.

    this%accum_max_dr      = 0.0_DP

    this%pbc          = (/ 1, 1, 1 /)
    this%locally_pbc  = (/ .true., .true., .true. /)

    this%border = 0.0_DP

    allocate(this%data)
    call init(this%data)

    call add_real3x3_attr( &
         this%data, &
         CELL_STR)
    call add_real3_attr( &
         this%data, &
         SHEAR_DX_STR)

    call add_integer( &
         this%data, &
         Z_STR, &
         alias=Z_ALIAS_STR, &
         tag=F_CONSTANT + F_TO_TRAJ )
    call add_integer( &
         this%data, &
         EL_STR, &
         F_CONSTANT + F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS )
#ifdef _MP
    call add_integer &
         (this%data, &
         INDEX_STR, &
         F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS )
#else
    call add_integer &
         (this%data, &
         INDEX_STR, &
         F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS )
#endif
    call add_real( &
         this%data, &
         M_STR, &
         F_CONSTANT + F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS )
#ifndef IMPLICIT_R
    call add_real3( &
         this%data, &
         R_STR, &
         F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS, &
         "angstroms", &
         length_to_A )
    call add_real3( &
         this%data, &
         R_NON_CYC_STR, &
         F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS, &
         "angstroms", &
         length_to_A )
#else
    call add_real3( &
         this%data, &
         R_NON_CYC_STR, &
         F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS, &
         "angstroms", &
         length_to_A )
#endif
    call add_real3( &
         this%data, &
         R_CONT_STR, &
         F_RESTART + F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS, &
         "angstroms", &
         length_to_A )
    call add_integer( &
         this%data, &
         G_STR, &
         F_CONSTANT + F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS )

  endsubroutine particles_init

  
  !**********************************************************************
  ! Allocate particle information
  !**********************************************************************
  subroutine particles_init_from_particles(this, from, error)
    implicit none

    type(particles_t), intent(inout)  :: this
    type(particles_t), intent(in)     :: from
    integer, intent(inout), optional  :: error

    ! ---

    this%initialized                    = .true.
    this%orthorhombic_cell_is_required  = from%orthorhombic_cell_is_required
    
    this%pbc          = (/ 1, 1, 1 /)
    this%locally_pbc  = (/ .true., .true., .true. /)

    this%border = 0.0_DP

    call init(this%data, from%data)

    call set_cell(this, from%Abox, from%pbc /= 0, error=error)

  endsubroutine particles_init_from_particles


  !**********************************************************************
  ! Allocate particle information
  !**********************************************************************
  logical function particles_initialized(p)
    implicit none

    type(particles_t), intent(in)  :: p

    ! ---

    particles_initialized  = p%initialized

  endfunction particles_initialized


  !>
  !! Allocate particle information
  !!
  !! Allocate particle information. This is also where all "per atom" data (particles%data)
  !! is allocated, so all data needed by other routines (such as molecules%next) should be
  !! registered.
  !!
  !! This means that one should call particles_init and others, such as dynamics_init and
  !! molecules_init, before calling particles_allocate.
  !<
  subroutine particles_allocate(this, nat, totnat, allow_def, error)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(in)               :: nat
    integer, intent(in), optional     :: totnat
    logical, intent(in), optional     :: allow_def
    integer, intent(inout), optional  :: error

    ! ---
    
    integer  :: i
    
    ! ---

    call allocate(this%data, nat, allow_def)

    this%nat            = nat
    this%natloc         = nat
    this%maxnatloc      = nat
    this%totnat         = nat
    this%dof            = 3*nat-3

    if (present(totnat)) then
       this%totnat  = totnat
    endif

    call particles_assign_ptrs(this)

    allocate(this%sym(nat))

    allocate(this%global2local(this%totnat))

    allocate(this%sort_index(nat))

#ifdef _MP
    allocate(this%from_rank(nat))
#endif

    do i = 1, nat
       this%index(i)         = i
       this%global2local(i)  = i
    enddo

    this%sym           = "H"
    this%Z             = 1
    this%g             = 1

#ifdef _MP
    this%global2local  = 0
#endif

    call set_cell(this, (/ 1.0_DP, 1.0_DP, 1.0_DP /), error=error)

    call update_elements(this)

  endsubroutine particles_allocate


  !>
  !! Check if the particles object has already been allocated
  !<
  function particles_allocated(this)
    implicit none

    type(particles_t), intent(in)  :: this
    logical                        :: particles_allocated

    ! ---

    particles_allocated  = allocated(this%data)

  endfunction particles_allocated


  !>
  !! Set the number of total particles in the simulation
  !!
  !! Set the number of total particles on all processors in this simulation.
  !! In particular, this will resize the *global2local* array.
  !<
  subroutine particles_set_total_nat(this, totnat)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(in)               :: totnat

    ! ---

    call prlog("- set_total_nat -")
    call prlog("     totnat  = " // totnat)

    call resize(this%global2local, totnat)
    this%totnat = totnat

    call prlog

  endsubroutine


  !>
  !! Destructor
  !!
  !! Remove this particles object from memory
  !<
  subroutine particles_del(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    this%initialized = .false.

    call del(this%data)
    deallocate(this%data)

    deallocate(this%sym)

    deallocate(this%global2local)
    deallocate(this%sort_index)

#ifdef _MP
    deallocate(this%from_rank)
#endif

  endsubroutine particles_del


  !**********************************************************************
  ! Assign shortcuts (i.e. r, v) to field in the *data* object
  !**********************************************************************
  subroutine particles_assign_ptrs(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    call attr_by_name(this%data, CELL_STR, this%Abox)
    call attr_by_name(this%data, SHEAR_DX_STR, this%shear_dx)

    call ptr_by_name(this%data, Z_STR, this%Z)
    call ptr_by_name(this%data, EL_STR, this%el)
    call ptr_by_name(this%data, INDEX_STR, this%index)
    call ptr_by_name(this%data, M_STR, this%m)
#ifndef IMPLICIT_R
    call ptr_by_name(this%data, R_STR, this%r)
#endif
    call ptr_by_name(this%data, R_NON_CYC_STR, this%r_non_cyc)
    call ptr_by_name(this%data, R_CONT_STR, this%r_cont)
    call ptr_by_name(this%data, G_STR, this%g)

  endsubroutine particles_assign_ptrs


  !**********************************************************************
  ! Copy particle f to t
  !**********************************************************************
  subroutine particles_move(this, t, f)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(in)               :: f
    integer, intent(in)               :: t

    ! ---

    integer  :: i

    ! --- 

    i  = this%index(f)

    this%sym(t)  = this%sym(f)
    call copy(this%data, t, f)

    this%global2local(i)  = t

  endsubroutine particles_move


  !**********************************************************************
  ! Deallocate particle information
  !**********************************************************************
  subroutine swap_particles(this, i, j)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(in)               :: i
    integer, intent(in)               :: j

    ! ---

    call swap(this%sym(i), this%sym(j))

    call swap(this%data, i, j)

    this%global2local(this%index(i)) = i
    this%global2local(this%index(j)) = j

  endsubroutine swap_particles


  !**********************************************************************
  ! Set masses and symbol from atomic number
  !**********************************************************************
  subroutine particles_set_from_Z(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    integer :: i

    ! ---

    do i = 1, this%nat
       if (this%Z(i) <= MAX_Z) then
          this%sym(i)  = ElementName(this%Z(i))
          this%m(i)    = ElementMass(this%Z(i))
       endif
    enddo

  endsubroutine particles_set_from_Z


  !**********************************************************************
  ! Compute statistics
  !**********************************************************************
  subroutine particles_update_elements(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    integer  :: i

#ifdef _MP
    type(MPI_context) :: mpi
#endif
    
    ! ---

    this%nZ  = 0

    do i = 1, this%natloc
       this%nZ(this%Z(i)) = this%nZ(this%Z(i))+1
    enddo

#ifdef _MP
    call initialise(mpi)
    call sum_in_place(mpi, this%nZ)
    call finalise(mpi)
#endif

    this%nel   = 0
    this%el2Z  = -1
    this%Z2el  = -1
    do i = 1, MAX_Z
       if (this%nZ(i) > 0) then
          this%nel             = this%nel+1
          this%Z2el(i)         = this%nel
          this%el2Z(this%nel)  = i
       endif
    enddo

    do i = 1, this%natloc
       this%el(i)  = this%Z2el(this%Z(i))
    enddo

  endsubroutine particles_update_elements

  
  !**********************************************************************
  ! Move all atoms that are outside the box inside.
  !**********************************************************************
  subroutine particles_inbox(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    integer            :: k, j

    real(DP), pointer  :: v(:, :)

    ! ---

    if (this%cell_is_orthorhombic) then
       if (any(this%shear_dx /= 0.0_DP)) then

          do k = 1, this%natloc
             
             do while (PNC(this, k, 3) < 0.0_DP)
                PNC3(this, k)    = PNC3(this, k)   + this%shear_dx
                PNC(this, k, 3)  = PNC(this, k, 3) + this%Abox(3, 3)
             enddo

             do while (PNC(this, k, 3) >= this%Abox(3, 3))
                PNC3(this, k)    = PNC3(this, k)   - this%shear_dx
                PNC(this, k, 3)  = PNC(this, k, 3) - this%Abox(3, 3)
             enddo

          enddo

       endif
       
       if (any(this%shear_dv /= 0.0_DP) .and. exists(this%data, V_STR)) then

          call ptr_by_name(this%data, V_STR, v)

          do k = 1, this%natloc

             do while (PNC(this, k, 3) < 0.0_DP)
                VEC3(v, k)       = VEC3(v, k) + this%shear_dv
             enddo

             do while (PNC(this, k, 3) >= this%Abox(3, 3))
                VEC3(v, k)       = VEC3(v, k) - this%shear_dv
             enddo

          enddo

       endif

       do j = 1, 3
          
          if (this%locally_pbc(j)) then
             do k = 1, this%natloc

                do while (PNC(this, k, j) < 0.0_DP)
                   PNC(this, k, j) = PNC(this, k, j) + this%Abox(j, j)
                enddo

                do while (PNC(this, k, j) >= this%Abox(j, j))
                   PNC(this, k, j) = PNC(this, k, j) - this%Abox(j, j)
                enddo
             
             enddo
          endif

       enddo

    else

       do j = 1, this%nat
          PNC3(this, j)  = cyclic_in_cell(this, PNC3(this, j))
       enddo

    endif

    ! Note: POS3 has different pbcity than PNC3
    call pnc2pos(this)

  endsubroutine particles_inbox


  !**********************************************************************
  ! Copy the non-cyclic coordinates to the one which are always inside
  ! the box.
  !**********************************************************************
  subroutine particles_pnc2pos(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

#ifndef IMPLICIT_R

    integer   :: k, j

    ! ---

    if (this%cell_is_orthorhombic) then

       do j = 1, this%natloc
          POS3(this, j)  = PNC3(this, j)
       enddo

       if (any(this%shear_dx /= 0.0_DP)) then

          do k = 1, this%natloc

             do while (POS(this, k, 3) < 0.0_DP)
                POS3(this, k)    = POS3(this, k)   + this%shear_dx
                POS(this, k, 3)  = POS(this, k, 3) + this%Abox(3, 3)
             enddo

             do while (POS(this, k, 3) >= this%Abox(3, 3))
                POS3(this, k)    = POS3(this, k)   - this%shear_dx
                POS(this, k, 3)  = POS(this, k, 3) - this%Abox(3, 3)
             enddo

          enddo

       endif

       do j = 1, 3
          
          if (this%pbc(j)) then
             do k = 1, this%natloc

                do while (POS(this, k, j) < 0.0_DP)
                   POS(this, k, j) = POS(this, k, j) + this%Abox(j, j)
                enddo

                do while (POS(this, k, j) >= this%Abox(j, j))
                   POS(this, k, j) = POS(this, k, j) - this%Abox(j, j)
                enddo
                
             enddo
          endif

       enddo

    else

       do j = 1, this%nat
          POS3(this, j)  = cyclic_in_cell(this, PNC3(this, j))
       enddo

    endif

#endif

  endsubroutine particles_pnc2pos


  !**********************************************************************
  ! Calculate the kinetic contribution to the pressure tensor
  !**********************************************************************    
  subroutine particles_compute_kinetic_energy_and_virial(this, v, f, wpot, ekin, fmax, wkin, pressure, mpi)
    implicit none

    type(particles_t), intent(inout)         :: this
    real(DP), intent(in)                     :: v(3, this%maxnatloc)
    real(DP), intent(in)                     :: f(3, this%maxnatloc)
    real(DP), intent(in)                     :: wpot(3, 3)
    real(DP), intent(out), optional          :: ekin
    real(DP), intent(out), optional          :: fmax
    real(DP), intent(out), optional          :: wkin(3, 3)
    real(DP), intent(out), optional          :: pressure(3, 3)
    type(MPI_context), intent(in), optional  :: mpi

    ! ---

    real(DP)  :: wkin_loc(3, 3)

    integer   :: i

    ! ---

    wkin_loc  = 0.0_DP
    do i = 1, this%natloc
       if (this%g(i) > 0) &
            wkin_loc = wkin_loc + this%m(i)*outer_product(VEC3(v, i), VEC3(v, i))
    enddo

    if (present(mpi)) then
       call sum_in_place(mpi, wkin_loc)
    endif

    if (present(wkin))  wkin  = wkin_loc

    if (present(ekin)) then
       ekin  = tr(3, wkin_loc)/2
    endif

    if (present(pressure)) then
       pressure  = ( wkin_loc - wpot )/volume(this)
    endif

    if (present(fmax)) then
       fmax  = 0.0_DP
       do i = 1, this%natloc
          if (this%g(i) > 0) then
             fmax = max(fmax, sqrt(dot_product(VEC3(f, i), VEC3(f, i))))
          endif
       enddo

       if (present(mpi)) then
          fmax = max(mpi, fmax)
       endif
    endif

  endsubroutine particles_compute_kinetic_energy_and_virial


  !**********************************************************************
  ! Dump information on particle *i* to log file
  !**********************************************************************
  subroutine particles_dump_info(this, i, cell)
    implicit none

    type(particles_t), intent(in)  :: this
    integer, intent(in)            :: i
    integer, intent(in), optional  :: cell(3)
    
    ! ---

    real(DP)  :: s(3)

    ! ---

    s    = matmul(this%Bbox, PNC3(this, i))
    s    = s - floor(s)

    write (ilog, *)
    write (ilog, '(A)')          "---"
    write (ilog, '(A, I15)')     "nat                = ", this%nat
    write (ilog, '(A, I15)')     "natloc             = ", this%natloc
    write (ilog, '(A, I15)')     "i                  = ", i
    write (ilog, '(A, I15)')     "index              = ", this%index(i)
    write (ilog, '(A, i15)')     "Z                  = ", this%Z(i)
    write (ilog, '(A, A)')       "symbol             = ", this%sym(i)
    write (ilog, '(A)')          "---"
    write (ilog, '(A, 3ES15.8)') "r                  = ", POS3(this, i)
    write (ilog, '(A, 3ES15.8)') "r_non_cyc          = ", PNC3(this, i)
    write (ilog, '(A, 3ES15.8)') "r_cont             = ", VEC3(this%r_cont, i)
    write (ilog, '(A, 3ES15.8)') "s                  = ", s
    if (present(cell)) then
       write (ilog, '(A, 3I15)') "cell               = ", cell
    endif
    write (ilog, '(A)')          "---"
    write (ilog, '(A, 3("/",F15.8,1X,"\",1X))')  "box vectors        = ", this%Abox(1, :)
    write (ilog, '(A, 3("|",F15.8,1X,"|",1X))')  "                     ", this%Abox(2, :)
    write (ilog, '(A, 3("\",F15.8,1X,"/",1X))')  "                     ", this%Abox(3, :)
    write (ilog, '(A)')          "---"
    write (ilog, '(A, 3F15.8)')  "lower              = ", this%lower
    write (ilog, '(A, 3F15.8)')  "upper              = ", this%upper
    write (ilog, '(A, 3F15.8)')  "lower_with_border  = ", this%lower_with_border
    write (ilog, '(A, 3F15.8)')  "upper_with_border  = ", this%upper_with_border
    write (ilog, '(A)')          "---"

  endsubroutine particles_dump_info


  !**********************************************************************
  ! Initialize units
  !**********************************************************************
  subroutine units_init(sou) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    ! ---

    integer(c_int), value  :: sou

    ! ---

    call prlog("- units_init -")

    system_of_units = sou

    if (system_of_units == eV_A) then
       call prlog("     Units are eV/A")

       length_to_A      = 1.0_DP
       length_to_Bohr   = 1.0/Bohr
       energy_to_eV     = 1.0_DP
       K_to_energy      = Boltzmann_K
       time_to_fs       = sqrt(eVA_to_fs_sq)
       velocity_to_Afs  = 1.0_DP/time_to_fs
       pressure_to_GPa  = elem_charge * 1d21
       mass_to_g_mol    = 1.0_DP

       length_str    = "A"
       energy_str    = "eV"
       time_str      = "10fs"
       force_str     = "eV/A"
       pressure_str  = "eV/A^3"
       mass_str      = "g/mol"

    else if (system_of_units == eV_A_fs) then
       call prlog("     Units are eV/A/fs")

       length_to_A      = 1.0_DP
       length_to_Bohr   = 1.0_DP/Bohr
       energy_to_eV     = 1.0_DP
       K_to_energy      = Boltzmann_K
       time_to_fs       = 1.0
       velocity_to_Afs  = 1.0_DP/time_to_fs
       pressure_to_GPa  = elem_charge * 1d21
       mass_to_g_mol    = 6.0221415_DP*1.60217646_DP*1d-3

       length_str    = "A"
       energy_str    = "eV"
       time_str      = "fs"
       force_str     = "eV/A"
       pressure_str  = "eV/A^3"
       mass_str      = "a.u."

    else if (system_of_units == H_Bohr) then
       call prlog("     Units are H/Bohr")

       length_to_A      = Bohr
       length_to_Bohr   = 1.0_DP
       energy_to_eV     = Hartree/Bohr
       K_to_energy      = Boltzmann_K/Hartree
       time_to_fs       = sqrt(evA_to_fs_sq) * Bohr/sqrt(Hartree)
       velocity_to_Afs  = 1.0_DP/time_to_fs
       pressure_to_GPa  = Hartree * elem_charge * 1d21 / ( Bohr**3 )
       mass_to_g_mol    = 1.0_DP

       length_str    = "Bohr"
       energy_str    = "H"
       time_str      = "fs"
       force_str     = "H/Bohr"
       pressure_str  = "H/Bohr^3"
       mass_str      = "g/mol"

    else
       write (*, '(A,I5)')  "[main_loop] Unknown system of units: ", &
            system_of_units
       stop
    endif

    ! Scale masses
    ElementMass = ElementMass_in_g_mol/mass_to_g_mol

    call prlog

  endsubroutine units_init


  !*******************************************************************
  ! Heap sort
  !   Quickly (n*log(n)) sorts particles according to sort index
  !*******************************************************************
  subroutine sort_particles(this, sort_index)
    implicit none

    type(particles_t), intent(inout)  :: this
    real(DP), intent(inout)           :: sort_index(this%nat)

    ! ---

    integer   :: i, j, root

    ! ---

    do i = this%natloc/2, 1, -1

       j = i
       call siftdown(j, this%natloc)

    enddo

    do i = this%natloc, 2, -1

       call swap(sort_index(1), sort_index(i))

       call swap_particles(this, 1, i)

       root = 1
       j = i - 1
       call siftdown(root, j)

    enddo

  contains

    subroutine siftdown(root, bottom)
      implicit none

      integer, intent(inout)  :: root, bottom
      
      ! ---

      integer                 :: done, maxchild

      ! ---

      done = 0

      do while ((root*2 <= bottom) .and. done == 0)

         if (root*2 == bottom) then
            maxchild = root * 2
         else if (sort_index(root*2) > sort_index(root*2+1)) then
            maxchild = root * 2
         else
            maxchild = root*2 + 1
         endif

         if (sort_index(root) < sort_index(maxchild)) then

            call swap(sort_index(root), sort_index(maxchild))

            call swap_particles(this, root, maxchild)

            root = maxchild
         else
            done = 1
         endif

      enddo

    endsubroutine siftdown

  endsubroutine sort_particles


  !*******************************************************************
  ! Check if global2local list is correct (DEBUG)
  !*******************************************************************
  subroutine check_global2local(this)
    implicit none

    type(particles_t), intent(in)  :: this

    ! ---

    integer  :: i, j, l
    logical  :: found

    ! ---

    do i = 1, this%totnat
       l = this%global2local(i)

       if (l > 0) then

          found = .false.

          do j = 1, this%nat
             if (this%index(j) == i) then

                if (l == j) then
                   found = .true.
                else
                   stop "l /= j"
                endif

             endif
          enddo

          if (.not. found) then
             stop "not found 1"
          endif

       endif

    enddo

    do i = 1, this%nat
       if (this%global2local(this%index(i)) /= i) then
          write (ilog, *)  i, this%index(i), this%global2local(this%index(i))
          stop "this%global2local(this%index(i)) /= i"
       endif
    enddo

  endsubroutine check_global2local


  !**********************************************************************
  ! Join two particle objects, the cell is taken from the first argument
  !**********************************************************************
  function particles_add(in_this, in_that) result(out)
    implicit none

    type(particles_t), intent(in)     :: in_this
    type(particles_t), intent(in)     :: in_that
    type(particles_t)                 :: out

    ! ---

    integer            :: i
    logical            :: del_this, del_that
    type(particles_t)  :: this, that

    ! ---

    this      = in_this
    that      = in_that

    del_this  = .false.
    del_that  = .false.

!    write (*, *)  out%sym(1)
!    write (*, *)  this%sym(1)

    if (associated(this%sym, out%sym)) then
       write (*, *)  "this == out"
       del_this  = .true.
    endif

    if (associated(that%sym, out%sym)) then
       write (*, *)  "this == out"
       del_that  = .true.
    endif

    call init(out, this)
    call allocate(out, this%natloc+that%natloc)

    call set_cell(out, this%Abox, this%pbc /= 0)

    out%sym(1:this%natloc)                          = this%sym(1:this%natloc)
    out%sym(this%natloc+1:this%natloc+that%natloc)  = that%sym(1:that%natloc)

    call copy(out%data, 1,             this%natloc,             this%data, 1, this%natloc)
    call copy(out%data, this%natloc+1, this%natloc+that%natloc, that%data, 1, that%natloc)

    do i = 1, that%natloc
       out%index(this%natloc+i)                    = that%index(i) + this%natloc
       out%global2local(out%index(this%natloc+i))  = i             + this%natloc
    enddo

    if (del_this) then
       call del(this)
    endif

    if (del_that) then
       call del(that)
    endif

  endfunction particles_add


  !**********************************************************************
  ! Create a supercell from the current particle configuration
  !**********************************************************************
  function particles_mul(this, arg) result(out)
    implicit none

    type(particles_t), intent(in)     :: this
    integer, intent(in)               :: arg(3)
    type(particles_t)                 :: out

    ! ---

    real(DP)  :: dr(3), cell(3, 3)
    integer   :: i, j, k, l, m

    ! ---

    call init(out, this)
    call allocate(out, this%natloc*arg(1)*arg(2)*arg(3))

    cell(:, 1)  = this%Abox(:, 1)*arg(1)
    cell(:, 2)  = this%Abox(:, 2)*arg(2)
    cell(:, 3)  = this%Abox(:, 3)*arg(3)
    call set_cell(out, cell, this%pbc /= 0)

    l = 1
    do i = 0, arg(1)-1
       do j = 0, arg(2)-1
          do k = 0, arg(3)-1
             out%sym(l:l+this%natloc-1)  = this%sym(1:this%natloc)
             call copy(out%data, l, l+this%natloc-1, this%data, 1, this%natloc)

             out%index(l:l+this%natloc-1)         = this%index(1:this%natloc) + l-1
             out%global2local(l:l+this%natloc-1)  = this%global2local(1:this%natloc) + l-1

             dr                                   = matmul(this%Abox, (/ i, j, k /))

             do m = l, l+this%natloc-1
#ifndef IMPLICIT_R
                POS3(out, m)  = POS3(out, m) + dr
#endif
                PNC3(out, m)  = PNC3(out, m) + dr
             enddo

             l  = l + this%natloc
          enddo 
       enddo
    enddo

  endfunction particles_mul


  !**********************************************************************
  ! Remove a particle consistently
  !**********************************************************************
  subroutine particles_remove(this, i, error)
    implicit none

    type(particles_t), intent(inout)   :: this
    integer,           intent(in)      :: i
    integer, optional, intent(out)     :: error

    ! ---

    integer  :: removed_global_index

    ! ---

    INIT_ERROR(error)

    removed_global_index  = this%index(i)

    ASSERT(this%global2local(removed_global_index) == i, "particles_remove: Particle not on this processor.", error)

    this%sym(i)  = this%sym(this%natloc)
    call copy(this%data, i, this%natloc)
    this%global2local(this%index(i))  = i

    ! Update index and global2local
    if (removed_global_index /= this%natloc) then
       this%index(this%global2local(this%natloc))  = removed_global_index
       this%global2local(removed_global_index)     = this%global2local(this%natloc)
    endif

    this%natloc  = this%natloc-1
    this%nat     = this%natloc

  endsubroutine particles_remove


  !**********************************************************************
  ! Call to require an orthorhombic cell
  !**********************************************************************
  subroutine particles_require_orthorhombic_cell(this, error)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(inout), optional  :: error

    ! ---

    this%orthorhombic_cell_is_required  = .true.
 
   if (.not. this%cell_is_orthorhombic) then
       RAISE_ERROR("Orthorhombic cell is requested, however cell is already non-orthorhombic.", error)
    endif

  endsubroutine particles_require_orthorhombic_cell


  !**********************************************************************
  ! Align the eigenvectors of the inertia tensor along the x-, y- and
  ! z-axis.
  !**********************************************************************
!!$  subroutine particles_align(this)
!!$    implicit none
!!$
!!$    type(particles_t), intent(inout)   :: this
!!$
!!$    ! ---
!!$
!!$    real(DP)  :: it(3, 3), Ixx, Iyy, Izz, Ixy, Ixz, Iyz
!!$
!!$    integer   :: i
!!$
!!$    ! ---
!!$
!!$    call particles_center(this)
!!$
!!$    it(:, :)  = 0.0_DP
!!$    do i = 1, p%natloc
!!$       it(
!!$    enddo
!!$
!!$  endsubroutine particles_align


  !**********************************************************************
  ! Center around zero
  !**********************************************************************
  subroutine particles_center(this, vacuum, cell, error)
    implicit none

    type(particles_t), intent(inout)  :: this
    real(DP), intent(in), optional    :: vacuum(3)
    real(DP), intent(in), optional    :: cell(3, 3)
    integer, intent(inout), optional  :: error

    ! ---

    real(DP)  :: com(3), box(3)

    integer   :: i

    ! ---

    com  = 0.0_DP
    do i = 1, this%natloc
       com  = com + this%m(i)*POS3(this, i)
    enddo
    com  = com/sum(this%m(1:this%natloc))

    if (present(vacuum)) then
       box(1)  = 2*maxval(POS(this, 1:this%natloc, 1)) + 2*vacuum(1)
       box(2)  = 2*maxval(POS(this, 1:this%natloc, 2)) + 2*vacuum(2)
       box(3)  = 2*maxval(POS(this, 1:this%natloc, 3)) + 2*vacuum(3)

       com  = com - box/2

       call set_cell(this, box, error=error)
       PASS_ERROR(error)
    endif

    if (present(cell)) then
       com  = com - (/ cell(1, 1)/2, cell(2, 2)/2, cell(3, 3)/2 /)
    endif

    do i = 1, this%natloc
#ifndef IMPLICIT_R
       POS3(this, i)  = POS3(this, i) - com
#endif
       PNC3(this, i)  = PNC3(this, i) - com
    enddo

  endsubroutine particles_center


  !>
  !! Notify the particles object of a change
  !!
  !! This function has to be called every time a change is made to the Particles object.
  !! For example, the neighbor list will only update if it detects a change to the
  !! Particles object.
  !!
  !! Internally, a counter is increased by one every time this function is called.
  !<
  subroutine particles_I_changed_positions(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    this%pos_rev  = this%pos_rev + 1

  endsubroutine particles_I_changed_positions


  !>
  !! Check if a change to the particles object has occured
  !!
  !! Internally, this compares the counter to a reference.
  !<
  logical function particles_have_positions_changed(this, last_rev)
    implicit none

    type(particles_t), intent(in)  :: this
    integer, intent(inout)         :: last_rev

    ! ---

    particles_have_positions_changed  = last_rev /= this%pos_rev
    last_rev                          = this%pos_rev

  endfunction particles_have_positions_changed


  !>
  !! Notify the particles object of a change
  !!
  !! This function has to be called every time a change is made to the Particles object.
  !! For example, the neighbor list will only update if it detects a change to the
  !! Particles object.
  !!
  !! Internally, a counter is increased by one every time this function is called.
  !<
  subroutine particles_I_changed_other(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    this%other_rev  = this%other_rev + 1

  endsubroutine particles_I_changed_other


  !>
  !! Check if a change to the particles object has occured
  !!
  !! Internally, this compares the counter to a reference.
  !<
  logical function particles_has_other_changed(this, last_rev)
    implicit none

    type(particles_t), intent(in)  :: this
    integer, intent(inout)         :: last_rev

    ! ---

    particles_has_other_changed  = last_rev /= this%other_rev
    last_rev                     = this%other_rev

  endfunction particles_has_other_changed



  !**********************************************************************
  ! Center around zero
  !**********************************************************************
  subroutine particles_translate(this, off)
    implicit none

    type(particles_t), intent(inout)   :: this
    real(DP), intent(in)               :: off(3)

    ! ---

#ifndef IMPLICIT_R
    POS(this, 1:this%natloc, 1)  = POS(this, 1:this%natloc, 1) + off(1)
    POS(this, 1:this%natloc, 2)  = POS(this, 1:this%natloc, 2) + off(2)
    POS(this, 1:this%natloc, 3)  = POS(this, 1:this%natloc, 3) + off(3)
#endif
    PNC(this, 1:this%natloc, 1)  = PNC(this, 1:this%natloc, 1) + off(1)
    PNC(this, 1:this%natloc, 2)  = PNC(this, 1:this%natloc, 2) + off(2)
    PNC(this, 1:this%natloc, 3)  = PNC(this, 1:this%natloc, 3) + off(3)

  endsubroutine particles_translate


  !**********************************************************************
  ! The volume of the current box
  !**********************************************************************
  real(DP) function particles_volume(p)
    implicit none

    type(particles_t), intent(in)  :: p

    ! ---

    real(DP)  :: vbox
    real(DP)  :: cross(3)
    integer   :: i

    ! ---

    cross(1)=p%Abox(2,2)*p%Abox(3,3)-p%Abox(3,2)*p%Abox(2,3)
    cross(2)=p%Abox(3,2)*p%Abox(1,3)-p%Abox(1,2)*p%Abox(3,3)
    cross(3)=p%Abox(1,2)*p%Abox(2,3)-p%Abox(2,2)*p%Abox(1,3)
    vbox=0d0
    do i=1,3
       vbox=vbox+p%Abox(i,1)*cross(i)
    enddo

    particles_volume  = vbox

  endfunction particles_volume


  !**********************************************************************
  ! Project r into a distance
  !**********************************************************************
  recursive function cyclic_in_bounds(p, r) result(cyc)
    implicit none

    type(particles_t), intent(in)  :: p

    real(DP), intent(in)  :: r(3)
    
    real(DP)              :: cyc(3)

    ! ---

    real(DP)  :: s(3)

    s    = matmul(p%Bbox, r)
    s    = s - nint(s)
    cyc  = matmul(p%Abox, s)

  endfunction cyclic_in_bounds


  !**********************************************************************
  ! Project r into the box
  !**********************************************************************
  function cyclic_in_cell(this, r) result(cyc)
    implicit none

    type(particles_t), intent(in)  :: this
    real(DP), intent(in)           :: r(3)
    
    real(DP)                       :: cyc(3)

    ! ---

    real(DP)  :: d(3), s(3)

    ! ---

    d    = this%shear_dx*floor(dot_product(this%Bbox(3, 1:3), r))
    s    = matmul(this%Bbox, r-d)
    s    = s - floor(s)
    cyc  = matmul(this%Abox, s)
    
  endfunction cyclic_in_cell


  !**********************************************************************
  ! Project r into the box
  !**********************************************************************
  function cyclic_in_cell_vec(this, r) result(cyc)
    implicit none

    type(particles_t), intent(in)  :: this
    real(DP), intent(in)           :: r(:, :)
    
    real(DP)                       :: cyc(3, size(r, 2))

    ! ---

    integer   :: i
    real(DP)  :: h(size(r,2)), d(3, size(r,2)), s(3, size(r, 2))

    ! ---

    h         = floor(matmul(this%Bbox(3, 1:3), r))
    forall(i=1:3)
      d(i,:)  = this%shear_dx(i)*h
    endforall
    s         = matmul(this%Bbox, r-d)
    s         = s - floor(s)
    cyc       = matmul(this%Abox, s)

  endfunction cyclic_in_cell_vec


  !**********************************************************************
  ! Project r into the box
  !**********************************************************************
  function cyclic_in_cellc(this, r, c) result(p)
    implicit none

    type(particles_t), intent(in)  :: this
    real(DP), intent(in)           :: r(3)
    integer, intent(in)            :: c

    real(DP)                       :: p
    
    ! ---

    real(DP)  :: d(3), s(3)

    ! ---

    d    = this%shear_dx*floor(dot_product(this%Bbox(3, 1:3), r))
    s    = matmul(this%Bbox, r-d)
    s    = s - floor(s)
    p    = dot_product(this%Abox(c, 1:3), s)

  endfunction cyclic_in_cellc


  !**********************************************************************
  ! Project r into the box
  !**********************************************************************
  function cyclic_in_cellc_vec(this, r, c) result(p)
    implicit none

    type(particles_t), intent(in)  :: this
    real(DP), intent(in)           :: r(:, :)
    integer, intent(in)            :: c

    real(DP)                       :: p(size(r, 2))
    
    ! ---

    integer   :: i
    real(DP)  :: h(size(r, 2)), d(3, size(r, 2)), s(3, size(r, 2)), cyc(3, size(r, 2))

    ! ---

    h           = floor(matmul(this%Bbox(3, 1:3), r))
    forall(i=1:3)
       d(i, :)  = this%shear_dx(i)*h
    endforall
    s           = matmul(this%Bbox, r-d)
    s           = s - floor(s)
    cyc         = matmul(this%Abox, s)
    p           = cyc(c, 1:size(r, 2))

  endfunction cyclic_in_cellc_vec


  !>
  !! Assign pointers to data
  !>
  subroutine particles_request_border(this, border)
    implicit none

    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: border

    ! ---

    this%border = max(border, this%border)

  endsubroutine particles_request_border


  !>
  !! Get effective box and reciprocal box, with consideration of Lees-Edwards
  !! boundary conditions.
  !<
  subroutine particles_get_true_cell(this, cell, rec_cell, error)
    implicit none

    type(particles_t),  intent(in)  :: this
    real(DP),           intent(out) :: cell(3,3)
    real(DP), optional, intent(out) :: rec_cell(3,3)
    integer,  optional, intent(out) :: error
    
    ! ---

    real(DP) :: A(3,3)
    integer  :: i, ipiv(3), info

    ! ---

    INIT_ERROR(error)

    if (any(this%shear_dx /= 0.0_DP)) then
       cell = this%Abox
       cell(1,3) = this%shear_dx(1)
       cell(2,3) = this%shear_dx(2)

       if (present(rec_cell)) then
          rec_cell = 0.0_DP
          do i = 1, 3
             rec_cell(i, i) = 1.0_DP
          enddo
          A  = cell
          call dgesv(3, 3, A, 3, ipiv, rec_cell, 3, info)

          if (info /= 0) then
             RAISE_ERROR("Failed to determine the reciprocal lattice. Cell = " // cell(:, 1) // ", " // cell(:, 2) // ", " // cell(:, 3), error)
          endif
       endif
    else
       cell = this%Abox
       if (present(rec_cell)) then
          rec_cell = this%Bbox
       endif
    endif

  endsubroutine particles_get_true_cell

endmodule particles
