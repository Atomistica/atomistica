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
!! of this module.
!<
module particles
  use, intrinsic :: iso_c_binding

  use libAtoms_module

  use logging
  use misc

  use data

  implicit none

  private

  !>
  !! Highest element number stored in the periodic table module
  !<
  integer, parameter :: MAX_Z = ubound(ElementName, 1)

  public :: Z_STR, EL_STR, R_NON_CYC_STR
  public :: SHEAR_DX_STR, CELL_STR, MAX_Z

  character(MAX_NAME_STR), parameter  :: Z_STR            = "Z"
  character(MAX_NAME_STR), parameter  :: Z_ALIAS_STR      = "atom_types"
  character(MAX_NAME_STR), parameter  :: EL_STR           = "internal_element_number"
  character(MAX_NAME_STR), parameter  :: R_NON_CYC_STR    = "coordinates"                  ! ... are allowed to leave the cell between neighbor list updates

  character(MAX_NAME_STR), parameter  :: SHEAR_DX_STR     = "shear_dx"
  character(MAX_NAME_STR), parameter  :: CELL_STR         = "cell"
  character(MAX_NAME_STR), parameter  :: PBC_STR          = "pbc"

  character(MAX_NAME_STR), parameter  :: V_STR            = "velocities"

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

     integer                :: pos_rev = 0          !> Have the positions been changed?
     integer                :: cell_rev = 0          !> Has the cell been changed?
     integer                :: other_rev = 0          !> Has anything else been changed?

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

     integer, pointer       :: pbc(:)

     !
     ! Lees-Edwards boundary conditions
     !

     real(DP), pointer      :: shear_dx(:)
     real(DP)               :: shear_dv(3)

     !
     ! Accumulated distance moved (not actually used in the Python interface)
     !

     real(DP)               :: accum_max_dr

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
     ! All particle data is managed by the *data* field. The other fields are
     ! pointers to the entries of data.
     !

     type(data_t), pointer  :: data

     integer, pointer       :: Z(:)                 ! element number
     integer, pointer       :: el(:)                ! element number

     ! These positions are always local and may be outside the global box.
     ! These differ on for different processes.
     real(DP), pointer      :: r_non_cyc(:, :)      ! displacement from last binning

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

     type(C_PTR)            :: tag

  endtype particles_t

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

  public :: I_changed_positions
  interface I_changed_positions
     module procedure particles_I_changed_positions
  endinterface

  public :: have_positions_changed
  interface have_positions_changed
     module procedure particles_have_positions_changed
  endinterface

  public :: I_changed_cell
  interface I_changed_cell
     module procedure particles_I_changed_cell
  endinterface

  public :: has_cell_changed
  interface has_cell_changed
     module procedure particles_has_cell_changed
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

  public :: require_orthorhombic_cell
  interface require_orthorhombic_cell
     module procedure particles_require_orthorhombic_cell
  endinterface

  public :: set_cell
  interface set_cell
     module procedure particles_set_cell, particles_set_cell_orthorhombic
  endinterface

  public :: get_true_cell
  interface get_true_cell
    module procedure particles_get_true_cell
  endinterface

  public :: set_lees_edwards
  interface set_lees_edwards
     module procedure particles_set_lees_edwards
  endinterface

  public :: volume
  interface volume
     module procedure particles_volume
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

    real(DP), parameter :: TOL = 1e-9

    ! ---

    real(DP) :: A(3,3), fac(3, 3)
    integer  :: i, in, ipiv(3)

    ! --

    if (present(pbc)) then
       where (pbc)
          this%pbc = 1
       elsewhere
          this%pbc = 0
       endwhere
    endif

    if (present(scale_atoms) .and. scale_atoms) then
       fac = matmul(Abox, this%Bbox)
       !$omp  parallel do default(none) &
       !$omp& shared(this) firstprivate(fac)
       do i = 1, this%natloc
          PNC3(this, i) = matmul(fac, PNC3(this, i))
       enddo
    endif

    this%Abox  = Abox

    this%Bbox  = 0.0_DP
    do i = 1, 3
       this%Bbox(i, i)  = 1.0_DP
    enddo

    this%cell_is_orthorhombic = &
         abs(dot_product(this%Abox(1, :), this%Abox(2, :))) < TOL .and. &
         abs(dot_product(this%Abox(2, :), this%Abox(3, :))) < TOL .and. &
         abs(dot_product(this%Abox(3, :), this%Abox(1, :))) < TOL

    if (.not. all(this%pbc /= 0)) then
       call require_orthorhombic_cell(this, error)
       PASS_ERROR(error)
    endif

    if (.not. this%cell_is_orthorhombic .and. this%orthorhombic_cell_is_required) then
       RAISE_ERROR("This cell is non-orthorhombic, however, an orthorhombic cell was required. Cell = " // this%Abox(:, 1) // ", " // this%Abox(:, 2) // ", " // this%Abox(:, 3), error)
    endif

    if (this%cell_is_orthorhombic) then
       A = 0.0_DP
       A(1,1) = this%Abox(1,1)
       A(2,2) = this%Abox(2,2)
       A(3,3) = this%Abox(3,3)
       this%Abox = A
    endif

    A  = this%Abox
    call dgesv(3, 3, A, 3, ipiv, this%Bbox, 3, in)

    if (in /= 0) then
       RAISE_ERROR("Failed to determine the reciprocal lattice. Cell = " // this%Abox(:, 1) // ", " // this%Abox(:, 2) // ", " // this%Abox(:, 3), error)
    endif

    this%lower  = (/ 0.0, 0.0, 0.0 /)
    this%upper  = (/ this%Abox(1, 1), this%Abox(2, 2), this%Abox(3, 3) /)

    this%lower_with_border = this%lower
    this%upper_with_border = this%upper

    call I_changed_cell(this)

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
    integer  :: ipiv(3), info

    ! ---

    INIT_ERROR(error)

    if (any(this%shear_dx /= 0.0_DP)) then
       cell = this%Abox
       cell(3,1) = this%shear_dx(1)
       cell(3,2) = this%shear_dx(2)

       if (present(rec_cell)) then
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

    this%border = 0.0_DP

    allocate(this%data)
    call init(this%data)

    call add_real3x3_attr( &
         this%data, &
         CELL_STR)
    call add_integer3_attr( &
         this%data, &
         PBC_STR)
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
    call add_real3( &
         this%data, &
         R_NON_CYC_STR, &
         F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS, &
         "angstroms")

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
    
    this%pbc = (/ 1, 1, 1 /)

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

    call allocate(this%data, nat, allow_def)

    this%nat            = nat
    this%natloc         = nat
    this%maxnatloc      = nat
    this%totnat         = nat

    if (present(totnat)) then
       this%totnat  = totnat
    endif

    call particles_assign_ptrs(this)

    this%Z             = 1

    call set_cell(this, (/ 1.0_DP, 1.0_DP, 1.0_DP /), &
         (/ .true., .true., .true. /), error=error)

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

  endsubroutine particles_del


  !**********************************************************************
  ! Assign shortcuts (i.e. r, v) to field in the *data* object
  !**********************************************************************
  subroutine particles_assign_ptrs(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    call attr_by_name(this%data, CELL_STR, this%Abox)
    call attr_by_name(this%data, PBC_STR, this%pbc)
    call attr_by_name(this%data, SHEAR_DX_STR, this%shear_dx)

    call ptr_by_name(this%data, Z_STR, this%Z)
    call ptr_by_name(this%data, EL_STR, this%el)
    call ptr_by_name(this%data, R_NON_CYC_STR, this%r_non_cyc)

  endsubroutine particles_assign_ptrs


  !**********************************************************************
  ! Compute statistics
  !**********************************************************************
  subroutine particles_update_elements(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    integer  :: i

    ! ---

    write (ilog, *)  "- particles_update_elements -"

    this%nZ  = 0

    do i = 1, this%natloc
       this%nZ(this%Z(i)) = this%nZ(this%Z(i))+1
    enddo

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

    write (ilog, '(4X,I2,A)')  this%nel, " elements found."
    do i = 1, this%nel
       write (ilog, '(4X,I2,A,A2,A,I6,A)')  i, " = ", ElementName(this%el2Z(i)), " (", this%nZ(this%el2Z(i)), " atoms found)"
    enddo

    write (ilog, *)

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
          
          if (this%pbc(j) /= 0) then
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

  endsubroutine particles_inbox


  !>
  !! Return the total cell volume
  !<
  real(DP) function particles_volume(p)
    implicit none

    type(particles_t), intent(in)  :: p

    ! ---

    real(DP)  :: vbox
    real(DP)  :: cross(3)
    integer   :: i

    ! ---

    cross(1) = p%Abox(2,2)*p%Abox(3,3)-p%Abox(3,2)*p%Abox(2,3)
    cross(2) = p%Abox(3,2)*p%Abox(1,3)-p%Abox(1,2)*p%Abox(3,3)
    cross(3) = p%Abox(1,2)*p%Abox(2,3)-p%Abox(2,2)*p%Abox(1,3)
    vbox = 0.0_DP
    do i = 1, 3
       vbox = vbox + p%Abox(i,1)*cross(i)
    enddo

    particles_volume = vbox

  endfunction particles_volume


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
    write (ilog, '(A, i15)')     "Z                  = ", this%Z(i)
    write (ilog, '(A)')          "---"
    write (ilog, '(A, 3ES15.8)') "r                  = ", POS3(this, i)
    write (ilog, '(A, 3ES15.8)') "r_non_cyc          = ", PNC3(this, i)
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
  ! Call to require an orthorhombic cell
  !**********************************************************************
  subroutine particles_require_orthorhombic_cell(this, error)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(inout), optional  :: error

    ! ---

    this%orthorhombic_cell_is_required  = .true.

    if (.not. this%cell_is_orthorhombic) then
       RAISE_ERROR("Orthorhombic cell is requested, however cell is already non-orthorhombic. Cell = " // this%Abox(:, 1) // ", " // this%Abox(:, 2) // ", " // this%Abox(:, 3), error)
    endif

  endsubroutine particles_require_orthorhombic_cell


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
  subroutine particles_I_changed_cell(this)
    implicit none

    type(particles_t), intent(inout)  :: this

    ! ---

    this%cell_rev  = this%cell_rev + 1

  endsubroutine particles_I_changed_cell


  !>
  !! Check if a change to the particles object has occured
  !!
  !! Internally, this compares the counter to a reference.
  !<
  logical function particles_has_cell_changed(this, last_rev)
    implicit none

    type(particles_t), intent(in)  :: this
    integer, intent(inout)         :: last_rev

    ! ---

    particles_has_cell_changed  = last_rev /= this%cell_rev
    last_rev                    = this%cell_rev

  endfunction particles_has_cell_changed


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
    last_rev                    = this%other_rev

  endfunction particles_has_other_changed


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
    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: border

    ! ---

    if (border > this%border) then

       call prlog("- particles_request_border -")
       call prlog("     old border  = " // this%border)
       call prlog("     request     = " // border)
       call prlog

       this%border = border

    endif

  endsubroutine particles_request_border

endmodule particles
