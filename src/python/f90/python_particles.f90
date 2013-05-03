!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
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

     logical                :: periodic(3)
     logical                :: locally_periodic(3)

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

  public :: set_lees_edwards
  interface set_lees_edwards
     module procedure particles_set_lees_edwards
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
  subroutine particles_set_cell(this, Abox, pbc, scale_atoms, ierror)
    implicit none

    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: Abox(3, 3)
    logical, optional, intent(in)    :: pbc(3)
    logical, optional, intent(in)    :: scale_atoms
    integer, optional, intent(inout) :: ierror
    
    ! ---

    real(DP) :: A(3,3), fac(3, 3)
    integer  :: i, in, ipiv(3)

    ! --

!    call info("- particles_set_cell")

    if (present(pbc)) then
       this%periodic  = pbc
    endif

    if (.not. all(this%periodic)) then
       call require_orthorhombic_cell(this)
    endif

    if (present(scale_atoms) .and. scale_atoms) then
       fac = matmul(Abox, this%Bbox)
       !$omp  parallel do default(none) &
       !$omp& shared(this) firstprivate(fac)
       do i = 1, this%natloc
          PNC3(this, i) = matmul(fac, PNC3(this, i))
       enddo
    endif

!    this%Abox  = ( Abox + transpose(Abox) )/2
    this%Abox  = Abox

    this%Bbox  = 0.0_DP
    do i = 1, 3
       this%Bbox(i, i)  = 1.0_DP
    enddo

    this%cell_is_orthorhombic = &
         dot_product(this%Abox(1, :), this%Abox(2, :)) == 0.0_DP .and. &
         dot_product(this%Abox(2, :), this%Abox(3, :)) == 0.0_DP .and. &
         dot_product(this%Abox(3, :), this%Abox(1, :)) == 0.0_DP

!    if (.not. this%cell_is_orthorhombic) then
!       call info("     Cell is not orthorhombic.")
!    endif

!    call info("     " // this%Abox(1, :))
!    call info("     " // this%Abox(2, :))
!    call info("     " // this%Abox(3, :))

    A  = this%Abox
    call dgesv(3, 3, A, 3, ipiv, this%Bbox, 3, in)

    if (in /= 0) then
       RAISE_ERROR("Failed to determine the reciprocal lattice. Cell = " // this%Abox(:, 1) // ", " // this%Abox(:, 2) // ", " // this%Abox(:, 3), ierror)
    endif

    if (.not. all(this%periodic)) then
       call require_orthorhombic_cell(this)
    endif

    if (.not. this%cell_is_orthorhombic .and. this%orthorhombic_cell_is_required) then
       RAISE_ERROR("This cell is non-orthorhombic, however, an orthorhombic cell was required.", ierror)
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
  subroutine particles_set_cell_orthorhombic(this, cell, pbc, scale_atoms, ierror)
    implicit none

    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: cell(3)
    logical, optional, intent(in)    :: pbc(3)
    logical, optional, intent(in)    :: scale_atoms
    integer, optional, intent(inout) :: ierror

    ! ---

    real(DP)  :: cell3x3(3, 3)

    ! ---

    cell3x3        = 0.0_DP
    cell3x3(1, 1)  = cell(1)
    cell3x3(2, 2)  = cell(2)
    cell3x3(3, 3)  = cell(3)

    call particles_set_cell(this, cell3x3, pbc=pbc, scale_atoms=scale_atoms, ierror=ierror)

  endsubroutine particles_set_cell_orthorhombic


  !**********************************************************************
  ! Set Lees-Edwards boundary conditions
  !**********************************************************************
  subroutine particles_set_lees_edwards(this, dx, dv, ierror)
    implicit none

    type(particles_t), intent(inout)  :: this
    real(DP), intent(in)              :: dx(3)
    real(DP), intent(in), optional    :: dv(3)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer   :: i

    real(DP)  :: old_dx(3)

    ! ---

    call require_orthorhombic_cell(this, ierror)
    PASS_ERROR(ierror)

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

    this%periodic          = (/ .true., .true., .true. /)
    this%locally_periodic  = (/ .true., .true., .true. /)

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
    call add_real3( &
         this%data, &
         R_NON_CYC_STR, &
         F_TO_TRAJ + F_COMMUNICATE + F_COMM_GHOSTS, &
         "angstroms")

  endsubroutine particles_init

  
  !**********************************************************************
  ! Allocate particle information
  !**********************************************************************
  subroutine particles_init_from_particles(this, from, ierror)
    implicit none

    type(particles_t), intent(inout)  :: this
    type(particles_t), intent(in)     :: from
    integer, intent(inout), optional  :: ierror

    ! ---

    this%initialized                    = .true.
    this%orthorhombic_cell_is_required  = from%orthorhombic_cell_is_required
    
    this%periodic          = (/ .true., .true., .true. /)
    this%locally_periodic  = (/ .true., .true., .true. /)

    this%border = 0.0_DP

    call init(this%data, from%data)

    call set_cell(this, from%Abox, from%periodic, ierror=ierror)

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
  subroutine particles_allocate(this, nat, totnat, allow_def, ierror)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(in)               :: nat
    integer, intent(in), optional     :: totnat
    logical, intent(in), optional     :: allow_def
    integer, intent(inout), optional  :: ierror

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

    call set_cell(this, (/ 1.0_DP, 1.0_DP, 1.0_DP /), ierror=ierror)

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
          
          if (this%locally_periodic(j)) then
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
  subroutine particles_require_orthorhombic_cell(this, ierror)
    implicit none

    type(particles_t), intent(inout)  :: this
    integer, intent(inout), optional  :: ierror

    ! ---

    this%orthorhombic_cell_is_required  = .true.

    if (.not. this%cell_is_orthorhombic) then
       RAISE_ERROR("Orthorhombic cell is requested, however cell is already non-orthorhombic.", ierror)
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

    this%border = max(border, this%border)

  endsubroutine particles_request_border

endmodule particles
