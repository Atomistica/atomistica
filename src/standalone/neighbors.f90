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
!! Binning and neighbor list module
!<
module neighbors
  use, intrinsic :: iso_c_binding

  use supplib

  use particles

#ifdef _OPENMP
  use omp_lib
#endif

#ifdef _MP
  use communicator, ONLY: mod_communicator, communicate_ghosts, communicate_particles, request_border
#endif

  implicit none

  private

  integer, parameter  :: FIXED_VERLET_SHELL = 0
  integer, parameter  :: FIXED_CUTOFF = 1

  public :: NEIGHPTR_T
  integer, parameter :: NEIGHPTR_T = C_INTPTR_T

  public :: neighbors_t
  type neighbors_t

     !
     ! Particles object
     !

     type(particles_t), pointer  :: p      => NULL()       !< Associated particles object
     integer                     :: p_rev  = -1            !< Number of changes reference counter

     !
     ! Current status
     !

     logical                     :: initialized = .false.  !< Has this neighbor list been initialized? FIXME! Is this necessary?

     !
     ! Configuration
     !

     integer                     :: avgn                   !< Average number of neighbors
     integer                     :: mode                   !< Fixed cutoff or fixed Verlet shell

     real(DP)                    :: interaction_range      !< Maximum interaction range of potentials using this neighbor list
     real(DP)                    :: verlet_shell           !< Size of the Verlet shell

     real(DP)                    :: cutoff                 !< Cut-off for this neighbor list (i.e., interaction_range + verlet_shell)

     real(DP)                    :: requested_bin_size     !< Bin size that has been requested
     real(DP)                    :: bin_size               !< Actual bin size (i.e. such that is matches the box size)

     logical                     :: sort_before_update     !< Sort the particles before updating the neighbor list

     !
     ! Binning stuff
     !

     integer, allocatable        :: binning_seed(:, :, :)
     integer, allocatable        :: binning_last(:, :, :)

     integer, allocatable        :: next_particle(:)

     !
     ! Binning information
     ! If the system size is to small, n_cells gives the number of
     ! repetitions of the unit cell to consider.
     !

     real(DP)                    :: box_size(3)
     real(DP)                    :: Abox(3, 3)

     integer                     :: n_cells_tot
     integer                     :: n_cells(3)
     real(DP)                    :: cell_size(3, 3)
     real(DP)                    :: rec_cell_size(3, 3)

     !
     ! Neighbor stuff
     !

     integer                          :: n_d
     integer, allocatable             :: d(:, :)

     integer(NEIGHPTR_T), allocatable :: seed(:)         !< Seed for the neighbor list for the first set of particles
     integer(NEIGHPTR_T), allocatable :: last(:)         !< End of the neighbor list for the first set of particles

     integer                          :: neighbors_size  !< Size of the neighbor list
     integer, allocatable             :: neighbors(:)    !< Neighbor list for the second set of particles

     integer, allocatable             :: dc(:, :)        !< Which cell did the neighbor come from?

     !
     ! Other
     !

     integer                     :: it              !< Number of iteration since last construction of the neighbor list
     integer                     :: nupdate = 0     !< Number of total updates

     !
     ! Statistics
     !

     real(DP)                    :: avgnn           !< Average number of neighbors
     
     !
     ! Tag - this is used to attach the umbrella Python instance
     !
     
     type(C_PTR) :: tag

  endtype neighbors_t


  public :: init
  interface init
     module procedure neighbors_init, neighbors_copy
  endinterface

  public :: del
  interface del
     module procedure neighbors_del
  endinterface

  public :: set
  interface set
     module procedure neighbors_set
  endinterface

  public :: request_interaction_range
  interface request_interaction_range
     module procedure neighbors_request_interaction_range
  endinterface

  public :: update
  interface update
     module procedure neighbors_update
  endinterface

  public :: find_neighbor
  interface find_neighbor
     module procedure neighbors_find_neighbor
  endinterface

  public :: get_number_of_all_neighbors
  interface get_number_of_all_neighbors
     module procedure neighbors_get_number_of_all_neighbors
  endinterface

  public :: pack
  interface pack
     module procedure neighbors_pack_scalar
  endinterface

!--- Internal

  interface set_particles
     module procedure neighbors_set_particles
  endinterface

  interface binning_init
     module procedure neighbors_binning_init
  endinterface

  interface binning_del
     module procedure neighbors_binning_del
  endinterface

  interface binning_update
     module procedure neighbors_binning_update
  endinterface

contains

  !<
  !! Initialize the neighbor list
  !>
  subroutine neighbors_init(this, avgn, cutoff, verlet_shell, bin_size, sort, error)
    implicit none

    type(neighbors_t),  intent(inout)  :: this            !< Neighbor list object
    integer,  optional, intent(in)     :: avgn            !< Average number of neighbors
    real(DP), optional, intent(in)     :: cutoff          !< Cutoff
    real(DP), optional, intent(in)     :: verlet_shell    !< Verlet shell thickness
    real(DP), optional, intent(in)     :: bin_size        !< Binning size
    logical,  optional, intent(in)     :: sort            !< Sort list before update?
    integer,  optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    this%nupdate = 0

    call del(this)

    !
    ! Initialize
    !

    this%initialized         = .false.
    this%p                   => NULL()
    this%p_rev               = -1

    !
    ! Default values
    !

    this%avgn                = 100
    this%mode                = FIXED_VERLET_SHELL
    this%interaction_range   = 0.0_DP
    this%verlet_shell        = 0.0_DP
    this%cutoff              = 0.0_DP
    this%requested_bin_size  = -1.0_DP
    this%bin_size            = 0.0_DP
    this%sort_before_update  = .false.

    this%neighbors_size      = 0

    !
    ! Set values if present
    !

    call set(this, avgn, cutoff, verlet_shell, bin_size, sort, error=error)
    PASS_ERROR(error)

  endsubroutine neighbors_init


  !>
  !! Create a copy of a neighbor list object
  !<
  subroutine neighbors_copy(this, that, error)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(neighbors_t), intent(in)     :: that
    integer, optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    call init(this, that%avgn, that%cutoff, that%verlet_shell, that%bin_size, &
       that%sort_before_update, error=error)
    PASS_ERROR(error)

    this%seed       = that%seed
    this%last       = that%last
    this%neighbors  = that%neighbors
    this%dc         = that%dc

  endsubroutine neighbors_copy


  !>
  !! Destroy the neighbor list, i.e. free all memory
  !<
  subroutine neighbors_del(this)
    implicit none

    type(neighbors_t), intent(inout)  :: this

    ! ---

    if (this%nupdate > 0) then
       call prlog("- neighbors_del -")
       call prlog("     Average number of neighbors per atom = " // (this%avgnn/this%nupdate))
       call prlog
    endif

    this%nupdate             = 0
    this%avgnn               = 0.0_DP

    if (allocated(this%seed))       deallocate(this%seed)
    if (allocated(this%last))       deallocate(this%last)
    if (allocated(this%neighbors))  deallocate(this%neighbors)
    if (allocated(this%dc))         deallocate(this%dc)

    if (allocated(this%d))          deallocate(this%d)

    call binning_del(this)

    this%initialized  = .false.

  endsubroutine neighbors_del


  !>
  !! Set neighbor list parameters
  !<
  subroutine neighbors_set(this, avgn, cutoff, verlet_shell, bin_size, sort, error)
    implicit none

    type(neighbors_t),  intent(inout)  :: this
    integer,  optional, intent(in)     :: avgn
    real(DP), optional, intent(in)     :: cutoff
    real(DP), optional, intent(in)     :: verlet_shell
    real(DP), optional, intent(in)     :: bin_size
    logical,  optional, intent(in)     :: sort
    integer,  optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    call del(this)

    if (present(avgn)) then
       this%avgn  = avgn
    endif

    if (present(cutoff) .and. cutoff > 0.0_DP) then
       this%mode = FIXED_CUTOFF
       this%cutoff = cutoff
    endif

    if (present(verlet_shell) .and. verlet_shell > 0.0_DP) then
       if (present(cutoff) .and. cutoff > 0.0_DP) then
          RAISE_ERROR("Please specify either *cutoff* or *verlet_shell*, not both.", error)
       endif
       this%mode = FIXED_VERLET_SHELL
       this%verlet_shell = verlet_shell
    endif

    if (present(bin_size)) then
       this%requested_bin_size  = bin_size
    endif

    if (present(sort)) then
       this%sort_before_update  = sort
    endif

  endsubroutine neighbors_set


  !>
  !! Request an interaction range. This is called by the respective interatomic potentials to
  !! register the interaction range they require.
  !<
  subroutine neighbors_request_interaction_range(this, cutoff, Z1, Z2)
    implicit none

    type(neighbors_t), intent(inout) :: this
    real(DP),          intent(in)    :: cutoff
    integer, optional, intent(in)    :: Z1, Z2

    ! ---

    call del(this)

    if (cutoff > this%interaction_range) then

       this%interaction_range  = cutoff

       call prlog("- neighbors_request_interaction_range -")
       call prlog("     request                = " // cutoff)
       call prlog("     new interaction range  = " // this%interaction_range)
       call prlog

    endif

  endsubroutine neighbors_request_interaction_range


  !>
  !! Connect this neighbor list to a particles object.
  !<
  subroutine neighbors_set_particles(this, p)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(particles_t), target         :: p

    ! ---

    call del(this)

    this%p      => p
    this%p_rev  = p%pos_rev-1

  endsubroutine neighbors_set_particles


  !>
  !! Update the neighbor list, this will only happen if particles have been
  !! moved farther than the Verlet shell.
  !<
  subroutine neighbors_update(this, p, error)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(particles_t), target         :: p
    integer, intent(inout), optional  :: error

    ! ---
  
    INIT_ERROR(error)

    !call timer_start('neighbors_update')

    if (.not. associated(this%p, p)) then
       call set_particles(this, p)
    endif

    if (.not. this%initialized .or. have_positions_changed(p, this%p_rev)) then
       call refresh_neighbor_list(this, p, error)
       PASS_ERROR(error)
    endif

    !call timer_stop('neighbors_update')

  endsubroutine neighbors_update


  subroutine refresh_neighbor_list(this, p, error)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: error

    ! ---

    logical  :: update_now

    ! ---
 
    INIT_ERROR(error)

    if (.not. this%initialized) then
       call prlog("- neighbors_update -")
       call prlog("     Initializing neighbor list.")

       if (this%sort_before_update) then
          call prlog("     Sorting particles before neighbor list update.")
       endif

       if (this%mode == FIXED_VERLET_SHELL) then
          this%cutoff  = this%interaction_range + this%verlet_shell
       else if (this%mode == FIXED_CUTOFF) then
          if (this%cutoff < this%interaction_range) then
             RAISE_ERROR("Cutoff " // this%cutoff // " smaller than the current interaction range " // this%interaction_range // ". Please increase.", error)
          endif
          this%verlet_shell = this%cutoff - this%interaction_range
       else
          RAISE_ERROR("Internal error: mode = " // this%mode, error)
       endif

       call prlog("     interaction_range  = "//this%interaction_range)
       call prlog("     verlet_shell       = "//this%verlet_shell)
       call prlog("     cutoff             = "//this%cutoff)

       if (this%cutoff <= 0.0_DP) then
          RAISE_ERROR("Cutoff needs to be larger than zero.", error)
       endif

       call prlog("     avgn               = "//this%avgn)

       this%neighbors_size  = p%maxnatloc * this%avgn

       allocate(this%seed(p%maxnatloc+1))
       allocate(this%last(p%maxnatloc+1))
       allocate(this%neighbors(this%neighbors_size))
       allocate(this%dc(3, this%neighbors_size))

       call log_memory_start("neighbors_update")

       call log_memory_estimate(this%seed)
       call log_memory_estimate(this%last)
       call log_memory_estimate(this%neighbors)
       call log_memory_estimate(this%dc)

       call log_memory_stop("neighbors_update")

       if (this%requested_bin_size > 0.0_DP) then
          this%bin_size  = this%requested_bin_size
       else
          this%bin_size  = this%cutoff
       endif

       p%accum_max_dr  = this%verlet_shell + 1.0_DP

       this%it              = 0
       this%initialized     = .true.

       call prlog

#ifdef _MP
       call request_border(mod_communicator, p, this%interaction_range, &
            verlet_shell = this%verlet_shell, &
            error = error)
       PASS_ERROR(error)
#endif

       call binning_init(this, p, error)
       PASS_ERROR(error)
    else

       if (any(this%Abox /= p%Abox)) then
#ifdef _MP
          call request_border(mod_communicator, p, this%interaction_range, &
               verlet_shell = this%verlet_shell, &
               error = error)
          PASS_ERROR(error)
#endif

          call binning_init(this, p, error)
          PASS_ERROR(error)
       endif

    endif

    !
    ! Update the neighbor list
    !

    call pnc2pos(p)

    ! Factor of 2* is because one particle can move right
    ! while the other particles moves opposite, hence the
    ! distance changes by 2*accum_max_dr.
    update_now  = 2*p%accum_max_dr >= this%verlet_shell

#ifdef _MP
    update_now = any(mod_communicator%mpi, update_now, error=error)
    PASS_ERROR(error)
#endif

    if (update_now) then

       if (this%sort_before_update) then
          call timer_start("sort_particles")
          call sort(p, p%sort_index)
          call timer_stop("sort_particles")
       endif

       call inbox(p)

#ifdef _MP
       DEBUG_WRITE("- communicate_particles -")
       call communicate_particles(mod_communicator, p)
       DEBUG_WRITE("- communicate_ghosts -")
       call communicate_ghosts(mod_communicator, p, .true.)
#endif

       this%it         = 0
       p%accum_max_dr  = 1d-6

       call binning_update(this, p, error)
       PASS_ERROR(error)
       call fill_neighbor_list(this, p, error)
       PASS_ERROR(error)

    else

#ifdef _MP
       DEBUG_WRITE("- communicate_ghosts -")
       call communicate_ghosts(mod_communicator, p, .false.)
#endif

       this%it  = this%it + 1

    endif
    
  endsubroutine refresh_neighbor_list


  !>
  !! Find all neighbors for these particles using binning. Do not call, used internally.
  !<
  recursive subroutine fill_neighbor_list(this, p, error)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(particles_t), intent(in)     :: p
    integer, intent(inout), optional  :: error

    ! ---

    real(DP)  :: Abox(3, 3), shear_dx(3)
    logical   :: locally_pbc(3)

    integer   :: i, j, k, x, nn
    integer   :: cell(3), cell2(3), cur_cell(3)
    integer   :: cur

    real(DP)  :: delta_r(3), abs_delta_r_sq, off(3)
    
    real(DP)  :: cutoff_sq

    integer   :: c(3)
    
    integer   :: chunk_len

    logical   :: lebc, any_c_not_zero

    integer   :: error_loc

#ifdef _OPENMP
    integer   :: chunk_start
#endif

    ! ---

    INIT_ERROR(error)

    call timer_start("fill_neighbor_list")

    Abox        = p%Abox
    locally_pbc = p%locally_pbc
    shear_dx    = p%shear_dx

#ifdef _OPENMP
    chunk_len = size(this%neighbors)/omp_get_max_threads()
#else
    chunk_len = size(this%neighbors)
#endif

    cutoff_sq = this%cutoff**2

    lebc  = .false.
    if (any(shear_dx /= 0.0_DP)) then
       lebc = .true.
    endif

    error_loc = ERROR_NONE
    nn = 0

#ifdef _OPENMP
    !$omp  parallel default(none) &
    !$omp& private(abs_delta_r_sq, cell, chunk_start, c, cell2, cur, cur_cell, delta_r, i, j, off, x, any_c_not_zero) &
    !$omp& firstprivate(chunk_len, cutoff_sq, lebc, Abox, shear_dx, locally_pbc) &
    !$omp& shared(this, p) &
    !$omp& reduction(+:error_loc) reduction(+:nn)

    chunk_start = 1 + omp_get_thread_num()*chunk_len
    cur = chunk_start
    
!    !omp critical
!    write (*, *)  cur, omp_get_thread_num()
!    !omp end critical
#else
    cur = 1
#endif

    !$omp do
    i_loop: do i = 1, p%nat
       ! Compute the 3-index of the current cell
       cell = floor( matmul(this%rec_cell_size, PNC3(p, i) - p%lower_with_border) + 1 )

       this%seed(i) = cur

       ! Loop over all (precomputed) cell distances in x-, y- and z-direction
       xyz_loop: do x = 1, this%n_d
          cur_cell  = cell + this%d(1:3, x)
          c         = 0

          ! Determine whether the neighboring cell is outside of the simulation
          ! domain and needs to be remapped into it. The variable c counts the
          ! distance to that cell in number of simulation cells. This allows
          ! the use of very small boxes, e.g. one atom simulations.
          do k = 1, 3
             if (locally_pbc(k)) then
                do while (cur_cell(k) < 1)
                   cur_cell(k) = cur_cell(k)+this%n_cells(k)
                   c(k)        = c(k)-1
                enddo
                do while (cur_cell(k) > this%n_cells(k))
                   cur_cell(k) = cur_cell(k)-this%n_cells(k)
                   c(k)        = c(k)+1
                enddo
             endif
          enddo

          ! Do we have Lees-Edwards boundary conditions and a jump across the
          ! simulation domain boundaries in z-direction? Then we need to
          ! re-evaluate the cell index considering the boundary conditions
          if (lebc .and. c(3) /= 0) then
             ! Compute 3-index of the cell for the (Lees-Edwards) shifted grid:
             cell2     = floor( matmul(this%rec_cell_size, PNC3(p, i) - c(3)*shear_dx - p%lower_with_border) + 1 )
             cur_cell  = cell2 + this%d(1:3, x)

             c         = 0

             do k = 1, 3
                if (locally_pbc(k)) then
                   do while (cur_cell(k) < 1)
                      cur_cell(k) = cur_cell(k)+this%n_cells(k)
                      c(k)        = c(k)-1
                   enddo
                   do while (cur_cell(k) > this%n_cells(k))
                      cur_cell(k) = cur_cell(k)-this%n_cells(k)
                      c(k)        = c(k)+1
                   enddo
                endif
             enddo
          endif

          cell_exists: if (.not. (any(cur_cell < 1) .or. any(cur_cell > this%n_cells)) .and. error_loc == ERROR_NONE) then
             any_c_not_zero  = any(c /= 0)
             off             = matmul(Abox, c) + c(3)*shear_dx
             j               = this%binning_seed(cur_cell(1), cur_cell(2), cur_cell(3))

             do while (j /= -1)
                if (i /= j .or. any_c_not_zero) then
                   delta_r = PNC3(p, i) - PNC3(p, j) - off

                   abs_delta_r_sq = dot_product(delta_r, delta_r)

                   if (abs_delta_r_sq < cutoff_sq) then
#ifdef _OPENMP
                      if (cur - chunk_start >= chunk_len) then
                         RAISE_DELAYED_ERROR("Neighbor list overflow. Current neighbor list position is " // cur // " while the size of this chunk runs from " // chunk_start // " to " // chunk_len // ".", error_loc)
#else
                      if (cur >= chunk_len) then
                         RAISE_ERROR("Neighbor list overflow. Current neighbor list position is " // cur // " while the size of this chunk runs from 1 to " // chunk_len // ".", error)
#endif
                      else
                         this%neighbors(cur)  = j
                         VEC3(this%dc, cur)   = -c

                         cur                  = cur + 1
                         nn                   = nn + 1
                      endif
                   endif
                endif

                j  = this%next_particle(j)
             enddo

          endif cell_exists

       enddo xyz_loop

       this%last(i)  = cur-1

       this%neighbors(cur)  = 0
       cur  = cur+1
    enddo i_loop
    !$omp end do
    !$omp end parallel

    INVOKE_DELAYED_ERROR(error_loc, error)

    this%seed(p%nat+1) = cur

    this%nupdate = this%nupdate + 1
    this%avgnn   = this%avgnn + real(nn, DP)/p%nat

    call timer_stop("fill_neighbor_list")

  endsubroutine fill_neighbor_list


  !
  ! Binning
  !

  !>
  !! Initialize global cell-subdivision, i.e., estimate a cell size
  !! from the given average density
  !<
  subroutine neighbors_binning_init(this, p, error)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(particles_t), intent(in)     :: p
    integer, intent(inout), optional  :: error

    ! ---

    real(DP)  :: cutoff_sq, nx(3), ny(3), nz(3), cv
    integer   :: x, y, z, dx, dy, dz, dy2, dz2

    ! ---

    INIT_ERROR(error)

    this%Abox         = p%Abox

    if (p%cell_is_orthorhombic) then
       this%box_size  = p%upper_with_border - p%lower_with_border
    else
       forall(x=1:3)
          this%box_size(x)  = sqrt(dot_product(this%Abox(:, x), this%Abox(:, x)))
       endforall
    endif

    this%n_cells  = int(this%box_size / this%bin_size)

    ! Enforce three cells minimum
    where (this%n_cells < 3)
       this%n_cells  = 3
    endwhere
    this%n_cells_tot  = this%n_cells(1)*this%n_cells(2)*this%n_cells(3)

    cutoff_sq         = this%cutoff**2

    !
    ! Otherwise, enable cell subdivision
    !

    if (p%cell_is_orthorhombic) then
       this%cell_size      = diagonal_matrix(this%box_size / this%n_cells)
       this%rec_cell_size  = diagonal_matrix(this%n_cells / this%box_size)
    else
       forall(x=1:3)
          this%cell_size(:, x)      = this%Abox(:, x) / this%n_cells(x)
          this%rec_cell_size(x, :)  = p%Bbox(x, :) * this%n_cells(x)
       endforall

    endif

    if (allocated(this%binning_seed)) then
       ! Number of cells changed?
       if (any(shape(this%binning_seed) /= this%n_cells)) then
          deallocate(this%binning_seed)
          deallocate(this%binning_last)
          allocate(this%binning_seed(this%n_cells(1), this%n_cells(2), this%n_cells(3)))
          allocate(this%binning_last(this%n_cells(1), this%n_cells(2), this%n_cells(3)))
       endif
    else
       call prlog("- neighbors_binning_init -")
       call prlog("     Binning enabled.")
       call prlog("     cutoff     = "//this%bin_size)
       call prlog("     box_size   = ( "//this%box_size//" )")
       call prlog("     n_cells    = ( "//this%n_cells//" )")
       !call prlog("     cell_size  = ( "//this%cell_size//" )")
       call prlog
       
       allocate(this%binning_seed(this%n_cells(1), this%n_cells(2), this%n_cells(3)))
       allocate(this%binning_last(this%n_cells(1), this%n_cells(2), this%n_cells(3)))
       allocate(this%next_particle(p%maxnatloc))
    endif

    !
    ! Create cell list for neighbor search
    !

    ! Compute the surface normal vectors
    nx  = cross_product(this%cell_size(:, 2), this%cell_size(:, 3))
    ny  = cross_product(this%cell_size(:, 3), this%cell_size(:, 1))
    nz  = cross_product(this%cell_size(:, 1), this%cell_size(:, 2))

    ! The cell volume
    cv  = dot_product(this%cell_size(:, 1), nx)

    ! Adjust the length of the surface normal vectors such that they point to
    ! the opposite surface
    nx  = cv * nx / dot_product(nx, nx)
    ny  = cv * ny / dot_product(ny, ny)
    nz  = cv * nz / dot_product(nz, nz)

    ! Now dx, dy, dz needs to be adjusted such that a sphere of size cutoff
    ! fits into the box spanned by (dx,dy,dz)*(nx,ny,nz)
    dx  = int(this%cutoff/sqrt(dot_product(nx, nx))) + 1
    dy  = int(this%cutoff/sqrt(dot_product(ny, ny))) + 1
    dz  = int(this%cutoff/sqrt(dot_product(nz, nz))) + 1

    if (allocated(this%d) .and. size(this%d, 2) < (2*dx+1)*(2*dy+1)*(2*dz+1)) then
       deallocate(this%d)
    endif
       
    if (.not. allocated(this%d)) then
       allocate(this%d(3, (2*dx+1)*(2*dy+1)*(2*dz+1)))
    endif

    this%n_d  = 0

    if (p%cell_is_orthorhombic) then

       x_loop2: do x = -dx, dx
          if (x == 0) then
             dy2  = int(this%cutoff*this%rec_cell_size(2, 2))+1
          else
             dy2  = int(sqrt(cutoff_sq - ((abs(x)-1)*this%cell_size(1, 1))**2)*this%rec_cell_size(2, 2))+1
          endif
          dy2 = min(dy, dy2)

          y_loop2: do y = -dy2, dy2
             if (x == 0) then
                if (y == 0) then
                   dz2  = int(this%cutoff*this%rec_cell_size(3, 3))+1
                else
                   dz2  = int(sqrt(cutoff_sq - ((abs(y)-1)*this%cell_size(2, 2))**2)*this%rec_cell_size(3, 3))+1
                endif
             else
                if (y == 0) then
                   dz2  = int(sqrt(cutoff_sq - ((abs(x)-1)*this%cell_size(1, 1))**2)*this%rec_cell_size(3, 3))+1
                else
                   dz2  = int(sqrt(cutoff_sq - ((abs(x)-1)*this%cell_size(1, 1))**2 - ((abs(y)-1)*this%cell_size(2, 2))**2)*this%rec_cell_size(3, 3))+1
                endif
             endif
             dz2 = min(dz, dz2)

             z_loop2: do z = -dz2, dz2

                this%n_d = this%n_d+1
                this%d(1:3, this%n_d) = (/ x, y, z /)

             enddo z_loop2

          enddo y_loop2
       enddo x_loop2

    else

       x2_loop2: do x = -dx, dx
          y2_loop2: do y = -dy, dy
             z2_loop2: do z = -dz, dz

                this%n_d = this%n_d+1
                this%d(1:3, this%n_d) = (/ x, y, z /)

             enddo z2_loop2
          enddo y2_loop2
       enddo x2_loop2
          
    endif

    call prlog

  endsubroutine neighbors_binning_init


  !>
  !! Destructor
  !!
  !! Delete the binning structure
  !<
  subroutine neighbors_binning_del(this)
    implicit none

    type(neighbors_t), intent(inout)  :: this
 
    ! ---

    if (allocated(this%binning_seed)) then
       deallocate(this%binning_seed)
    endif
    if (allocated(this%binning_last)) then
       deallocate(this%binning_last)
    endif
    if (allocated(this%next_particle)) then
       deallocate(this%next_particle)
    endif

  endsubroutine neighbors_binning_del


  !>
  !! Bin the particles into the corresponding binning structure for subsequent neighbors search.
  !!
  !! Bin the particles into the corresponding binning structure for subsequent neighbors search.
  !<
  subroutine neighbors_binning_update(this, p, error)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: error

    ! ---

    integer   :: i, j, cell(3)

    ! ---

    INIT_ERROR(error)

    this%next_particle  = -1

    this%binning_seed   = -1
    this%binning_last   = -1

    do i = 1, p%nat
       cell  = floor(matmul(this%rec_cell_size, PNC3(p, i)-p%lower_with_border))+1

       ! Fix for the "Particle outside simulation domain" problem
       do j = 1, 3
          if (p%locally_pbc(j)) then
             if (cell(j) < 1) then
                cell(j) = cell(j) + this%n_cells(j)
                PNC3(p, i) = PNC3(p, i) + p%Abox(1:3, j)
             else if (cell(j) > this%n_cells(j)) then
                cell(j) = cell(j) - this%n_cells(j)
                PNC3(p, i) = PNC3(p, i) - p%Abox(1:3, j)
             endif
          endif
       enddo

       if (any(cell < 1) .or. any(cell > this%n_cells)) then
          call particles_dump_info(p, i, cell)
          RAISE_ERROR("Particle outside simulation domain.", error)
       endif

       if (this%binning_seed(cell(1), cell(2), cell(3)) == -1) then
          this%binning_seed(cell(1), cell(2), cell(3))  = i
          this%binning_last(cell(1), cell(2), cell(3))  = i
       else
          this%next_particle(this%binning_last(cell(1), cell(2), cell(3)))  = i
          this%binning_last(cell(1), cell(2), cell(3))                      = i
       endif

       p%sort_index(i)  = xyz2index(cell(1), cell(2), cell(3), this%n_cells)
    enddo

  endsubroutine neighbors_binning_update


  !>
  !! Search for the pair \param i - \param j and return the neighbor index
  !!
  !! This method searches for the pair \param i - \param j and return the neighbor index.
  !! Returned will be both, the \param i - \param j and \param j - \param i index
  !! in the parameters \param n1 and \param n2.
  !<
  subroutine neighbors_find_neighbor(this, i, j, n1, n2)
    implicit none

    type(neighbors_t), intent(inout)  :: this
    integer, intent(in)               :: i
    integer, intent(in)               :: j
    integer, intent(out)              :: n1
    integer, intent(out)              :: n2

    ! ---

    integer  :: n

    ! ---

    n1  = -1
    n2  = -1

    do n = this%seed(i), this%last(i)
       if (this%neighbors(n) == j) then
          n1  = n
       endif
    enddo

    do n = this%seed(j), this%last(j)
       if (this%neighbors(n) == i) then
          n2  = n
       endif
    enddo

  endsubroutine neighbors_find_neighbor


  !>
  !! Return total number of neighbors
  !<
  function neighbors_get_number_of_all_neighbors(this) result(s)
    use, intrinsic :: iso_c_binding

    implicit none

    type(neighbors_t), intent(in) :: this
    integer                       :: s

    ! ---

    integer :: i

    ! ---

    s = 0
    do i = 1, this%p%nat
       s = s + this%last(i)-this%seed(i)+1
    enddo

  endfunction neighbors_get_number_of_all_neighbors


  !>
  !! Bring a list of scalar per bond information into order
  !<
  subroutine neighbors_pack_scalar(this, r1, r2)
    implicit none

    type(neighbors_t), intent(in)  :: this
    real(DP),          intent(in)  :: r1(*)
    real(DP),          intent(out) :: r2(*)

    ! ---

    integer :: i, ni, j

    ! ---

    j = 0
    do i = 1, this%p%nat
       do ni = this%seed(i), this%last(i)
          j = j + 1
          r2(j) = r1(ni)
       enddo
    enddo

  endsubroutine neighbors_pack_scalar

endmodule neighbors
