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

!**********************************************************************
! 3D domain decomposition
!**********************************************************************

#include "macros.inc"

#ifdef _MP

module communicator
  use supplib

  use particles

  implicit none

  private

  public :: communicator_t
  type communicator_t

     !
     ! Decomposition type
     !

     integer   :: decomposition(3)    = (/ 2, 2, 2 /)

     real(DP)  :: requested_border    = 0.0_DP
     real(DP)  :: border              = 0.0_DP
     real(DP)  :: verlet_shell        = 0.0_DP

     logical   :: communicate_forces  = .false.

     !
     ! MPI-Communicator
     !

     type(MPI_context)  :: mpi

     !
     ! Periodicity
     !

     logical    :: pbc(3)

     !
     ! Neighboring processeses in x-, y- and z- direction
     !

     integer    :: r(3), l(3)

     !
     ! Offsets
     !

     real(DP)   :: off_r(3), off_l(3)

     !
     ! Sizes
     !

     integer    :: n_particle_data
     integer    :: n_ghost_data
     integer    :: n_force_data

     !
     ! Lists containing particle information
     !

     real(DP), pointer  :: send_l(:)  => NULL()
     real(DP), pointer  :: send_r(:)  => NULL()
     real(DP), pointer  :: recv_l(:)  => NULL()
     real(DP), pointer  :: recv_r(:)  => NULL()

     !
     ! Lists containing pointers to ghost particles
     !

     integer            :: n_ghosts_r(3), n_ghosts_l(3)

     integer, pointer   :: ghosts_r(:)  => NULL()
     integer, pointer   :: ghosts_l(:)  => NULL()

     !
     ! Statistics
     !

     integer            :: n_send_p_tot
     integer            :: n_recv_p_tot
     integer            :: n_send_g_tot
     integer            :: n_recv_g_tot
     integer            :: nit_p
     integer            :: nit_g
     
  endtype communicator_t


  !
  ! The global parallelization module
  !

  public :: mod_communicator
  type(communicator_t), target, save :: mod_communicator

  public :: init
  interface init
     module procedure communicator_init
  endinterface

  public :: del
  interface del
     module procedure communicator_del
  endinterface

  public :: allocate
  interface allocate
     module procedure communicator_allocate
  endinterface

  public :: request_border
  interface request_border
     module procedure communicator_request_border
  endinterface

  public :: communicate_particles
  interface communicate_particles
     module procedure communicator_communicate_particles
  endinterface

  public :: communicate_ghosts
  interface communicate_ghosts
     module procedure communicator_communicate_ghosts
  endinterface

  public :: communicate_forces
  interface communicate_forces
     module procedure communicator_communicate_forces
  endinterface

  public :: register
  interface register
     module procedure communicator_register
  endinterface register

contains

  !>
  !! Constructor
  !!
  !! Initialize the parallelization module
  !<
  subroutine communicator_init(this, p, decomposition, verlet_shell, context, error)
    implicit none

    type(communicator_t),        intent(inout) :: this
    type(particles_t), target,   intent(inout) :: p
    integer,           optional, intent(in)    :: decomposition(3)
    real(DP),          optional, intent(in)    :: verlet_shell
    type(MPI_Context), optional, intent(in)    :: context
    integer,           optional, intent(inout) :: error

    ! ---

    integer   :: d
    real(DP)  :: l(3)
    logical   :: periods_for_mpi(3)

    ! ---

    call prlog("- communicator_init -")

    if (present(decomposition)) then
       this%decomposition  = decomposition
    endif

    if (present(verlet_shell)) then
       this%verlet_shell   = verlet_shell
    endif

    call prlog("     decomposition = ( "//this%decomposition//" )")

    if (this%decomposition(1)*this%decomposition(2)*this%decomposition(3) /= mpi_n_procs()) then
       RAISE_ERROR("Decomposition geometry requires " // this%decomposition(1)*this%decomposition(2)*this%decomposition(3) // " processes, however, MPI returns " // mpi_n_procs() // " processes.", error)
    endif

    this%pbc              = p%pbc /= 0
    this%requested_border = 0.0_DP

    periods_for_mpi = (/ .true., .true., .true. /)
    call initialise(this%mpi, &
         dims     = this%decomposition, &
         periods  = periods_for_mpi, &
         context  = context, &
         error    = error)
    PASS_ERROR(error)

    call prlog("     coords        = ( "//this%mpi%my_coords//" )")

    do d = 1, 3
       call cart_shift( &
            this%mpi, d-1, 1, this%l(d), this%r(d), error=error)
       PASS_ERROR(error)
    enddo

    l        = (/ p%Abox(1, 1), p%Abox(2, 2), p%Abox(3, 3) /)

    p%lower  = this%mpi%my_coords     * l/this%decomposition
    p%upper  = (this%mpi%my_coords+1) * l/this%decomposition

    call prlog("     lower         = ( "//p%lower//" )")
    call prlog("     upper         = ( "//p%upper//" )")

    this%off_r  = 0.0_DP
    this%off_l  = 0.0_DP

    do d = 1, 3
       if (p%pbc(d) /= 0 .and. this%decomposition(d) > 1) then
          ! No pbcity in binning because we explicitly copy the atoms
          ! from the other processors
          p%locally_pbc(d) = .false.
       else
          this%pbc(d) = .false.
       endif

       if (this%mpi%my_coords(d) == 0) then
          this%off_l(d) = -l(d)
       else if (this%mpi%my_coords(d) == this%decomposition(d)-1) then
          this%off_r(d) = l(d)
       endif
    enddo

    call prlog("     pbc (global)  = ( "//(p%pbc /= 0)//" )")
    call prlog("     pbc (par.)    = ( "//this%pbc//" )")
    call prlog("     pbc (local)   = ( "//p%locally_pbc//" )")

    call prlog("     off_l         = ( "//this%off_l//" )")
    call prlog("     off_r         = ( "//this%off_r//" )")

    this%n_send_p_tot  = 0
    this%n_recv_p_tot  = 0
    this%n_send_g_tot  = 0
    this%n_recv_g_tot  = 0
    this%nit_p         = 0
    this%nit_g         = 0

    call prlog

  endsubroutine communicator_init


  !>
  !! Allocate buffers
  !!
  !! Allocate buffers
  !<
  subroutine communicator_allocate(this, p)
    implicit none

    type(communicator_t), intent(inout)  :: this
    type(particles_t), intent(inout)    :: p

    ! ---

    integer  :: s

    ! ---

    call prlog("- communicator_allocate -")

    call log_memory_start("communicator_allocate")

    call size_by_tag(p%data, F_COMMUNICATE, this%n_particle_data)
    call size_by_tag(p%data, F_COMM_GHOSTS, this%n_ghost_data)
    call size_by_tag(p%data, F_COMM_FORCES, this%n_force_data)
    ! Additional rank information
    this%n_ghost_data  = this%n_ghost_data + 1

    call prlog("     n_particle_data  = " // this%n_particle_data)
    call prlog("     n_ghost_data     = " // this%n_ghost_data)
    call prlog("     n_force_data     = " // this%n_force_data)

    s = max(this%n_particle_data, this%n_ghost_data) * p%maxnatloc

    allocate(this%send_l(s))
    allocate(this%send_r(s))
    allocate(this%recv_l(s))
    allocate(this%recv_r(s))

    allocate(this%ghosts_r(p%maxnatloc))
    allocate(this%ghosts_l(p%maxnatloc))

    call log_memory_estimate(this%send_l)
    call log_memory_estimate(this%send_r)
    call log_memory_estimate(this%recv_l)
    call log_memory_estimate(this%recv_r)

    call log_memory_estimate(this%ghosts_r)
    call log_memory_estimate(this%ghosts_l)

    call log_memory_stop("communicator_allocate")

    !
    ! Copy ghost particles (for first integration step)
    !

    call communicate_ghosts(mod_communicator, p, .true.)
    if (this%communicate_forces) then
       call communicate_forces(mod_communicator, p)
    endif

    call prlog

  endsubroutine communicator_allocate


  !>
  !! Set the communication border
  !!
  !! Set the communication border
  !<
  subroutine communicator_request_border(this, p, border, verlet_shell, error)
    implicit none

    type(communicator_t), intent(inout) :: this
    type(particles_t),   intent(inout) :: p
    real(DP),            intent(in)    :: border
    real(DP), optional,  intent(in)    :: verlet_shell
    integer,  optional,  intent(inout) :: error

    ! ---

    integer  :: d

    ! ---

    call prlog("- communicator_request_border -")

    if (present(verlet_shell)) then
       this%verlet_shell  = verlet_shell
    endif

    this%requested_border  = max(this%requested_border, border)
    this%border            = this%requested_border + this%verlet_shell

    call prlog("     requested_border  = "//this%requested_border)
    call prlog("     verlet_shell      = "//this%verlet_shell)
    call prlog("     border            = "//this%border)

    if (any(this%pbc .and. (this%decomposition .gt. 1) .and. (p%upper - p%lower < 2*this%border))) then
       RAISE_ERROR("Domain smaller than twice the border. This does not work (yet).", error)
    else if (any(p%upper - p%lower < 2*this%border)) then
       call prlog("    (Attention: Domain smaller than twice the border in at least one direction)")
    endif

    do d = 1, 3
       if (this%pbc(d) .or. this%mpi%my_coords(d) /= 0) then
          p%lower_with_border(d)  = p%lower(d) - this%border
       else
          p%lower_with_border(d)  = p%lower(d)
       endif
       if (this%pbc(d) .or. this%mpi%my_coords(d) /= this%decomposition(d)-1) then
          p%upper_with_border(d)  = p%upper(d) + this%border
       else
          p%upper_with_border(d)  = p%upper(d)
       endif
    enddo

    call prlog("     lower_with_border  = ( "//p%lower_with_border//" )")
    call prlog("     upper_with_border  = ( "//p%upper_with_border//" )")

    call prlog   
 
  endsubroutine communicator_request_border

  
  !>
  !! Destructor
  !!
  !! Delete all communicators
  !<
  subroutine communicator_del(this)
    implicit none

    type(communicator_t), intent(inout)  :: this

    ! ---

    call prlog("- communicator_del -")

    deallocate(this%send_l)
    deallocate(this%send_r)
    deallocate(this%recv_l)
    deallocate(this%recv_r)

    deallocate(this%ghosts_l)
    deallocate(this%ghosts_r)

    call prlog("     Average number of particles sent/received per iteration:")
    call prlog("     Particles send  = "//(1.0_DP*this%n_send_p_tot)/this%nit_p)
    call prlog("     Particles recv  = "//(1.0_DP*this%n_recv_p_tot)/this%nit_p)
    call prlog("     Ghosts send     = "//(1.0_DP*this%n_send_g_tot)/this%nit_g)
    call prlog("     Ghosts recv     = "//(1.0_DP*this%n_recv_g_tot)/this%nit_g)
    call prlog

    call finalise(this%mpi)

  endsubroutine communicator_del


  !**********************************************************************
  ! Add particle data to the send buffer
  !**********************************************************************
  subroutine copy_to_send_buffer(p, i, n, buffer)
    implicit none

    type(particles_t), intent(inout)  :: p
    integer, intent(in)               :: i
    integer, intent(inout)            :: n
    real(DP), intent(inout)           :: buffer(:)

    ! ---

    call pack_buffer(p%data, F_COMMUNICATE, i, n, buffer)

    p%global2local(p%index(i)) = 0 ! This one is gone

!    write (ilog, '(5X,A,3I5,6F20.10)')  "send: ", p%index(i), i, p%Z(i), POS3(p, i), PNC3(p, i)

  endsubroutine copy_to_send_buffer


  !**********************************************************************
  ! Copy particle data from the receive buffer
  !**********************************************************************
  subroutine copy_from_recv_buffer(p, n, buffer, off)
    implicit none

    type(particles_t), intent(inout)  :: p
    integer, intent(in)               :: n
    real(DP), intent(in)              :: buffer(:)
    real(DP), intent(in)              :: off(3)

    ! ---

    integer  :: i

    ! ---

!    do i = 1, n
    i = 0
    do while (i < n)
       p%natloc = p%natloc+1

       call unpack_buffer(p%data, F_COMMUNICATE, i, buffer, p%natloc)

       p%Z(p%natloc)                      = p%el2Z(p%el(p%natloc))
       p%global2local(p%index(p%natloc))  = p%natloc

!       write (ilog, '(5X,A,3I5,9F20.10)')  "recv1: ", p%index(p%natloc), p%global2local(p%index(p%natloc)), p%Z(p%natloc), POS3(p, p%natloc), PNC3(p, p%natloc), off(:)

!       POS3(p, p%natloc)  = POS3(p, p%natloc) + off(:)
       PNC3(p, p%natloc)  = PNC3(p, p%natloc) + off

!       write (ilog, '(5X,A,3I5,9F20.10)')  "recv2: ", p%index(p%natloc), p%global2local(p%index(p%natloc)), p%Z(p%natloc), POS3(p, p%natloc), PNC3(p, p%natloc), off(:)
    enddo

  endsubroutine copy_from_recv_buffer


  !**********************************************************************
  ! Communicate particles which left the domains
  ! to the neighboring domains (former order routine)
  !**********************************************************************
  subroutine communicator_communicate_particles(this, p, error)
    implicit none

    type(communicator_t), intent(inout)  :: this
    type(particles_t), intent(inout)    :: p
    integer, intent(out), optional      :: error

    ! ---

    ! 
    ! General and auxiliary variables      
    !

    integer   :: i, d
    integer   :: oldnatloc

    ! 
    ! Structure variables, mpi and system structure
    !

    integer   :: n_send_l, n_send_r, n_recv_l, n_recv_r

    real(DP)  :: off_l(3), off_r(3)

    ! ---

    INIT_ERROR(error)

    call timer_start("communicator_communicate_particles")

    this%nit_p = this%nit_p + 1

    do i = p%natloc+1, p%nat
       p%global2local(p%index(i)) = 0
    enddo

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    do d = 1, 3

       if (this%decomposition(d) > 1) then

          oldnatloc  = p%natloc
          p%natloc   = 0

          n_send_r   = 0
          n_send_l   = 0

          do i = 1, oldnatloc

             if (PNC(p, i, d) >= p%upper(d)) then
                ! Send to the right

                call copy_to_send_buffer(p, i, n_send_r, this%send_r)
!                n_send_r = n_send_r + 1

             else if (PNC(p, i, d) < p%lower(d)) then
                ! Send to the left

                call copy_to_send_buffer(p, i, n_send_l, this%send_l)
!                n_send_l = n_send_l + 1

             else
                ! Keep on this processor and reorder

                p%natloc = p%natloc+1

                if (p%natloc /= i) then
                   call move(p, p%natloc, i)
                endif

             endif

          enddo

          !write (ilog, *)  d, "r: ", n_send_r
          !write (ilog, *)  d, "l: ", n_send_l

          this%n_send_p_tot = this%n_send_p_tot + n_send_r + n_send_l

         call sendrecv(this%mpi, &
               this%send_r(1:n_send_r), this%r(d), 0, &
               this%recv_l(1:this%n_particle_data*p%maxnatloc), this%l(d), 0, &
               n_recv_l, &
               error = error)
          PASS_ERROR(error)

          call sendrecv(this%mpi, &
               this%send_l(1:n_send_l), this%l(d), 1, &
               this%recv_r(1:this%n_particle_data*p%maxnatloc), this%r(d), 1, &
               n_recv_r, &
               error = error)
          PASS_ERROR(error)

          this%n_recv_p_tot = this%n_recv_p_tot + n_recv_r/this%n_particle_data + n_recv_l/this%n_particle_data

          off_l    = 0.0_DP
          ! This will be done by inbox
          off_l(d) = this%off_l(d)

          off_r   = 0.0_DP
          ! This will be done by inbox
          off_r(d) = this%off_r(d)

          call copy_from_recv_buffer(p, n_recv_l, this%recv_l, off_l)
          call copy_from_recv_buffer(p, n_recv_r, this%recv_r, off_r)

       endif

    enddo

    p%nat = p%natloc

    call timer_stop("communicator_communicate_particles")

  endsubroutine communicator_communicate_particles


  !**********************************************************************
  ! Copy particle data to the (ghost) send buffer
  !**********************************************************************
  subroutine copy_to_send_ghosts(mpi, p, i, n, buffer)
    implicit none

    type(MPI_context), intent(in)     :: mpi
    type(particles_t), intent(inout)  :: p
    integer, intent(in)               :: i
    integer, intent(inout)            :: n
    real(DP), intent(inout)           :: buffer(:)

    ! ---

    call pack_buffer(p%data, F_COMM_GHOSTS, i, n, buffer)

    n  = n + 1
    if (i > p%natloc) then
       buffer(n)  = p%from_rank(i)
    else
       buffer(n)  = mpi%my_proc
    endif

!    if (p%index(i) == 914635) then
!       write (ilog, '(5X,A,2I10,6F20.10)')  "g-send: ", p%index(i), p%el(i), POS3(p, i), PNC3(p, i)
!    endif

  endsubroutine copy_to_send_ghosts


  !**********************************************************************
  ! Copy particle data from the (ghost) receive buffer
  !**********************************************************************
  subroutine copy_from_recv_ghosts(p, n, buffer, off)
    implicit none

    type(particles_t), intent(inout)  :: p
    integer, intent(in)               :: n
    real(DP), intent(in)              :: buffer(:)
    real(DP), intent(in)              :: off(3)

    ! ---

    integer  :: i

    ! ---

    i = 0
    do while (i < n)
       p%nat = p%nat+1

       call unpack_buffer(p%data, F_COMM_GHOSTS, i, buffer, p%nat)

       i  = i + 1
       p%from_rank(p%nat)  = buffer(i)

       p%Z(p%nat)                     = p%el2Z(p%el(p%nat))
       p%global2local(p%index(p%nat)) = p%nat

       PNC3(p, p%nat)  = PNC3(p, p%nat) + off
#ifndef IMPLICIT_R
       POS3(p, p%nat)  = in_cell(p, PNC3(p, p%nat))
#endif

#ifdef DEBUG
       if (.not. ( all(POS3(p, p%nat) >= 0.0_DP) .and. all(POS3(p, p%nat) < (/ p%Abox(1, 1), p%Abox(2, 2), p%Abox(3, 3) /)) )) then
          call particles_dump_info(p, p%nat)
          EXIT_ON_ERROR("Particle outside of the simulation domain.", i)
       endif
#endif

!       if (p%index(p%nat) == 914635) then
!          write (ilog, '(5X,A,2I5,9F20.10)')  "g-recv: ", p%index(p%nat), p%el(p%nat), POS3(p, p%nat), PNC3(p, p%nat), off(:)
!       endif
    enddo

  endsubroutine copy_from_recv_ghosts


  !**********************************************************************
  ! Communicate ghost particles to neighboring domains
  ! (former prebinning routine). *bwidth* is the width of
  ! the border.
  !**********************************************************************
  subroutine communicator_communicate_ghosts(this, p, new_list, error)
    implicit none

    type(communicator_t), intent(inout)  :: this
    type(particles_t), intent(inout)    :: p
    logical, intent(in)                 :: new_list
    integer, intent(inout), optional    :: error

    ! ---

    real(DP)  :: upper(3), lower(3)
    integer   :: i, d, list_off_r, list_off_l, n_send_r, n_send_l, n_recv_r, n_recv_l

    real(DP)  :: off_l(3), off_r(3)

    ! ---

    INIT_ERROR(error)

    call timer_start("communicator_communicate_ghosts")

    this%nit_g = this%nit_g + 1

    do d = 1, 3

       if (this%pbc(d) .or. this%mpi%my_coords(d) /= 0) then
          lower(d)  = p%lower(d) + this%border
       else
          lower(d)  = p%lower(d)
       endif

       if (this%pbc(d) .or. this%mpi%my_coords(d) /= this%decomposition(d)-1) then
          upper(d)  = p%upper(d) - this%border
       else
          upper(d)  = p%upper(d)
       endif

    enddo

    do i = p%natloc+1, p%nat
       p%global2local(p%index(i)) = 0
    enddo

    p%nat  = p%natloc

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    list_off_r = 0
    list_off_l = 0
    do d = 1, 3

       if (this%decomposition(d) > 1) then

          n_send_r  = 0
          n_send_l  = 0

          if (new_list) then

             this%n_ghosts_r(d)  = 0
             this%n_ghosts_l(d)  = 0

             do i = 1, p%nat
                if (PNC(p, i, d) >= upper(d)) then
                   call copy_to_send_ghosts(this%mpi, p, i, n_send_r, this%send_r)

                   this%n_ghosts_r(d)                            = this%n_ghosts_r(d)+1
                   this%ghosts_r(list_off_r+this%n_ghosts_r(d))  = p%index(i)

                else if (PNC(p, i, d) < lower(d)) then
                   call copy_to_send_ghosts(this%mpi, p, i, n_send_l, this%send_l)

                   this%n_ghosts_l(d)                            = this%n_ghosts_l(d)+1
                   this%ghosts_l(list_off_l+this%n_ghosts_l(d))  = p%index(i)

                endif
             enddo

          else

             do i = 1, this%n_ghosts_r(d)
                call copy_to_send_ghosts(this%mpi, p, p%global2local(this%ghosts_r(list_off_r+i)), n_send_r, this%send_r)
             enddo

             do i = 1, this%n_ghosts_l(d)
                call copy_to_send_ghosts(this%mpi, p, p%global2local(this%ghosts_l(list_off_l+i)), n_send_l, this%send_l)
             enddo

          endif

          this%n_send_g_tot = this%n_send_g_tot + this%n_ghosts_r(d) + this%n_ghosts_l(d)

          call sendrecv(this%mpi, &
               this%send_r(1:n_send_r), this%r(d), 0, &
               this%recv_l(1:this%n_ghost_data*p%maxnatloc), this%l(d), 0, &
               n_recv_l, &
               error = error)
          PASS_ERROR(error)

          call sendrecv(this%mpi, &
               this%send_l(1:n_send_l), this%l(d), 1, &
               this%recv_r(1:this%n_ghost_data*p%maxnatloc), this%r(d), 1, &
               n_recv_r, &
               error = error)
          PASS_ERROR(error)

          this%n_recv_g_tot = this%n_recv_g_tot + n_recv_r/this%n_ghost_data + n_recv_l/this%n_ghost_data

          off_l    = 0.0_DP
          off_l(d) = this%off_l(d)

          off_r    = 0.0_DP
          off_r(d) = this%off_r(d)

          call copy_from_recv_ghosts(p, n_recv_l, this%recv_l, off_l)
          call copy_from_recv_ghosts(p, n_recv_r, this%recv_r, off_r)

          list_off_r = list_off_r + this%n_ghosts_r(d)
          list_off_l = list_off_l + this%n_ghosts_l(d)

       endif

    enddo

    call timer_stop("communicator_communicate_ghosts")

  endsubroutine communicator_communicate_ghosts


  !**********************************************************************
  ! Copy forces to the (ghost) send buffer
  !**********************************************************************
  subroutine copy_forces_to_send_ghosts(p, i, n, buffer)
    implicit none

    type(particles_t), intent(inout)  :: p
    integer, intent(in)               :: i
    integer, intent(in)               :: n
    real(DP), intent(inout)           :: buffer(:)

    ! ---

    integer  :: m

    ! ---

    m = n
    call pack_buffer(p%data, F_COMM_FORCES, i, m, buffer)

    !write (ilog, '(A,I5,3F20.10)')  "Fsend: ", p%index(i), FOR3(p, i)

  endsubroutine copy_forces_to_send_ghosts
  

  !**********************************************************************
  ! Copy particle data from the (ghost) receive buffer
  !**********************************************************************
  subroutine copy_forces_from_recv_ghosts(p, cur, n, buffer)
    implicit none

    type(particles_t), intent(inout)  :: p
    integer, intent(inout)            :: cur
    integer, intent(in)               :: n
    real(DP), intent(in)              :: buffer(:)

    ! ---

    integer  :: i

    ! ---

    i = 0
    do while (i < n)
       cur = cur+1

       call unpack_buffer(p%data, F_COMM_FORCES, i, buffer, cur)
    enddo

  endsubroutine copy_forces_from_recv_ghosts


  !**********************************************************************
  ! Communicate forces of ghost particles back.
  ! This is needed for rigid object (i.e., water) or to reduce the
  ! border size in BOPs.
  !**********************************************************************
  subroutine communicator_communicate_forces(this, p, error)
    implicit none

    type(communicator_t), intent(inout)  :: this
    type(particles_t), intent(inout)    :: p
    integer, intent(inout), optional    :: error

    ! ---

    real(DP)  :: upper(3), lower(3)
    integer   :: i, d, list_off_r, list_off_l, n_recv_r, n_recv_l, cur

    ! ---

    call timer_start("communicator_communicate_forces")

    do d = 1, 3

       if (this%pbc(d) .or. this%mpi%my_coords(d) /= 0) then
          lower(d)  = p%lower(d) + this%border
       else
          lower(d)  = p%lower(d)
       endif

       if (this%pbc(d) .or. this%mpi%my_coords(d) /= this%decomposition(d)-1) then
          upper(d)  = p%upper(d) - this%border
       else
          upper(d)  = p%upper(d)
       endif

    enddo

    !
    ! Loop over dimensions and distribute particle in the
    ! respective direction
    !

    list_off_r = 0
    list_off_l = 0
    cur = p%natloc
    do d = 1, 3

       if (this%decomposition(d) > 1) then

          do i = 1, this%n_ghosts_r(d)
             call copy_forces_to_send_ghosts(p, p%global2local(this%ghosts_r(list_off_r+i)), (i-1)*this%n_force_data, this%send_r)
          enddo

          do i = 1, this%n_ghosts_l(d)
             call copy_forces_to_send_ghosts(p, p%global2local(this%ghosts_l(list_off_l+i)), (i-1)*this%n_force_data, this%send_l)
          enddo

          call sendrecv(this%mpi, &
               this%send_r(1:this%n_force_data*this%n_ghosts_r(d)), this%r(d), 0, &
               this%recv_l(1:this%n_force_data*p%maxnatloc), this%l(d), 0, &
               n_recv_l, &
               error = error)
          PASS_ERROR(error)

          call sendrecv(this%mpi, &
               this%send_l(1:this%n_force_data*this%n_ghosts_l(d)), this%l(d), 1, &
               this%recv_r(1:this%n_force_data*p%maxnatloc), this%r(d), 1, &
               n_recv_r, &
               error = error)
          PASS_ERROR(error)

          call copy_forces_from_recv_ghosts(p, cur, n_recv_l, this%recv_l)
          call copy_forces_from_recv_ghosts(p, cur, n_recv_r, this%recv_r)

          list_off_r = list_off_r + this%n_ghosts_r(d)
          list_off_l = list_off_l + this%n_ghosts_l(d)

       endif

    enddo

    call timer_stop("communicator_communicate_forces")

  endsubroutine communicator_communicate_forces


  !>
  !! Register object introspection
  !!
  !! Expose state of communicator object through a dictionary object
  !<
  subroutine communicator_register(this, cfg)
    use, intrinsic :: iso_c_binding

    implicit none

    type(communicator_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)             :: cfg

    ! ---

    type(c_ptr)  :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("DomainDecomposition"), &
         CSTR("Domain decomposition module."))

    call ptrdict_register_intpoint_property(m, c_loc(this%decomposition(1)), &
         CSTR("decomposition"), &
         CSTR("Number of domains in each direction, i.e. type of the decomposition."))

  endsubroutine communicator_register

endmodule communicator

#endif
