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
!   shared
!   classtype:sliding_p_t classname:SlidingP interface:potentials
! @endmeta

!>
!! Sliding friction barostat
!!
!! Apply a constant pressure in z-direction and slide the
!! system perependicular at a specified angle. The full
!! details of this barostat are described in:
!!
!! L. Pastewka, S. Moser, M. Moseler,
!! Tribol. Lett. 39, 49 (2010)
!!   doi: 10.1007/s11249-009-9566-8
!!
!! If logging is turned on, a file (sliding_p.out) will contain the
!! forces exerted by the other atoms in the system on the moving
!! rigid slab on top (marked as px, py, pz). Here, p_i = F_i / A,
!! where A is the lateral size of the box (A = Lx*Ly). Additionally,
!! the heigh of the box (h) is included, defined as the difference
!! between the highest and lowest atomic z-coordinate.
!!
!! So, for sliding along x, the shear stress is immediately given by
!! \f[
!!   \tau = \frac{F}{A} = |px|,
!! \f]
!! and the friction coefficient is obtained as
!! \f[
!!   \mu = \frac{\tau}{|pz|}.
!! \f]
!<

#include "macros.inc"

module sliding_p
  use supplib

  use io
  use logging
  use timer

  use data
  use particles
  use neighbors
  use dynamics

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  character(MAX_NAME_STR), parameter  :: STRESS_STR  = "SlidingP::stress"
  character(MAX_NAME_STR), parameter  :: HEIGHT_STR  = "SlidingP::height"

  public :: sliding_p_t
  type sliding_p_t

     real(DP)  :: Pz     = 1.0_DP        !< Normal pressure

     real(DP)  :: abs_v  = -1.0_DP       !< Sliding velocity (absolute value)
     real(DP)  :: angle  = 0.0_DP        !< Angle of sliding with respect to the x-axis
     real(DP)  :: v(3)                   !< Velocity vector

     real(DP)  :: h0     = -1.0_DP       !< Initial sample height
     real(DP)  :: C11    = 1.0_DP        !< Young modulus for compression in z-direction

     real(DP)  :: fC     = -1.0_DP       !< Recurrence frequency, autocomputed if < 0

     real(DP)  :: p      = 0.2_DP        !< Transmission function at recurrence frequency

     real(DP)  :: divM    = 1.0_DP       !< Divides M, the total mass of the top atoms, by this factor
     integer   :: calc    = 1            !< 1=use p and v to calculate M and gamma; 2=directly read divM and gamma from input file

     integer   :: bot    = -1            !< Group of fixed atoms
     integer   :: n_bot
     logical(BOOL)   :: reset_top_vel  = .false.   !< Reset the velocity of the top moving atoms
     integer   :: top    = -1            !< Group of moving atoms
     integer   :: n_top
     real(DP)  :: M                      !< Total mass of top atoms (selected by group number "top"); they move as if having the (individually set) mass M

     real(DP)  :: gamma   = 0.1_DP       !< Dissipation constant
     real(DP)  :: mass                   !< Total mass of moving atoms

     logical(BOOL)  :: log  = .false.          !< Log presure to a file
     integer  :: un

     real(DP)   :: v_top(3)              !< Current velocity of the top slab

     real(DP), pointer  :: v_arr(:, :)   !< Velocities

     real(DP), pointer  :: stress(:) => NULL()
     real(DP), pointer  :: height => NULL()

  endtype sliding_p_t

  !
  ! Values for "calc": How to calculate/get M and gamma
  !

  integer, parameter  :: PARAM_CALCULATEDAMPING = 1
  integer, parameter  :: PARAM_SETDAMPING       = 2
  integer, parameter  :: PARAM_MASSLESSDAMPING  = 3
  integer, parameter  :: PARAM_CONSTANTHEIGHT   = 4

  public :: register_data
  interface register_data
     module procedure sliding_p_register_data
  endinterface register_data

  public :: init
  interface init
     module procedure sliding_p_init
  endinterface

  public :: set
  interface set
     module procedure sliding_p_init
  endinterface

  public :: del
  interface del
     module procedure sliding_p_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure sliding_p_bind_to
  endinterface

  public :: adjust_velocities_and_forces
  interface adjust_velocities_and_forces
     module procedure sliding_p_adjust_velocities_and_forces
  endinterface

  public :: energy_and_forces_with_dyn
  interface energy_and_forces_with_dyn
     module procedure sliding_p_energy_and_forces
  endinterface

  public :: register
  interface register
    module procedure sliding_p_register
  endinterface

contains

  !>
  !! Register any data columns with the Atoms object
  !<
  subroutine sliding_p_register_data(this, p, ierror)
    implicit none

    type(sliding_p_t),           intent(inout)  :: this
    type(particles_t),           intent(inout)  :: p
    integer,           optional, intent(out)    :: ierror

    ! ---

    INIT_ERROR(ierror)

    call add_real3_attr(p%data, STRESS_STR, F_TO_TRAJ+F_TO_ENER, ierror=ierror)
    PASS_ERROR(ierror)
    call add_real_attr(p%data, HEIGHT_STR, F_TO_TRAJ+F_TO_ENER, ierror=ierror)
    PASS_ERROR(ierror)

  endsubroutine sliding_p_register_data


  !>
  !! Constructor
  !!
  !! Construct a SlidingP object
  !<
  subroutine sliding_p_init(this, Pz, abs_v, angle, h0, C11, fC, p, divM, gamma, calc, bot, top, reset_top_vel, log)
    implicit none

    type(sliding_p_t), intent(inout)  :: this
    real(DP), intent(in), optional    :: Pz
    real(DP), intent(in), optional    :: abs_v
    real(DP), intent(in), optional    :: angle
    real(DP), intent(in), optional    :: h0
    real(DP), intent(in), optional    :: C11
    real(DP), intent(in), optional    :: fC
    real(DP), intent(in), optional    :: p
    real(DP), intent(in), optional    :: divM
    real(DP), intent(in), optional    :: gamma
    integer, intent(in), optional     :: calc
    integer, intent(in), optional     :: bot
    integer, intent(in), optional     :: top
    logical, intent(in), optional     :: reset_top_vel
    logical, intent(in), optional     :: log

    ! ---

    call prlog("- sliding_p_init -")

    ASSIGN_PROPERTY(Pz)
    ASSIGN_PROPERTY(abs_v)
    ASSIGN_PROPERTY(angle)

    ASSIGN_PROPERTY(h0)
    ASSIGN_PROPERTY(C11)

    ASSIGN_PROPERTY(fC)
    ASSIGN_PROPERTY(p)

    ASSIGN_PROPERTY(divM)
    ASSIGN_PROPERTY(gamma)
    ASSIGN_PROPERTY(calc)

    ASSIGN_PROPERTY(bot)
    ASSIGN_PROPERTY(top)
    ASSIGN_PROPERTY(reset_top_vel)

    ASSIGN_PROPERTY(log)

    call prlog

  endsubroutine sliding_p_init
  

  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine sliding_p_del(this)
    implicit none

    type(sliding_p_t), intent(inout)  :: this

    ! ---

#ifdef _MP
    if (mod_communicator%mpi%my_proc == ROOT) then
#endif

    if (this%log) then
       call fclose(this%un)
    endif

#ifdef _MP
    endif
#endif

  endsubroutine sliding_p_del


  !>
  !! Initialize a sliding_p object
  !!
  !! Initialize a sliding_p object
  !<
  subroutine sliding_p_bind_to(this, p, nl, ierror)
    implicit none

    type(sliding_p_t), intent(inout)  :: this
    type(particles_t), intent(in)     :: p
    type(neighbors_t), intent(in)     :: nl
    integer, intent(inout), optional  :: ierror

    ! ---

    integer   :: i

    real(DP)  :: l, lx, ly, a, k, h1, h2, h, fC, wC, curM

    ! ---

    call prlog("- sliding_p_bind_to -")

    call prlog("     Pz       = " // this%Pz)
    call prlog("     abs_v    = " // this%abs_v)
    call prlog("     angle    = " // this%angle)
    call prlog("     h0       = " // this%h0)
    call prlog("     C11      = " // this%C11)
    call prlog("     p        = " // this%p)
    select case (this%calc)
     case (PARAM_CALCULATEDAMPING)
      call prlog("     calc     = " // this%calc // " (Use 'v' and 'p' to calculate 'gamma' and 'M'.)")
     case (PARAM_SETDAMPING)
      call prlog("     calc     = " // this%calc // " (Calculate 'M' via 'divM' and 'gamma' is set.)")
     case (PARAM_MASSLESSDAMPING)
      call prlog("     calc     = " // this%calc // " (Calculate 'M' via 'divM' and 'gamma' is set, multiplicated by mass of top atoms.)")
     case (PARAM_CONSTANTHEIGHT)
      call prlog("     calc     = " // this%calc // " (Constant height calculation.)")
     case default
       RAISE_ERROR("The parameter 'calc' (=" // this%calc // ") is not valid! It must be either 1 (set 'p' and 'v') or 2 or 3 (set 'divM' and 'gamma')!", ierror)
    end select
    call prlog("     bot      = " // this%bot)
    call prlog("     top      = " // this%top)

    this%n_bot  = count(p%g(1:p%natloc) == this%bot)
    this%n_top  = count(p%g(1:p%natloc) == this%top)

#ifdef _MP
    call sum_in_place(mod_communicator%mpi, this%n_bot)
    call sum_in_place(mod_communicator%mpi, this%n_top)
#endif

    if (this%n_bot == 0) then
       RAISE_ERROR("No *bottom* atoms found.", ierror)
    endif

    if (this%n_top == 0) then
       RAISE_ERROR("No *top* atoms found.", ierror)
    endif

    call prlog("     * n_bot  = " // this%n_bot)
    call prlog("     * n_top  = " // this%n_top)

    call ptr_by_name(p%data, V_STR, this%v_arr)

    h1  = minval(POS(p, 1:p%natloc, 3), p%g(1:p%natloc) == this%top)
    h2  = maxval(POS(p, 1:p%natloc, 3), p%g(1:p%natloc) == this%bot)

#ifdef _MP
    h1 = min(mod_communicator%mpi, h1)
    h2 = max(mod_communicator%mpi, h2)
#endif

    h   = h1 - h2

    if (this%h0 > 0.0_DP) then
       h  = this%h0
    endif

    if (h < 0.0_DP) then
       RAISE_ERROR("Computed height of the simulation cell is < 0. Maybe you switched top and bottom groups?", ierror)
    endif

    call prlog("     * h      = " // h)

    a   = this%angle * PI / 180
    lx  = p%Abox(1, 1)
    ly  = p%Abox(2, 2)

    if (this%fC > 0.0_DP) then
       fC  = this%fC
    else

       l   = sqrt( (lx*cos(a))**2 + (ly*sin(a))**2 )
       fC  = abs(this%abs_v)/l

       call prlog("     * fC     = " // fC)

    endif

    k       = this%C11 * lx * ly / h
    wC      = 2*PI*fC

    curM    = sum(p%m(1:p%natloc), p%g(1:p%natloc) == this%top)

    select case (this%calc)
     case (PARAM_CALCULATEDAMPING)
      this%M  = k / wC**2 * sqrt( 1.0_DP/this%p**2 - 1 )
     case (PARAM_SETDAMPING:PARAM_MASSLESSDAMPING)
      if (this%divM /= 0.0_DP) then
       this%M  = curM / this%divM
      else
       this%M  = curM / 1.00_DP
      end if
     case default
       this%M  = curM / 1.00_DP ! This case should not happen (not defined variable "calc")
    end select

#ifdef _MP
    call sum_in_place(mod_communicator%mpi, curM)
#endif

    call prlog("     * M      = " // this%M // " ( current mass of " // curM // " times " // (this%M/curM) // " )")

    select case (this%calc)
     case (PARAM_CALCULATEDAMPING)
      this%gamma  = sqrt(2*this%M*k)
     case (PARAM_SETDAMPING)
      ! this%gamma already read and set from input file
     case (PARAM_MASSLESSDAMPING)
      ! this%gamma already read and set from input file, but will be changed now, according to the mass of top atoms, so that the damping will be independant of the mass
      this%gamma = this%gamma * this%M
     case default
      this%gamma = 0.0_DP ! This case should not happen (not defined variable "calc")
    end select

    select case (this%calc)
     case (PARAM_CALCULATEDAMPING, PARAM_MASSLESSDAMPING)
      call prlog("     * gamma  = " // this%gamma)
     case (PARAM_SETDAMPING)
      ! this%gamma already read and set from input file
      call prlog("     gamma  = " // this%gamma)
     case default
      call prlog("     gamma  = " // this%gamma)
    end select

    this%v  = (/ this%abs_v * cos(a), this%abs_v * sin(a), 0.0_DP /)

    if (this%abs_v <= 0.0_DP) &
         this%v  = 0.0_DP

#ifdef _MP
    if (mpi_id() == ROOT) then
#endif

    if (this%log) then
       this%un  = fopen("sliding_p.out", F_WRITE)
       write (this%un, '(A6,14X,4A20)')  "#01:ti", "02:px", "03:py", "04:pz", "05:h"
    endif

#ifdef _MP
    endif
#endif

    call prlog("     * v      = " // this%v)

!    call adjust_velocities_and_forces(this, p, f)

    if (this%calc == PARAM_CONSTANTHEIGHT) then
       do i = 1, p%natloc
          if (p%g(i) == this%top) then
             VEC3(this%v_arr, i) = this%v
          endif
       enddo

       this%v_top = this%v
    else
       if (this%reset_top_vel) then
          call prlog("     Resetting slider velocity.")
          
          do i = 1, p%natloc
             if (p%g(i) == this%top) then
                VEC3(this%v_arr, i) = 0.0_DP
             endif
          enddo
       endif

       this%v_top = 0.0_DP
    endif

    call attr_by_name(p%data, STRESS_STR, this%stress)
    call attr_by_name(p%data, HEIGHT_STR, this%height)

    call prlog

  endsubroutine sliding_p_bind_to


  !>
  !! Invoke pressure coupling
  !!
  !! Invoke pressure coupling
  !<
  subroutine sliding_p_adjust_velocities_and_forces(this, p, f, ti)
    implicit none

    type(sliding_p_t), intent(inout)  :: this
    type(particles_t), intent(in)     :: p
    real(DP), intent(inout)           :: f(3, p%maxnatloc)
    real(DP), intent(in), optional    :: ti

    ! ---

    integer   :: i
    real(DP)  :: f_top(3), vz, az, h1, h2, h

    ! ---

    h1  = minval(POS(p, 1:p%natloc, 3), p%g(1:p%natloc) == this%top)
    h2  = maxval(POS(p, 1:p%natloc, 3), p%g(1:p%natloc) == this%bot)

#ifdef _MP
    h1 = min(mod_communicator%mpi, h1)
    h2 = max(mod_communicator%mpi, h2)
#endif

    h  = h1 - h2

    f_top  = 0.0_DP
    do i = 1, p%natloc
       if (p%g(i) == this%top) then
          f_top  = f_top + VEC3(f, i)
       endif
    enddo

#ifdef _MP
    call sum_in_place(mod_communicator%mpi, f_top)

    if (mpi_id() == ROOT) then
#endif

    this%stress = f_top / (p%Abox(1, 1)*p%Abox(2, 2))
    this%height = h

    if (this%log) then
       if (present(ti)) then
          write (this%un, '(5ES20.10)')  ti, this%stress, h
       else
          write (this%un, '(4ES20.10)')  this%stress, h
       endif
    endif

#ifdef _MP
    endif
#endif

    if (this%calc == PARAM_CONSTANTHEIGHT) then

       do i = 1, p%natloc
          if (p%g(i) == this%bot) then
             VEC3(this%v_arr, i) = 0.0_DP
             VEC3(f, i)          = 0.0_DP
          else if (p%g(i) == this%top) then
             VEC3(this%v_arr, i) = this%v
          endif
       enddo

       this%v_top = this%v

    else

       vz  = sum(VEC(this%v_arr, 1:p%natloc, 3), p%g == this%top) / this%n_top
#ifdef _MP
       call sum_in_place(mod_communicator%mpi, vz)
#endif
       az  = ( f_top(3) - this%Pz*p%Abox(1, 1)*p%Abox(2, 2) - this%gamma*vz )

       do i = 1, p%natloc
          if (p%g(i) == this%bot) then
             VEC3(this%v_arr, i)  = 0.0_DP
             VEC3(f, i)           = 0.0_DP
          else if (p%g(i) == this%top) then
             VEC(this%v_arr, i, 1:2)  = this%v(1:2)
             VEC(this%v_arr, i, 3)    = vz
             VEC3(f, i)               = (/ 0.0_DP, 0.0_DP, (p%m(i)*az)/this%M /)
          endif
       enddo

       this%v_top(1:2)  = this%v(1:2)
       this%v_top(3)    = vz

    endif

  endsubroutine sliding_p_adjust_velocities_and_forces


  !>
  !! Invoke pressure coupling
  !!
  !! Invoke pressure coupling
  !<
  subroutine sliding_p_energy_and_forces(this, dyn, nl, ierror)
    implicit none

    type(sliding_p_t), intent(inout)   :: this
    type(dynamics_t), target           :: dyn 
    type(neighbors_t), intent(in)      :: nl
    integer, intent(inout), optional   :: ierror

    ! ---

    call timer_start("sliding_p_force")

    call adjust_velocities_and_forces(this, dyn%p, dyn%f, dyn%ti)

    call timer_stop("sliding_p_force")

  endsubroutine sliding_p_energy_and_forces


  subroutine sliding_p_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(sliding_p_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)              :: cfg
    type(c_ptr), intent(out)             :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("SlidingP"), &
         CSTR("Invoke constant sliding velocity and adjustment of normal load."))

    call ptrdict_register_real_property(m, c_loc(this%Pz), CSTR("Pz"), &
         CSTR("Normal load."))

    call ptrdict_register_real_property(m, c_loc(this%abs_v), CSTR("v"), &
         CSTR("Sliding velocity."))
    call ptrdict_register_real_property(m, c_loc(this%angle), CSTR("angle"), &
         CSTR("Sliding angle relative to the x-axis (in degrees)."))

    call ptrdict_register_real_property(m, c_loc(this%h0), CSTR("h0"), &
         CSTR("Height of the sample."))
    call ptrdict_register_real_property(m, c_loc(this%C11), CSTR("C11"), &
         CSTR("C11 elastic constant."))

    call ptrdict_register_real_property(m, c_loc(this%fC), CSTR("fC"), &
         CSTR("Recurrence frequency fC (compute automatically if < 0)."))
    call ptrdict_register_real_property(m, c_loc(this%p), CSTR("p"), &
         CSTR("Value of the transfer function at fC. (used if calc = 1)"))

    call ptrdict_register_real_property(m, c_loc(this%divM), CSTR("divM"), &
         CSTR("Divides M, total mass of top atoms, by this factor. (used if calc = 2)"))
    call ptrdict_register_real_property(m, c_loc(this%gamma), CSTR("gamma"), &
         CSTR("Damping parameter for the harmonic oscillations of two bulks (top and bottom) hitting each other. (used if calc = 2)"))
    call ptrdict_register_integer_property(m, c_loc(this%calc), CSTR("calc"), &
         CSTR("Chose parameter calculation: (1 = use p and v, better for sliding; 2 = use divM and gamma, better for pressing)."))

    call ptrdict_register_integer_property(m, c_loc(this%bot), CSTR("bot"), &
         CSTR("Bottom atoms (group)."))
    call ptrdict_register_integer_property(m, c_loc(this%top), CSTR("top"), &
         CSTR("Top atoms (group). These atoms are moved at the constant velocity."))
    call ptrdict_register_boolean_property(m, c_loc(this%reset_top_vel), &
         CSTR("reset_top_vel"), &
         CSTR("Reset the velocity of the top moving atoms. This is necessay when changing the *p* parameters such that the higher mass of the top atoms does not crush the system."))

    call ptrdict_register_boolean_property(m, c_loc(this%log), CSTR("log"), &
         CSTR("Log forces to 'sliding_p.out'."))

  endsubroutine sliding_p_register

endmodule sliding_p
