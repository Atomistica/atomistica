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
!   dependencies:verlet_support.f90
!   classtype:sliding_t_t classname:SlidingT interface:integrators
! @endmeta

!>
!! Langevin dynamics relative to the sliding velocity
!!
!! Langevin dynamics relative to the velocities of a group of top atoms,
!! and relative to zero sliding velocity for the bottom atoms. This is intended
!! to be used in connection with SlidingP which sets the sliding velocities
!! that are used for thermalization.
!<

#include "macros.inc"
#include "filter.inc"

module sliding_t
  use supplib
  use rng

  use particles
  use filter

  use verlet_support

  use sliding_p

#ifdef _MP
  use communicator
#endif

  private

  public :: sliding_t_t
  type sliding_t_t

     integer    :: bot = 2
     integer    :: top = 3

     real(DP)   :: T            = -1.0
     real(DP)   :: dT           = 0.0

     real(DP)   :: dissipation  = 1.0
     real(DP)   :: tau          = -1.0

  endtype sliding_t_t


  public :: init
  interface init
     module procedure sliding_t_init
  endinterface

  public :: del
  interface del
     module procedure sliding_t_del
  endinterface

  public :: step1_with_barostat, sliding_t_step1
  interface step1_with_barostat
     module procedure sliding_t_step1
  endinterface

  public :: step2
  interface step2
     module procedure sliding_t_step2
  endinterface

  public :: register
  interface register
    module procedure sliding_t_register
  endinterface

contains

  !>
  !! Sliding_t init
  !!
  !! Sliding_t init. Either dissipation (=1/tau) or tau should be
  !! specified, not both. If nothing is given, the default value of dissipation
  !! is used. If this%tau is set, then that is used.
  !<
  subroutine sliding_t_init(this, p, T, dT, dissipation, tau, error)
    implicit none

    type(sliding_t_t), intent(inout)  :: this
    type(particles_t), intent(inout)  :: p
    real(DP), intent(in), optional    :: T
    real(DP), intent(in), optional    :: dT
    real(DP), intent(in), optional    :: dissipation
    real(DP), intent(in), optional    :: tau
    integer, intent(out), optional    :: error

    ! ---

    INIT_ERROR(error)

    ! Checks
    if (present(dissipation) .and. present(tau)) then
       RAISE_ERROR("Please specify either *dissipation* or *tau*.", error)
    endif

    ! Init
    if (present(T)) then
       this%T  = T
    endif
    if (present(dT)) then
       this%dT  = dT
    endif
    if (present(dissipation)) then
       this%dissipation  = dissipation
       this%tau = -1.0_DP
    endif
    if (present(tau)) then
       this%tau = tau
    endif

    if (this%tau > 0.0_DP) then
       this%dissipation  = 1.0_DP/this%tau
    endif

    call rng_init

    write (ilog, '(A)')  "- sliding_t_init -"

    write (ilog, '(5X,A,F20.10)')  "T           = ", this%T
    write (ilog, '(5X,A,F20.10)')  "dT          = ", this%dT
    write (ilog, '(5X,A,F20.10)')  "dissipation = ", this%dissipation
    write (ilog, '(5X,A,F20.10)')  " -> tau     = ", 1.0_DP/this%dissipation
    write (ilog, '(5X,A,I3)')      "bot         = ", this%bot
    write (ilog, '(5X,A,I3)')      "top         = ", this%top

  endsubroutine sliding_t_init


  !>
  !! Destructor
  !<
  subroutine sliding_t_del(this)
    implicit none

    type(sliding_t_t), intent(inout)  :: this

    ! ---

  endsubroutine sliding_t_del


  !>
  !! Position update and velocity estimation
  !<
  subroutine sliding_t_step1(this, barostat, p, v, f, dt, max_dt, max_dr, &
      max_dr_sq, error)
    implicit none

    type(sliding_t_t),  intent(inout) :: this
    type(sliding_p_t),  intent(in)    :: barostat(:)
    type(particles_t),  intent(inout) :: p
    real(DP),           intent(inout) :: v(3, p%maxnatloc)
    real(DP),           intent(in)    :: f(3, p%maxnatloc)
    real(DP),           intent(inout) :: dt
    real(DP), optional, intent(in)    :: max_dr
    real(DP), optional, intent(in)    :: max_dt
    real(DP), optional, intent(inout) :: max_dr_sq
    integer,  optional, intent(out)   :: error

    ! ---

    real(DP)  :: c0, c1, c2, gamdt, d2t

    integer   :: i, j

    real(DP)  :: etar, etav, sigmar, sigmav, covrv
    real(DP)  :: dr(3), hlp, T_au, v0(3)

    real(DP)  :: l_max_dr_sq

    ! ---

    INIT_ERROR(error)

    if (size(barostat) /= 1) then
       RAISE_ERROR('SlidingT needs a single instance of SlidingP.', error)
    endif

    call timer_start("sliding_t_step1")

    T_au  = this%T * K_to_energy

    !
    ! Adaptive time stepping
    !

    call timestep(p, v, f, dt, max_dt, max_dr)
    
    d2t = dt**2
 

    !
    ! Integrate
    !

    gamdt = this%dissipation*dt
    c0    = exp(-gamdt)
    c1    = (1.0-c0)/gamdt
    c2    = (1.0-c1)/gamdt

    hlp   = 2.d0-(3.d0+c0**2-4.d0*c0)/gamdt

    l_max_dr_sq  = 0.0_DP

    !$omp  parallel do default(none) &
    !$omp& shared(barostat, f, p, this, v) &
    !$omp& firstprivate(c0, c1, c2, dt, d2t, gamdt, hlp, T_au) &
    !$omp& private(covrv, dr, etar, etav, j, sigmar, sigmav, v0) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0) then

          if (p%g(i) == this%bot.or. p%g(i) == this%top) then

             v0  = 0.0_DP
             if (p%g(i) == this%top) then
                ! Thermalize relativ to the top velocity
                v0  = barostat(1)%v_top
             endif

             dr          = v0 * dt + c1 * ( VEC3(v, i) - v0 ) * dt + c2 * VEC3(f, i) / p%m(i) * d2t
             VEC3(v, i)  = v0 + c0 * ( VEC3(v, i) - v0 ) + (c1-c2) * VEC3(f, i) / p%m(i) * dt

             !
             ! The random part (Langevin)
             !
          
             if (hlp > 0.0_DP) then
                sigmar = sqrt(T_au/p%m(i)*d2t/gamdt*hlp)
                sigmav = sqrt(T_au/p%m(i)*(1.d0-c0**2))
                covrv  = T_au/p%m(i)*dt/gamdt*(1.d0-c0)**2

                do j = 1, 3
                   call gaucorr(etar, etav, sigmar, sigmav, covrv)
                   dr(j)         = dr(j) + etar
                   VEC(v, i, j)  = VEC(v, i, j) + etav
                enddo

             endif

          else

             ! Usual Verlet dynamics
             dr          = VEC3(v, i) * dt + 0.5_DP * VEC3(f, i) / p%m(i) * d2t
             VEC3(v, i)  = VEC3(v, i) + 0.5_DP * VEC3(f, i) / p%m(i) * dt

          endif

#ifndef IMPLICIT_R
             POS3(p, i) = POS3(p, i) + dr
#endif
             PNC3(p, i) = PNC3(p, i) + dr
             PCN3(p, i) = PCN3(p, i) + dr

          l_max_dr_sq = max(l_max_dr_sq, dot_product(dr, dr))

       endif

    enddo


    !
    ! Maximum particle displacement
    !

    p%accum_max_dr  = p%accum_max_dr + sqrt(l_max_dr_sq)

    if (present(max_dr_sq)) then
       max_dr_sq  = max(max_dr_sq, l_max_dr_sq)
    endif

    call I_changed_positions(p)

    call timer_stop("sliding_t_step1")

  endsubroutine sliding_t_step1


  !>
  !! Velocity correction
  !<
  subroutine sliding_t_step2(this, p, v, f, dt)
    implicit none

    type(sliding_t_t), intent(inout)  :: this
    type(particles_t), intent(inout)  :: p
    real(DP), intent(inout)           :: v(3, p%maxnatloc)
    real(DP), intent(in)              :: f(3, p%maxnatloc)
    real(DP), intent(in)              :: dt

    ! ---

    real(DP)  :: c0, c1, c2, gamdt

    integer   :: i

    ! ---

    call timer_start("sliding_t_step2")

    !
    ! Communicate forces back to this processor if required
    !

#ifdef _MP
    if (mod_communicator%communicate_forces) then
       DEBUG_WRITE("- communicate_forces -")
       call communicate_forces(mod_communicator, p)
    endif
#endif

    !
    ! Integrate
    !

    gamdt = this%dissipation*dt
    c0 = exp(-gamdt)
    c1 = (1.0-c0)/gamdt
    c2 = (1.0-c1)/gamdt

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, this, v) &
    !$omp& firstprivate(c0, c1, c2, dt)
    do i = 1, p%natloc

       if (p%g(i) > 0) then

          if (p%g(i) == this%bot.or. p%g(i) == this%top) then

             VEC3(v, i)  = VEC3(v, i) + c2 * VEC3(f, i) / p%m(i) * dt

         else

            VEC3(v, i)  = VEC3(v, i) + 0.5_DP * VEC3(f, i) / p%m(i) * dt

         endif

      endif

    enddo

    if (this%dT /= 0.0) then 
       this%T = this%T + dt*this%dT
    endif

    call timer_stop("sliding_t_step2")
    
  endsubroutine sliding_t_step2


  subroutine sliding_t_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(sliding_t_t), target, intent(inout)  :: this
    type(c_ptr),               intent(in)     :: cfg
    type(c_ptr),               intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("SlidingT"), &
         CSTR("Langevin thermostat for the outer borders of the system."))

    call ptrdict_register_integer_property(m, c_loc(this%bot), CSTR("bot"), &
         CSTR("Group of the bottom thermostat atoms."))
    call ptrdict_register_integer_property(m, c_loc(this%top), CSTR("top"), &
         CSTR("Group of the top thermostat atoms."))

    call ptrdict_register_real_property(m, c_loc(this%T), CSTR("T"), &
         CSTR("Temperature for the Langevin thermostat."))
    call ptrdict_register_real_property(m, c_loc(this%dT), CSTR("dT"), &
         CSTR("Temperature change per time step."))
    call ptrdict_register_real_property(m, c_loc(this%dissipation), &
         CSTR("dissipation"), &
         CSTR("Dissipation constant for the Langevin thermostat."))
    call ptrdict_register_real_property(m, c_loc(this%tau), CSTR("tau"), &
         CSTR("Relaxation time constant for the Langevin thermostat."))
    
  endsubroutine sliding_t_register

endmodule sliding_t
