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
!   classtype:verlet_langevin_t classname:Langevin interface:integrators
! @endmeta

!>
!! Velocity-verlet integrator with homogeneous Langevin dynamics for the whole
!! system
!!
!! Velocity-verlet integrator with homogeneous Langevin dynamics for the whole
!! system
!<

#include "macros.inc"
#include "filter.inc"

module verlet_langevin
  use supplib
  use rng

  use particles

  use filter
  use verlet_support

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  public :: verlet_langevin_t
  type verlet_langevin_t

     character(MAX_EL_STR)  :: elements     = "*"
     integer                :: els          = 0

     real(DP)               :: T            = -1.0_DP
     real(DP)               :: dT           = 0.0_DP

     real(DP)               :: dissipation  = 1.0_DP
     real(DP)               :: tau          = -1.0_DP

  endtype verlet_langevin_t


  public :: init
  interface init
     module procedure verlet_global_langevin_init
  endinterface

  public :: del
  interface del
     module procedure verlet_global_langevin_del
  endinterface

  public :: step1
  interface step1
     module procedure verlet_global_langevin_step1
  endinterface

  public :: step2
  interface step2
     module procedure verlet_global_langevin_step2
  endinterface

  public :: register
  interface register
    module procedure verlet_langevin_register
  endinterface

contains

  !>
  !! Global Langevin init
  !!
  !! Global Langevin init. Either dissipation (=1/tau) or tau should be
  !! specified, not both. If nothing is given, the default value of dissipation
  !! is used. If this%tau is set, then that is used.
  !<
  subroutine verlet_global_langevin_init(this, p, T, dT, dissipation, tau, error)
    implicit none

    type(verlet_langevin_t), intent(inout)  :: this
    type(particles_t), intent(inout)        :: p
    real(DP), intent(in), optional          :: T
    real(DP), intent(in), optional          :: dT
    real(DP), intent(in), optional          :: dissipation
    real(DP), intent(in), optional          :: tau
    integer, intent(out), optional          :: error
        
    ! ---

    INIT_ERROR(error)

    call prlog("- verlet_global_langevin_init -")

    ! Checks
    if (present(dissipation) .and. present(tau)) then
       RAISE_ERROR("Please specify either *dissipation* or *tau*.", error)
    end if

    ! Init
    if (present(T)) then
       this%T = T
    endif
    if (present(dT)) then
       this%dT = dT
    endif

    call prlog("     T           = "//this%T)
    call prlog("     dT          = "//this%dT)

    if (present(dissipation)) then
       this%dissipation  = dissipation
       this%tau = -1.0_DP
       if (present(tau)) then
          RAISE_ERROR("Please specify either *dissipation* or *tau*.", error)
       endif
    endif
    if (present(tau)) then
       this%tau = tau
    endif

    if (this%tau > 0.0_DP) then
       this%dissipation  = 1.0_DP/this%tau
       call prlog("tau           = "//1.0_DP/this%dissipation)
       call prlog("* dissipation = "//this%dissipation)
    else
       call prlog("dissipation   = "//this%dissipation)
       call prlog("* tau         = "//1.0_DP/this%dissipation)
    endif

    call rng_init

    call prlog

  endsubroutine verlet_global_langevin_init


  !**********************************************************************
  ! Delete a Verlet object
  !**********************************************************************
  subroutine verlet_global_langevin_del(this)
    implicit none

    type(verlet_langevin_t), intent(inout)  :: this

    ! ---

  endsubroutine verlet_global_langevin_del


  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine verlet_global_langevin_step1(this, p, v, f, dt, max_dt, max_dr, max_dr_sq)
    implicit none

    type(verlet_langevin_t), intent(inout)  :: this
    type(particles_t), intent(inout)               :: p
    real(DP), intent(inout)                        :: v(3, p%maxnatloc)
    real(DP), intent(in)                           :: f(3, p%maxnatloc)
    real(DP), intent(inout)                        :: dt
    real(DP), intent(in), optional                 :: max_dr
    real(DP), intent(in), optional                 :: max_dt
    real(DP), intent(inout), optional              :: max_dr_sq

    ! ---

    real(DP)  :: c0, c1, c2, gamdt, d2t

    integer   :: i, j

    real(DP)  :: etar, etav, sigmar, sigmav, covrv
    real(DP)  :: dr(3), hlp, T_au

    real(DP)  :: l_max_dr_sq

    ! ---

    call timer_start("verlet_global_langevin_step1")

    if (this%els == 0) then
       this%els  = filter_from_string(this%elements, p)
    endif

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
    !$omp& shared(f, p, this, v) &
    !$omp& firstprivate(c0, c1, c2, dt, d2t, gamdt, hlp, T_au) &
    !$omp& private(covrv, dr, etar, etav, j, sigmar, sigmav) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) then

          dr          = c1 * VEC3(v, i) * dt + c2 * VEC3(f, i) / p%m(i) * d2t
          VEC3(v, i)  = c0 * VEC3(v, i) + (c1-c2) * VEC3(f, i) / p%m(i) * dt

          !
          ! The random part (Langevin)
          !
          
          if (hlp > 0.0) then
             sigmar = sqrt(T_au/p%m(i)*d2t/gamdt*hlp)
             sigmav = sqrt(T_au/p%m(i)*(1.d0-c0**2))
             covrv  = T_au/p%m(i)*dt/gamdt*(1.d0-c0)**2

             do j = 1, 3
                call gaucorr(etar, etav, sigmar, sigmav, covrv)
                dr(j)         = dr(j) + etar
                VEC(v, i, j)  = VEC(v, i, j) + etav
             enddo

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

    call timer_stop("verlet_global_langevin_step1")

  endsubroutine verlet_global_langevin_step1


  !**********************************************************************
  ! Velocity correction
  !**********************************************************************
  subroutine verlet_global_langevin_step2(this, p, v, f, dt)
    implicit none

    type(verlet_langevin_t), intent(inout)  :: this
    type(particles_t), intent(inout)               :: p
    real(DP), intent(inout)                        :: v(3, p%maxnatloc)
    real(DP), intent(in)                           :: f(3, p%maxnatloc)
    real(DP), intent(in)                           :: dt

    ! ---

    real(DP)  :: c0, c1, c2, gamdt

    integer   :: i

    ! ---

    call timer_start("verlet_global_langevin_step2")

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

       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) &
            VEC3(v, i)  = VEC3(v, i) + c2 * VEC3(f, i) / p%m(i) * dt

    enddo

    if (this%dT /= 0.0) then 
       this%T = this%T + dt*this%dT
    endif

    !
    ! Update virial and kinetic energy
    !

!    call compute_kinetic_energy_and_virial(p)

    call timer_stop("verlet_global_langevin_step2")
    
  endsubroutine verlet_global_langevin_step2


  subroutine verlet_langevin_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(verlet_langevin_t), target, intent(inout)  :: this
    type(c_ptr),                     intent(in)     :: cfg
    type(c_ptr),                     intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Langevin"), &
         CSTR("Global Langevin thermostat."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements for which to enable this integrator."))

    call ptrdict_register_real_property(m, c_loc(this%T), CSTR("T"), &
         CSTR("Temperature for the Langevin thermostat."))
!    call ptrdict_register_real_property(m, this%Tend, CSTR("Tend"), &
!         CSTR("End temperature for thermostat (if > 0, exponential cooling is used)."))
    call ptrdict_register_real_property(m, c_loc(this%dT), CSTR("dT"), &
         CSTR("Temperature change per time step."))
    call ptrdict_register_real_property(m, c_loc(this%dissipation), &
         CSTR("dissipation"), &
         CSTR("Dissipation constant for the Langevin thermostat."))
    call ptrdict_register_real_property(m, c_loc(this%tau), CSTR("tau"), &
         CSTR("Relaxation time constant for the Langevin thermostat."))
    
  endsubroutine verlet_langevin_register

endmodule verlet_langevin
