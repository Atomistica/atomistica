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
!   classtype:verlet_langevin_1d_t classname:Langevin1D interface:integrators
! @endmeta

!>
!! Velocity-verlet integrator with homogeneous Langevin dynamics for the whole
!! system, but thermalization only in a single of the three cartesian
!! directions.
!!
!! Velocity-verlet integrator with homogeneous Langevin dynamics for the whole
!! system, but thermalization only in a single of the three cartesian
!! directions.
!<

#include "macros.inc"

module verlet_langevin_1d
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

  integer, parameter  :: n_dims       = 3
  integer, parameter  :: len_dim_str  = 15
  integer, parameter  :: ALL_DIMS     = 0

  ! This is need for xlf
  character(len_dim_str), parameter  :: STR_x            = CSTR("x")
  character(len_dim_str), parameter  :: STR_y            = CSTR("y")
  character(len_dim_str), parameter  :: STR_z            = CSTR("z")
  character(len_dim_str), parameter  :: dim_strs(n_dims) = &
       (/ STR_x, STR_y, STR_z /)

  public :: verlet_langevin_1d_t
  type verlet_langevin_1d_t

     integer   :: d            = 2
     integer   :: d2           = 2
     integer   :: d3           = 2

     real(DP)  :: T            = -1.0_DP
     real(DP)  :: dT           = 0.0_DP

     real(DP)  :: dissipation  = 1.0_DP
     real(DP)  :: tau          = -1.0_DP

  endtype verlet_langevin_1d_t


  public :: init
  interface init
     module procedure verlet_global_langevin_1d_init
  endinterface

  public :: del
  interface del
     module procedure verlet_global_langevin_1d_del
  endinterface

  public :: step1
  interface step1
     module procedure verlet_global_langevin_1d_step1
  endinterface

  public :: step2
  interface step2
     module procedure verlet_global_langevin_1d_step2
  endinterface

  public :: register
  interface register
    module procedure verlet_global_langevin_1d_register
  endinterface

contains

  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine verlet_global_langevin_1d_init(this, p, d, T, dT, dissipation, &
       tau, ierror)
    implicit none

    type(verlet_langevin_1d_t), intent(inout)  :: this
    type(particles_t),          intent(inout)  :: p
    integer,          optional, intent(in)     :: d
    real(DP),         optional, intent(in)     :: T
    real(DP),         optional, intent(in)     :: dT
    real(DP),         optional, intent(in)     :: dissipation
    real(DP),         optional, intent(in)     :: tau
    integer,          optional, intent(inout)  :: ierror

    ! ---

    INIT_ERROR(ierror)

    call prlog("- verlet_global_langevin_1d_init -")

    if (present(d)) then
       this%d = d
    endif
    if (present(T)) then
       this%T = T
    endif
    if (present(dT)) then
       this%dT = dT
    endif

    this%d2 = modulo(this%d+1, 3)+1
    this%d3 = modulo(this%d+2, 3)+1
    this%d  = modulo(this%d,   3)+1

    call prlog("d             = ( "//this%d//", "//this%d2//", "//this%d3//" )")
    call prlog("T             = "//this%T)
    call prlog("dT            = "//this%dT)

    if (present(dissipation)) then
       this%dissipation  = dissipation
       this%tau = -1.0_DP
       if (present(tau)) then
          RAISE_ERROR("Please specify either *dissipation* or *tau*.", ierror)
       endif
    endif
    if (present(tau)) then
       this%tau  = tau
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

  endsubroutine verlet_global_langevin_1d_init


  !**********************************************************************
  ! Delete a Verlet object
  !**********************************************************************
  subroutine verlet_global_langevin_1d_del(this)
    implicit none

    type(verlet_langevin_1d_t), intent(inout)  :: this

    ! ---

  endsubroutine verlet_global_langevin_1d_del


  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine verlet_global_langevin_1d_step1(this, p, v, f, dt, max_dt, &
       max_dr, max_dr_sq)
    implicit none

    type(verlet_langevin_1d_t), intent(inout)  :: this
    type(particles_t), intent(inout)                  :: p
    real(DP), intent(inout)                           :: v(3, p%maxnatloc)
    real(DP), intent(in)                              :: f(3, p%maxnatloc)
    real(DP), intent(inout)                           :: dt
    real(DP), intent(in), optional                    :: max_dr
    real(DP), intent(in), optional                    :: max_dt
    real(DP), intent(inout), optional                 :: max_dr_sq

    ! ---

    real(DP)           :: c0, c1, c2, gamdt, d2t

    integer            :: i

    real(DP)           :: etar, etav, sigmar, sigmav, covrv
    real(DP)           :: dr(3), hlp, l_max_dr_sq, T_au

    ! ---

    call timer_start("verlet_global_langevin_1d_step1")

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
    !$omp& private(covrv, dr, etar, etav, sigmar, sigmav) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0) then

          dr(this%d)               = c1 * VEC(v, i, this%d) * dt + c2 * VEC(f, i, this%d) / p%m(i) * d2t
          dr(this%d2)              = VEC(v, i, this%d2) * dt + 0.5 * VEC(f, i, this%d2) / p%m(i) * d2t
          dr(this%d3)              = VEC(v, i, this%d3) * dt + 0.5 * VEC(f, i, this%d3) / p%m(i) * d2t

          VEC(v, i, this%d)   = c0*VEC(v, i, this%d) + (c1-c2) * VEC(f, i, this%d) / p%m(i) * dt
          VEC(v, i, this%d2)  = VEC(v, i, this%d2) + 0.5 * VEC(f, i, this%d2) / p%m(i) * dt
          VEC(v, i, this%d3)  = VEC(v, i, this%d3) + 0.5 * VEC(f, i, this%d3) / p%m(i) * dt

          !
          ! The random part (Langevin)
          !
          
          if (hlp > 0.0) then
             sigmar = sqrt(T_au/p%m(i)*d2t/gamdt*hlp)
             sigmav = sqrt(T_au/p%m(i)*(1.d0-c0**2))
             covrv  = T_au/p%m(i)*dt/gamdt*(1.d0-c0)**2

             call gaucorr(etar, etav, sigmar, sigmav, covrv)
             dr(this%d)             = dr(this%d) + etar
             VEC(v, i, this%d) = VEC(v, i, this%d) + etav
          endif

#ifndef IMPLICIT_R
          POS3(p, i) = POS3(p, i) + dr(:)
#endif
          PNC3(p, i) = PNC3(p, i) + dr(:)
          PCN3(p, i) = PCN3(p, i) + dr(:)

          l_max_dr_sq = max(l_max_dr_sq, dot_product(dr(:), dr(:)))

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

    call timer_stop("verlet_global_langevin_1d_step1")

  endsubroutine verlet_global_langevin_1d_step1


  !**********************************************************************
  ! Velocity correction
  !**********************************************************************
  subroutine verlet_global_langevin_1d_step2(this, p, v, f, dt)
    implicit none

    type(verlet_langevin_1d_t), intent(inout)  :: this
    type(particles_t), intent(inout)                  :: p
    real(DP), intent(inout)                           :: v(3, p%maxnatloc)
    real(DP), intent(in)                              :: f(3, p%maxnatloc)
    real(DP), intent(in)                              :: dt

    ! ---

    real(DP)  :: c0, c1, c2, gamdt

    integer   :: i

    ! ---

    call timer_start("verlet_global_langevin_1d_step2")

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
    !$omp& shared(f, c0, c1, c2, dt, p, this, v)
    do i = 1, p%natloc

       if (p%g(i) > 0) then
          VEC(v, i, this%d)   = VEC(v, i, this%d) + c2 * VEC(f, i, this%d) / p%m(i) * dt
          VEC(v, i, this%d2)  = VEC(v, i, this%d2) + 0.5 * VEC(f, i, this%d2) / p%m(i) * dt
          VEC(v, i, this%d3)  = VEC(v, i, this%d3) + 0.5 * VEC(f, i, this%d3) / p%m(i) * dt
       endif

    enddo

    if (this%dT /= 0.0) then 
       this%T = this%T + this%dT
    endif

    !
    ! Update virial and kinetic energy
    !

!    call compute_kinetic_energy_and_virial(p)

    call timer_stop("verlet_global_langevin_1d_step2")
    
  endsubroutine verlet_global_langevin_1d_step2


  subroutine verlet_global_langevin_1d_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(verlet_langevin_1d_t), target, intent(inout)  :: this
    type(c_ptr),                        intent(in)     :: cfg
    type(c_ptr),                        intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Langevin1D"), &
         CSTR("Global Langevin thermostat with thermalization only in one direction."))

    call ptrdict_register_enum_property(m, c_loc(this%d), &
         n_dims, len_dim_str, dim_strs(:), &
         CSTR("d"), &
         CSTR("Dimension to thermalize: 'x', 'y' or 'z'"))

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
    
  endsubroutine verlet_global_langevin_1d_register

endmodule verlet_langevin_1d
