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
!   classtype:verlet_local_langevin_t classname:LocalLangevin
!   interface:integrators
! @endmeta

!>
!! Velocity-verlet integrator with a Langevin thermostat that allows
!! specification of a per-atom temperature and dissipation constant.
!!
!! Velocity-verlet integrator with a Langevin thermostat that allows
!! specification of a per-atom temperature and dissipation constant.
!<

#include "macros.inc"

module verlet_local_langevin
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

  public :: verlet_local_langevin_t
  type verlet_local_langevin_t

     real(DP), pointer  :: T(:)
     real(DP), pointer  :: dissipation(:)

  endtype verlet_local_langevin_t


  public :: init
  interface init
     module procedure verlet_local_langevin_init
  endinterface

  public :: del
  interface del
     module procedure verlet_local_langevin_del
  endinterface

  public :: step1
  interface step1
     module procedure verlet_local_langevin_step1
  endinterface

  public :: step2
  interface step2
     module procedure verlet_local_langevin_step2
  endinterface

  public :: register
  interface register
    module procedure verlet_local_langevin_register
  endinterface

contains

  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine verlet_local_langevin_init(this, p)
    implicit none

    type(verlet_local_langevin_t), intent(inout)  :: this
    type(particles_t), intent(inout)              :: p

    ! ---

    if (.not. rng_initialized) then
       call rng_init
    endif

    call add_real( &
         p%data, &
         T_STR, &
         F_VERBOSE_ONLY + F_COMMUNICATE, &
         "Kelvins" )
    call add_real( &
         p%data, &
         DISSIPATION_STR, &
         F_CONSTANT + F_VERBOSE_ONLY + F_COMMUNICATE )

    this%T            => NULL()
    this%dissipation  => NULL()

  endsubroutine verlet_local_langevin_init


  !**********************************************************************
  ! Called if the particles object changes
  !**********************************************************************
  subroutine verlet_local_langevin_del(this)
    implicit none

    type(verlet_local_langevin_t), intent(inout)  :: this

    ! ---

  endsubroutine verlet_local_langevin_del


  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine verlet_local_langevin_step1(this, p, v, f,dt, max_dt, max_dr, max_dr_sq)
    implicit none

    type(verlet_local_langevin_t), intent(inout)  :: this
    type(particles_t), intent(inout)              :: p
    real(DP), intent(inout)                       :: v(3, p%maxnatloc)
    real(DP), intent(in)                          :: f(3, p%maxnatloc)
    real(DP), intent(inout)                       :: dt
    real(DP), intent(in), optional                :: max_dt
    real(DP), intent(in), optional                :: max_dr
    real(DP), intent(inout), optional             :: max_dr_sq

    ! ---

    real(DP)           :: c0, c1, c2, gamdt, d2t

    integer            :: i, j

    real(DP)           :: etar, etav, sigmar, sigmav, covrv, cur_T
    real(DP)           :: dr(3), hlp, l_max_dr_sq

    ! ---

    call timer_start("verlet_local_langevin_step1")

    if (.not. associated(this%T)) &
         call ptr_by_name(p%data, T_STR, this%T)
    if (.not. associated(this%dissipation)) &
         call ptr_by_name(p%data, DISSIPATION_STR, this%dissipation)

    !
    ! Adaptive time stepping
    !

    call timestep(p, v, f, dt, max_dt, max_dr)

    d2t = dt**2

    !
    ! Integrate
    !

    l_max_dr_sq  = 0.0_DP

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, this, v) &
    !$omp& firstprivate(dt, d2t, K_to_energy) &
    !$omp& private(c0, c1, c2, covrv, cur_T, dr, etar, etav, gamdt, hlp, i, sigmar, sigmav) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0) then

          if (this%dissipation(i) > 0.0_DP) then
             gamdt = this%dissipation(i)*dt
             c0 = exp(-gamdt)
             c1 = (1.0_DP-c0)/gamdt
             c2 = (1.0_DP-c1)/gamdt
          else
             c0 = 1.0_DP
             c1 = 1.0_DP
             c2 = 0.5_DP
          endif

          dr          = c1 * VEC3(v, i) * dt + c2 * VEC3(f, i) / p%m(i) * d2t
          VEC3(v, i)  = c0 * VEC3(v, i) + (c1-c2) * VEC3(f, i) / p%m(i) * dt

          !
          ! The random part (Langevin)
          !
          
          if (this%dissipation(i) > 0.0_DP) then
             cur_T = this%T(i)*K_to_energy

             hlp = 2.d0-(3.d0+c0**2-4.d0*c0)/gamdt

             if (hlp > 0.0_DP .and. cur_T > 0.0_DP) then
                ! Only noise if T > 0, otherwise this will only dampen
                sigmar = sqrt(cur_T/p%m(i)*d2t/gamdt*hlp)
                sigmav = sqrt(cur_T/p%m(i)*(1.d0-c0**2))
                covrv  = cur_T/p%m(i)*dt/gamdt*(1.d0-c0)**2

                do j = 1, 3
                   call gaucorr(etar, etav, sigmar, sigmav, covrv)
                   dr(j)              = dr(j) + etar
                   VEC(v, i, j)  = VEC(v, i, j) + etav
                enddo
             endif
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

    call timer_stop("verlet_local_langevin_step1")

  endsubroutine verlet_local_langevin_step1


  !**********************************************************************
  ! Velocity correction
  !**********************************************************************
  subroutine verlet_local_langevin_step2(this, p, v, f, dt)
    implicit none

    type(verlet_local_langevin_t), intent(in)  :: this
    type(particles_t), intent(inout)           :: p
    real(DP), intent(inout)                    :: v(3, p%maxnatloc)
    real(DP), intent(in)                       :: f(3, p%maxnatloc)
    real(DP), intent(in)                       :: dt

    ! ---

    real(DP)  :: c0, c1, c2, gamdt

    integer   :: i

    ! ---

    call timer_start("verlet_local_langevin_step2")

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

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, this, v) &
    !$omp& firstprivate(dt) &
    !$omp& private(c0, c1, c2, gamdt)
    do i = 1, p%natloc
       
       if (p%g(i) > 0) then

          if (this%dissipation(i) > 0.0_DP) then
             gamdt = this%dissipation(i)*dt
             c0 = exp(-gamdt)
             c1 = (1.0_DP-c0)/gamdt
             c2 = (1.0_DP-c1)/gamdt
          else
             c2 = 0.5_DP
          endif

          VEC3(v, i) = VEC3(v, i) + c2 * VEC3(f, i) / p%m(i) * dt

       endif

    enddo

    !
    ! Update virial and kinetic energy
    !

!    call compute_kinetic_energy_and_virial(p)

    call timer_stop("verlet_local_langevin_step2")
    
  endsubroutine verlet_local_langevin_step2


  subroutine verlet_local_langevin_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(verlet_local_langevin_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                       :: cfg
    type(c_ptr), intent(out)                      :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("LocalLangevin"), &
         CSTR("Local Langevin thermostat."))

  endsubroutine verlet_local_langevin_register

endmodule verlet_local_langevin
