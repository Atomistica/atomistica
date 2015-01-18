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
!   classtype:verlet_t classname:Verlet interface:integrators
! @endmeta

!>
!! Standard Velocity-Verlet integrator
!!
!! Standard Velocity-Verlet integrator
!<

#include "macros.inc"
#include "filter.inc"

module verlet
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

  public :: verlet_t
  type verlet_t

     character(MAX_EL_STR)  :: elements  = "*"
     integer                :: els       = 0

  endtype verlet_t


  public :: del
  interface del
     module procedure verlet_del
  endinterface

  public :: step1
  interface step1
     module procedure verlet_step1
  endinterface

  public :: step2
  interface step2
     module procedure verlet_step2
  endinterface

  public :: register
  interface register
    module procedure verlet_register
  endinterface

contains

  !>
  !! Destructor
  !<
  subroutine verlet_del(this)
    implicit none

    type(verlet_t), intent(inout)  :: this

    ! ---

  endsubroutine verlet_del


  !>
  !! Position update and velocity estimation
  !<
  subroutine verlet_step1(this, p, v, f, dt, max_dt, max_dr, max_dr_sq)
    implicit none

    type(verlet_t), intent(inout)     :: this
    type(particles_t), intent(inout)  :: p
    real(DP), intent(inout)           :: v(3, p%maxnatloc)
    real(DP), intent(in)              :: f(3, p%maxnatloc)
    real(DP), intent(inout)           :: dt
    real(DP), intent(in), optional    :: max_dt
    real(DP), intent(in), optional    :: max_dr
    real(DP), intent(inout), optional :: max_dr_sq

    ! ---

    integer            :: i

    real(DP)           :: dr(3), l_max_dr_sq

    ! ---

    call timer_start("verlet_step1")

    if (this%els == 0) then
       this%els  = filter_from_string(this%elements, p)
    endif

    !
    ! Adaptive time stepping
    !

    call timestep(p, v, f, dt, max_dt, max_dr)
        
    !
    ! Integrate
    !

    l_max_dr_sq  = 0.0_DP

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, this, v) &
    !$omp& firstprivate(dt) &
    !$omp& private(dr) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) then

          VEC3(v, i)       = VEC3(v, i) + 0.5_DP * VEC3(f, i) / p%m(i) * dt
          dr               = VEC3(v, i) * dt

#ifndef IMPLICIT_R
          POS3(p, i)       = POS3(p, i) + dr
#endif
          PNC3(p, i)       = PNC3(p, i) + dr
          PCN3(p, i)       = PCN3(p, i) + dr

          l_max_dr_sq      = max(l_max_dr_sq, dot_product(dr, dr))

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

    call timer_stop("verlet_step1")

  endsubroutine verlet_step1


  !>
  !! Velocity correction
  !<
  subroutine verlet_step2(this, p, v, f, dt)
    implicit none

    type(verlet_t), intent(in)        :: this
    type(particles_t), intent(inout)  :: p
    real(DP), intent(inout)           :: v(3, p%maxnatloc)
    real(DP), intent(in)              :: f(3, p%maxnatloc)
    real(DP), intent(in)              :: dt

    ! ---

    integer   :: i

    ! ---

    call timer_start("verlet_step2")

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

    !$omp parallel do default(none) &
    !$omp shared(f, p, this, v) &
    !$omp firstprivate(dt)
    do i = 1, p%natloc

       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) &
            VEC3(v, i) = VEC3(v, i) + 0.5_DP * VEC3(f, i) / p%m(i) * dt

    enddo


    !
    ! Update virial and kinetic energy
    !

!    call compute_kinetic_energy_and_virial(p)

    call timer_stop("verlet_step2")
    
  endsubroutine verlet_step2


  subroutine verlet_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(verlet_t), target, intent(inout)  :: this
    type(c_ptr),            intent(in)     :: cfg
    type(c_ptr),            intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Verlet"), &
         CSTR("The Velocity-Verlet integrator."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, CSTR("elements"), &
         CSTR("Elements for which to enable this integrator."))

  endsubroutine verlet_register

endmodule verlet
