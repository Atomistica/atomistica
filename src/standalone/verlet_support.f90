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
! @endmeta

!>
!! Helper subroutines used by all Velocity-Verlet integrators
!!
!! Helper subroutines used by all Velocity-Verlet integrators
!<

#include "filter.inc"
#include "macros.inc"

module verlet_support
  use supplib

  use data
  use particles

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  character(MAX_NAME_STR), parameter  :: T_STR            = "temperatures"
  character(MAX_NAME_STR), parameter  :: DISSIPATION_STR  = "langevin_dissipation"

  public :: T_STR, DISSIPATION_STR
  public :: timestep, verlet_r, verlet_v

contains

  !>
  !! Adapt time step
  !<
  subroutine timestep(p, v, f, dt, max_dt, max_dr, error)
    implicit none

    type(particles_t),  intent(in)    :: p
    real(DP),           intent(in)    :: v(3, p%maxnatloc)
    real(DP),           intent(in)    :: f(3, p%maxnatloc)
    real(DP),           intent(inout) :: dt
    real(DP), optional, intent(in)    :: max_dt
    real(DP), optional, intent(in)    :: max_dr
    integer,  optional, intent(out)   :: error

    ! ---

    integer   :: i
    
    real(DP)  :: dr(3), dr_sq, max_dr_sq, d2t

    ! ---

    INIT_ERROR(error)

    if (present(max_dr) .and. present(max_dr)) then

       if (max_dr > 0.0_DP) then

          d2t = max_dt**2
    
          max_dr_sq = 0.0
          !$omp  parallel do default(none) &
          !$omp& shared(d2t, f, max_dt, p, v) &
          !$omp& private(dr, dr_sq) reduction(max:max_dr_sq)
          do i = 1, p%natloc
             if (p%g(i) > 0) then
                dr         = VEC3(v, i)*max_dt + sqrt(dot_product(VEC3(f, i), VEC3(f, i)))*VEC3(f, i)/p%m(i)*d2t
                dr_sq      = dot_product(dr, dr)
                max_dr_sq  = max(dr_sq, max_dr_sq)
             endif
          enddo
          
          if (max_dr_sq > 0.0) then
             dt = min(max_dr/sqrt(max_dr_sq)*max_dt, max_dt)
          else
             dt = max_dt
          endif

#ifdef _MP
          dt = min(mod_communicator%mpi, dt, error)
          PASS_ERROR(error)
#endif

       endif

    endif

  endsubroutine timestep


  !>
  !! Update velocities (Verlet half-step)
  !<
  subroutine verlet_v(els, p, v, f, dt)
    implicit none

    integer,           intent(in)    :: els
    type(particles_t), intent(in)    :: p
    real(DP),          intent(inout) :: v(3, p%maxnatloc)
    real(DP),          intent(in)    :: f(3, p%maxnatloc)
    real(DP),          intent(in)    :: dt

    ! ---

    integer  :: i

    ! ---

    !$omp  parallel do default(none) &
    !$omp& shared(els, f, p, v) &
    !$omp& firstprivate(dt)
    do i = 1, p%natloc

       if (p%g(i) > 0 .and. IS_EL(els, p, i)) then
          VEC3(v, i) = VEC3(v, i) + 0.5_DP * VEC3(f, i) / p%m(i) * dt
       endif

    enddo

  endsubroutine verlet_v


  !>
  !! Update positions (Verlet half-step)
  !<
  subroutine verlet_r(els, p, v, f, dt, l_max_dr_sq, fac)
    implicit none

    integer,            intent(in)    :: els
    type(particles_t),  intent(inout) :: p
    real(DP),           intent(in)    :: v(3, p%maxnatloc)
    real(DP),           intent(in)    :: f(3, p%maxnatloc)
    real(DP),           intent(in)    :: dt
    real(DP),           intent(out)   :: l_max_dr_sq
    real(DP), optional, intent(in)    :: fac(3)

    ! ---

    integer  :: i
    real(DP) :: dr(3), vfac(3)

    ! ---

    l_max_dr_sq  = 0.0_DP

    vfac = 1.0_DP
    if (present(fac)) then
       vfac = fac
    endif

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, v) &
    !$omp& firstprivate(dt, els, vfac) &
    !$omp& private(dr) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0 .and. IS_EL(els, p, i)) then

          dr               = vfac * VEC3(v, i) * dt

#ifndef IMPLICIT_R
          POS3(p, i)       = POS3(p, i) + dr
#endif
          PNC3(p, i)       = PNC3(p, i) + dr
          PCN3(p, i)       = PCN3(p, i) + dr

          l_max_dr_sq      = max(l_max_dr_sq, dot_product(dr, dr))

       endif

    enddo

  endsubroutine verlet_r

endmodule verlet_support
