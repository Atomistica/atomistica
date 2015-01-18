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
!   classtype:harmonic_hook_t classname:HarmonicHook interface:potentials
! @endmeta

!>
!! Harmonic hook
!!
!! Anchor a set of particles to a certain position in space.
!<

#include "macros.inc"
#include "filter.inc"

module harmonic_hook
  use, intrinsic :: iso_c_binding

  use libAtoms_module

  use logging

  use particles
  use neighbors
  use filter

  implicit none

  private

  public :: POINT_HOOK, LINE_HOOK, PLANE_HOOK
  public :: n_hook_type, len_hook_type_str
  public :: hook_type_strs, hook_type_cstrs

  integer, parameter  :: POINT_HOOK  = 0
  integer, parameter  :: LINE_HOOK   = 1
  integer, parameter  :: PLANE_HOOK  = 2

  integer, parameter  :: n_hook_type        = 3
  integer, parameter  :: len_hook_type_str  = 6
    
  character(*), parameter  :: STR_point  = "point"
  character(*), parameter  :: STR_line   = "line "
  character(*), parameter  :: STR_plane  = "plane"
  character(len_hook_type_str), parameter  :: hook_type_strs(0:n_hook_type-1) = &
       (/ STR_point, STR_line, STR_plane /)
  character(len_hook_type_str), parameter  :: hook_type_cstrs(n_hook_type) = &
       (/ CSTR(STR_point), CSTR(STR_line), CSTR(STR_plane) /)

  ! ---

  public :: harmonic_hook_t
  type harmonic_hook_t

     !
     ! Elements on which to act
     !

     character(MAX_EL_STR)  :: elements  = "*"
     integer                :: els  = 0

     !
     ! Potential parameters
     !

     integer   :: hook_type  = POINT_HOOK                     !< Type of hook

     real(DP)  :: k          = 0.01_DP                        !< Spring constant

     real(DP)  :: r0(3)      = (/ 0.0_DP, 0.0_DP, 0.0_DP /)   !< Anchor
     real(DP)  :: d(3)       = (/ 0.0_DP, 0.0_DP, 1.0_DP /)   !< Anchor direction

  endtype harmonic_hook_t


  public :: init
  interface init
     module procedure harmonic_hook_init
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure harmonic_hook_energy_and_forces
  endinterface

  public :: register
  interface register
    module procedure harmonic_hook_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine harmonic_hook_init(this)
    implicit none

    type(harmonic_hook_t), intent(inout)  :: this

    ! ---

    call prlog("- harmonic_hook_init -")

    this%els  = 0

    call prlog("     " // this%elements)
    call prlog("     hook_type  = " // hook_type_strs(this%hook_type))
    call prlog("     r0         = " // this%r0)

    if (this%hook_type == LINE_HOOK .or. this%hook_type == PLANE_HOOK) then
       this%d  = this%d / sqrt(dot_product(this%d, this%d))
       call prlog("     d          = " // this%d)
    endif

    call prlog

  endsubroutine harmonic_hook_init


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine harmonic_hook_energy_and_forces(this, p, nl, epot, for, wpot, epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)
    implicit none

    type(harmonic_hook_t), intent(inout)  :: this
    type(particles_t), intent(inout)      :: p
    type(neighbors_t), intent(in)         :: nl
    real(DP), intent(inout)               :: epot
    real(DP), intent(inout)               :: for(3, p%maxnatloc)
    real(DP), intent(inout)               :: wpot(3, 3)
    real(DP), intent(inout), optional     :: epot_per_at(p%maxnatloc)
    real(DP), intent(inout), optional     :: epot_per_bond(nl%neighbors_size)
    real(DP), intent(inout), optional     :: f_per_bond(3, nl%neighbors_size)
    real(DP), intent(inout), optional     :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), intent(inout), optional     :: wpot_per_bond(3, 3, nl%neighbors_size)
    integer, intent(inout), optional      :: ierror

    ! ---

    integer   :: i

    real(DP)  :: m, r0(3), f(3)

    ! ---

    if (this%els == 0) then
       this%els  = filter_from_string(this%elements, p)
    endif

    r0  = 0.0_DP
    m   = 0.0_DP
    !$omp  parallel do default(none) &
    !$omp& shared(p, this) &
    !$omp& reduction(+:r0) reduction(+:m)
    do i = 1, p%natloc
       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) then
          r0  = r0 + p%m(i)*POS3(p, i)
          m   = m + p%m(i)
       endif
    enddo
    m   = 1.0_DP/m
    r0  = m*r0

    f   = this%r0 - r0

    select case (this%hook_type)

       case (POINT_HOOK)

       case (LINE_HOOK)
          f  = f - dot_product(f, this%d)*this%d

       case (PLANE_HOOK)
          f  = dot_product(f, this%d)*this%d

       case default
          RAISE_ERROR("Internal error: Unknown hook_type.", ierror)

    endselect

    epot  = epot + 0.5_DP*this%k*dot_product(f, f)
    f     = f * this%k * m

    !$omp  parallel do default(none) &
    !$omp& shared(for, p, this) &
    !$omp& firstprivate(f)
    do i = 1, p%natloc
       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) then
          VEC3(for, i)  = VEC3(for, i) + p%m(i)*f
       endif
    enddo

  endsubroutine harmonic_hook_energy_and_forces


  subroutine harmonic_hook_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(harmonic_hook_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                  :: cfg
    type(c_ptr), intent(out)                 :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("HarmonicHook"), &
         CSTR("Anchor a set of atoms to a spring."))

    call ptrdict_register_string_property(m, c_loc(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("List of elements on which this potential should act."))

    call ptrdict_register_enum_property(m, c_loc(this%hook_type), &
         n_hook_type-1, len_hook_type_str, hook_type_cstrs, &
         CSTR("hook_type"), &
         CSTR("Hook type: 'point', 'line', 'plane'"))

    call ptrdict_register_real_property(m, c_loc(this%k), CSTR("k"), &
         CSTR("Spring constant."))
    call ptrdict_register_point_property(m, c_loc(this%r0(1)), CSTR("r0"), &
         CSTR("Hook position."))
    call ptrdict_register_point_property(m, c_loc(this%d(1)), CSTR("d"), &
         CSTR("Hook direction (for line and plane only)."))

  endsubroutine harmonic_hook_register

endmodule harmonic_hook
