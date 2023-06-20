!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
!! and others. See the AUTHORS file in the top-level Atomistica directory.
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
!   classtype:born_mayer_t classname:BornMayer interface:potentials
! @endmeta

!>
!! The Born-Mayer potential
!!
!! The Born-Mayer potential
!<

#include "macros.inc"
#include "filter.inc"

module born_mayer
  use libAtoms_module

  use ptrdict

  use logging
  use timer

  use particles
  use neighbors
  use filter

  implicit none

  private

  public :: born_mayer_t
  type born_mayer_t

     !
     ! Potential parameters
     !

     real(DP)  :: A = 1
     real(DP)  :: rho = 1
     real(DP)  :: cutoff = 1

     !
     ! Element on which to apply the force
     !

     character(MAX_EL_STR)  :: element1 = "C"
     character(MAX_EL_STR)  :: element2 = "C"
     integer                :: el1
     integer                :: el2

     !
     ! Shift the potential at the cut-off
     !

     real(DP) :: shift
     
  endtype born_mayer_t


  public :: init
  interface init
     module procedure born_mayer_init
  endinterface

  public :: del
  interface del
     module procedure born_mayer_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure born_mayer_bind_to
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure born_mayer_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure born_mayer_register
  endinterface

contains


  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine born_mayer_init(this, element1, element2, A, rho, cutoff)
    implicit none

    type(born_mayer_t), intent(inout)     :: this
    character(*), intent(in), optional    :: element1
    character(*), intent(in), optional    :: element2
    real(DP), intent(in), optional        :: A
    real(DP), intent(in), optional        :: rho
    real(DP), intent(in), optional        :: cutoff

    ! ---

    ASSIGN_PROPERTY(element1)
    ASSIGN_PROPERTY(element2)
    ASSIGN_PROPERTY(A)
    ASSIGN_PROPERTY(rho)
    ASSIGN_PROPERTY(cutoff)

  endsubroutine born_mayer_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine born_mayer_del(this)
    implicit none

    type(born_mayer_t), intent(inout)  :: this

    ! ---

  endsubroutine born_mayer_del


  !>
  !! Initialization
  !!
  !! Constructs the parameter sets
  !<
  subroutine born_mayer_bind_to(this, p, nl, ierror)
    implicit none

    type(born_mayer_t), intent(inout)  :: this
    type(particles_t), intent(in)      :: p
    type(neighbors_t), intent(inout)   :: nl
    integer, intent(inout), optional   :: ierror

    ! ---

    write (ilog, '(A)')  "- born_mayer_init -"

    write (ilog, '(5X,A,A,A)') trim(this%element1), "-", trim(this%element2)

    write (ilog, '(5X,A,F7.3)')  "A = ", this%A
    write (ilog, '(5X,A,F7.3)')  "rho = ", this%rho
    write (ilog, '(5X,A,F7.3)')  "cutoff = ", this%cutoff

    this%el1 = filter_from_string(this%element1, p)
    this%el2 = filter_from_string(this%element2, p)

    call request_interaction_range(nl, this%cutoff)
    this%shift = this%A * exp(-this%cutoff/this%rho)

    write (ilog, *)

  endsubroutine born_mayer_bind_to


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine born_mayer_energy_and_forces(this, p, nl, epot, for, wpot, epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)
    implicit none

    type(born_mayer_t), intent(inout)  :: this
    type(particles_t), intent(in)      :: p
    type(neighbors_t), intent(inout)   :: nl
    real(DP), intent(inout)            :: epot
    real(DP), intent(inout)            :: for(3, p%maxnatloc)
    real(DP), intent(inout)            :: wpot(3, 3)
    real(DP), intent(inout), optional  :: epot_per_at(p%maxnatloc)
    real(DP), intent(inout), optional  :: epot_per_bond(nl%neighbors_size)
    real(DP), intent(inout), optional  :: f_per_bond(3, nl%neighbors_size)
#ifdef LAMMPS
    real(DP), intent(inout), optional  :: wpot_per_at(6, p%maxnatloc)
    real(DP), intent(inout), optional  :: wpot_per_bond(6, nl%neighbors_size)
#else
    real(DP), intent(inout), optional  :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), intent(inout), optional  :: wpot_per_bond(3, 3, nl%neighbors_size)
#endif
    integer, intent(inout), optional   :: ierror

    ! ---

    integer   :: i, ni, j
    real(DP)  :: dr(3), abs_dr, exp_r, f(3), e

    ! ---

    call timer_start("born_mayer_energy_and_forces")

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    e = 0.0_DP
    do i = 1, p%natloc

       if (IS_EL(this%el1, p, i)) then

          do ni = nl%seed(i), nl%last(i)
             DISTJ_SQ(p, nl, i, ni, j, dr, abs_dr)

             if (i <= j .and. IS_EL(this%el2, p, j)) then
                if (abs_dr < this%cutoff**2) then
                   abs_dr  = sqrt(abs_dr)

                   exp_r = exp(-abs_dr/this%rho)

                   e = e + this%A*exp_r - this%shift
                   f = (this%A / this%rho)*exp_r*dr/abs_dr

                   VEC3(for, i) = VEC3(for, i) + f
                   VEC3(for, j) = VEC3(for, j) - f
                endif
             endif
          enddo

       else if (IS_EL(this%el2, p, i)) then

          do ni = nl%seed(i), nl%last(i)
             DISTJ_SQ(p, nl, i, ni, j, dr, abs_dr)

             if (i <= j .and. IS_EL(this%el1, p, j)) then
                if (abs_dr < this%cutoff**2) then
                   abs_dr  = sqrt(abs_dr)

                   exp_r = exp(-abs_dr/this%rho)

                   e = e + this%A*exp_r - this%shift
                   f = (this%A / this%rho)*exp_r*dr/abs_dr

                   VEC3(for, i) = VEC3(for, i) + f
                   VEC3(for, j) = VEC3(for, j) - f
                endif
             endif
          enddo

       endif

    enddo
    epot = epot + e

    call timer_stop("born_mayer_energy_and_forces")

  endsubroutine born_mayer_energy_and_forces

  subroutine born_mayer_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(born_mayer_t), target, intent(inout)     :: this
    type(c_ptr), intent(in)                  :: cfg
    type(c_ptr), intent(out)                 :: m

    ! ---

    m = ptrdict_register_section(cfg, "BornMayer" // char(0), &
         "The Born-Mayer potential." // char(0))

    call ptrdict_register_real_property(m, c_loc(this%A), "A" // char(0), &
         "Interaction energy." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%rho), "rho" // char(0), &
         "Interaction diameter." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%cutoff), "cutoff" // char(0), &
         "Potential cutoff: If smaller than zero, the cutoff is set such that the potential is only repulsive." // char(0))

    call ptrdict_register_string_property(m, c_loc(this%element1), MAX_EL_STR, "element1" // char(0), &
         "First element." // char(0))

    call ptrdict_register_string_property(m, c_loc(this%element2), MAX_EL_STR, "element2" // char(0), &
         "Second element." // char(0))

  endsubroutine born_mayer_register

endmodule born_mayer
