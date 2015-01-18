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
!   classtype:constant_force_t classname:ConstantForce interface:potentials
! @endmeta

!>
!! Constant force field
!!
!! Constant force field (e.g. gravitation)
!<

#include "macros.inc"

module constant_force
  use libAtoms_module

  use logging

  use particles
  use dynamics
  use neighbors

  implicit none

  private

  public :: constant_force_t
  type constant_force_t

     !
     ! Particle group
     !

     integer   :: g  = 1
     
     !
     ! Magnitude
     !

     real(DP)  :: f(3)     = [ 0.0_DP, 0.0_DP, 1.0_DP ]

     !
     ! Spatial variation
     !

     integer   :: n(3)     = [ 0, 0, 0 ]
     real(DP)  :: alpha(3) = [ 0.0_DP, 0.0_DP, 0.0_DP ]

     !
     ! Oszillation freqency
     !

     real(DP)  :: freq     = 0.0_DP

  endtype constant_force_t


  public :: init
  interface init
     module procedure constant_force_init
  endinterface

  public :: energy_and_forces_with_dyn
  interface energy_and_forces_with_dyn
     module procedure constant_force_energy_and_forces
  endinterface

  public :: register
  interface register
    module procedure constant_force_register
  endinterface

contains

  !>
  !! Initialize a ConstantForce object
  !<
  subroutine constant_force_init(this)
    implicit none

    type(constant_force_t), intent(inout)  :: this

    ! ---

    call prlog("- constant_force_init -")
    call prlog("     $Id$")

    call prlog("     group = " // this%g)
    call prlog("     f     = " // this%f)
    call prlog("     n     = " // this%n)
    call prlog("     alpha = " // this%alpha)
    call prlog("     freq  = " // this%freq)

    call prlog

  endsubroutine constant_force_init


  !>
  !! Compute the force
  !<
  subroutine constant_force_energy_and_forces(this, dyn, nl, ierror)
    implicit none

    type(constant_force_t), intent(in)    :: this
    type(dynamics_t),       intent(inout) :: dyn
    type(neighbors_t),      intent(in)    :: nl
    integer,      optional, intent(inout) :: ierror

    ! ---

    integer   :: i
    real(DP)  :: q(3), freq, c(3), a(3)

    ! ---

    q = 2*PI*this%n
    freq = 2*PI*this%freq
    c = ( dyn%p%Abox(1, :) + dyn%p%Abox(2, :) + dyn%p%Abox(3, :) )/3
    a = this%alpha*this%alpha

    do i = 1, dyn%p%nat
       if (dyn%p%g(i) == this%g) then
          VEC3(dyn%f, i) = VEC3(dyn%f, i) + this%f * &
               cos(dot_product(q, matmul(dyn%p%Bbox, PNC3(dyn%p, i))) - freq*dyn%ti) * &
               exp(-a * (POS3(dyn%p, i)-c)**2)
       endif
    enddo

  endsubroutine constant_force_energy_and_forces


  subroutine constant_force_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(constant_force_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                   :: cfg
    type(c_ptr), intent(out)                  :: m

    ! ---

    this%g     = 1
    this%f(:)  = 0.0

    m = ptrdict_register_section(cfg, CSTR("ConstantForce"), &
         CSTR("Constant force field (e.g. gravitation)."))

    call ptrdict_register_integer_property(m, c_loc(this%g), CSTR("group"), &
         CSTR("Group of particles on which this force should act."))

    call ptrdict_register_point_property(m, c_loc(this%f(1)), CSTR("f"), &
         CSTR("The force magnitude and direction"))
    call ptrdict_register_intpoint_property(m, c_loc(this%n(1)), CSTR("n"), &
         CSTR("Spatial variation: Fourier component"))
    call ptrdict_register_point_property(m, c_loc(this%alpha(1)), CSTR("alpha"), &
         CSTR("Spatial variation: Decay length (from middle of simulation cell)."))
    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Oscillation frequency of the force."))

  endsubroutine constant_force_register

endmodule constant_force
