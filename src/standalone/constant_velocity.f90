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
!   classtype:constant_velocity_t classname:ConstantVelocity
!   interface:callables
! @endmeta

!>
!! Constant velocity for a group of particles
!!
!! Constant velocity for a group of particles
!<

#include "macros.inc"

module constant_velocity
  use libAtoms_module

  use io
  use logging

  use particles
  use neighbors
  use dynamics

  implicit none

  private

  public :: constant_velocity_t
  type constant_velocity_t

     !
     ! Output forces?
     !

     real(DP)  :: out_freq  = -1.0_DP
     integer   :: un
     
     !
     ! Particle group
     !

     integer   :: g  = 1
     
     !
     ! Magnitude
     !

     real(DP)  :: v(3)  = 0.0_DP

     !
     ! Force for force output
     !

     real(DP)  :: f(3)

     !
     ! Time
     !

     real(DP)  :: ti

  endtype constant_velocity_t


  public :: init
  interface init
     module procedure constant_velocity_init
  endinterface

  public :: del
  interface del
     module procedure constant_velocity_del
  endinterface

  public :: invoke
  interface invoke
     module procedure constant_velocity_invoke
  endinterface

  public :: register
  interface register
    module procedure constant_velocity_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a constant_velocity object
  !<
  subroutine constant_velocity_init(this)
    implicit none

    type(constant_velocity_t), intent(inout)  :: this

    ! ---

    call prlog("- constant_velocity")
    call prlog("     v =  " // this%v)
    call prlog("     g =  " // this%g)

    if (this%out_freq > 0) then
       this%un  = fopen("constant_velocity_force.out", F_WRITE)
       this%f   = 0.0_DP
    endif

    this%ti  = 0.0_DP

    call prlog

  endsubroutine constant_velocity_init


  !>
  !! Destructor
  !!
  !! Delete a constant_velocity object
  !<
  subroutine constant_velocity_del(this)
    implicit none

    type(constant_velocity_t), intent(inout)  :: this

    ! ---

    if (this%out_freq > 0) then
       call fclose(this%un)
    endif

  endsubroutine constant_velocity_del


  !>
  !! Fix the velocity of a group of particles
  !!
  !! Fix the velocity of a group of particles
  !<
  subroutine constant_velocity_invoke(this, dyn, nl, ierror)
    implicit none

    type(constant_velocity_t), intent(inout)   :: this
    type(dynamics_t), intent(inout)            :: dyn
    type(neighbors_t), intent(in)              :: nl
    integer, intent(inout), optional           :: ierror

    ! ---

    integer   :: i
    real(DP)  :: f(3)

    ! ---

    this%ti  = this%ti + dyn%dt

    f  = 0.0_DP
    do i = 1, dyn%p%nat
       if (dyn%p%g(i) == this%g) then
          f  = f + VEC3(dyn%f, i)

          VEC3(dyn%v, i)  = this%v
          VEC3(dyn%f, i)  = 0.0_DP
       endif
    enddo

    this%f  = this%f + f*dyn%dt

    if (this%out_freq > 0) then
       if (this%ti >= this%out_freq) then
          write (this%un, '(I10,4ES20.10)')  dyn%it, dyn%ti, this%f/this%ti
          this%f   = 0.0_DP
          this%ti  = 0.0_DP
       endif
    endif

  endsubroutine constant_velocity_invoke


  subroutine constant_velocity_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(constant_velocity_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                   :: cfg
    type(c_ptr), intent(out)                  :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("ConstantVelocity"), &
         CSTR("Constant velocity field."))

    call ptrdict_register_integer_property(m, c_loc(this%g), CSTR("group"), &
         CSTR("Group of particles for which the velocity should be fixed."))

    call ptrdict_register_point_property(m, c_loc(this%v(1)), CSTR("v"), &
         CSTR("The velocity."))

    call ptrdict_register_real_property(m, c_loc(this%out_freq), CSTR("out_freq"), &
         CSTR("Interval in which to output the force necessary to drive that motion."))

  endsubroutine constant_velocity_register

endmodule constant_velocity
