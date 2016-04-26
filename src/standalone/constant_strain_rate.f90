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
!   classtype:constant_strain_rate_t classname:ConstantStrainRate interface:callables
! @endmeta

!>
!! Constant strain rate deformation
!!
!! Constant strain rate deformation
!<

#include "macros.inc"

module constant_strain_rate
  use supplib

  use particles
  use neighbors
  use dynamics
  
  implicit none

  private

  public :: constant_strain_rate_t
  type constant_strain_rate_t

     !
     ! Hydrostatic pressure components
     !

     real(DP)  :: strain_rate(3) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     real(DP)  :: velocity(3) = (/ 0.0_DP, 0.0_DP, 0.0_DP /)

  endtype constant_strain_rate_t


  public :: init
  interface init
     module procedure constant_strain_rate_init
  endinterface

  public :: invoke
  interface invoke
     module procedure constant_strain_rate_invoke
  endinterface

  public :: register
  interface register
    module procedure constant_strain_rate_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a constant_strain_rate object
  !<
  subroutine constant_strain_rate_init(this, strain_rate, velocity)
    implicit none

    type(constant_strain_rate_t), intent(inout) :: this
    real(DP),           optional, intent(in)    :: strain_rate(3)
    real(DP),           optional, intent(in)    :: velocity(3)

    ! ---

    ASSIGN_PROPERTY(strain_rate)
    ASSIGN_PROPERTY(velocity)

    call prlog("- constant_strain_rate_init -")
    call prlog("strain_rate = " // this%strain_rate(1) // " " // this%strain_rate(2) // " " // this%strain_rate(3))
    call prlog("velocity    = " // this%velocity(1) // " " // this%velocity(2) // " " // this%velocity(3))
    call prlog

  endsubroutine constant_strain_rate_init


  !>
  !! Adjuste the pressure
  !!
  !! Adjuste the pressure
  !<
  subroutine constant_strain_rate_invoke(this, dyn, nl, ierror)
    implicit none

    type(constant_strain_rate_t), intent(inout) :: this
    type(dynamics_t),             intent(inout) :: dyn
    type(neighbors_t),            intent(in)    :: nl
    integer,            optional, intent(inout) :: ierror

    ! ---

    integer  :: i
    real(DP) :: Abox(3, 3), s(3)

    call timer_start('constant_strain_rate_invoke')

    Abox = dyn%p%Abox
    s = 1.0_DP + this%strain_rate*dyn%dt
    forall(i = 1: 3)
       s(i) = s(i) + this%velocity(i)*dyn%dt/Abox(i, i)
    endforall

    do i = 1, 3
#ifndef IMPLICIT_R
       POS(dyn%p, :, i) = POS(dyn%p, :, i) * s(i)
#endif
       PNC(dyn%p, :, i) = PNC(dyn%p, :, i) * s(i)
       Abox(i, i)       = dyn%p%Abox(i, i) * s(i)
    enddo

    call set_cell(dyn%p, Abox, error=ierror)
    PASS_ERROR(ierror)

    call timer_stop('constant_strain_rate_invoke')

  endsubroutine constant_strain_rate_invoke


  !>
  !! Registry
  !!
  !! Registry
  !<
  subroutine constant_strain_rate_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(constant_strain_rate_t), target, intent(inout) :: this
    type(c_ptr),                          intent(in)    :: cfg
    type(c_ptr),                          intent(out)   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("ConstantStrainRate"), &
         CSTR("Impose constant strain rate or cell velocity."))

    call ptrdict_register_point_property(m, c_loc(this%strain_rate(1)), CSTR("strain_rate"), &
         CSTR("Strain rate."))

    call ptrdict_register_point_property(m, c_loc(this%velocity(1)), CSTR("velocity"), &
         CSTR("Velocity."))

  endsubroutine constant_strain_rate_register

endmodule constant_strain_rate
