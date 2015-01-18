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
!   classtype:remove_rotation_t classname:RemoveRotation interface:callables
! @endmeta

!>
!! Remove the rotation around an axis
!!
!! Remove the rotation around an axis
!<

#include "macros.inc"
#include "filter.inc"

module remove_rotation
  use libAtoms_module

  use particles
  use filter
  use neighbors
  use dynamics

#ifdef _MP
  use mpi
  use communicator
#endif

  implicit none

  private

  public :: remove_rotation_t
  type remove_rotation_t

     type(particles_t), pointer  :: p  => NULL()
     
     character(MAX_EL_STR)       :: elements  = "*"
     integer                     :: els
     
     real(DP)                    :: r0(3)  = (/ 0.0_DP, 0.0_DP, 0.0_DP /)   !< Center of rotation
     real(DP)                    :: d(3)   = (/ 0.0_DP, 0.0_DP, 1.0_DP /)   !< Rotation axis

  endtype remove_rotation_t


  public :: init
  interface init
     module procedure remove_rotation_init
  endinterface

  public :: invoke
  interface invoke
     module procedure remove_rotation_invoke
  endinterface

  public :: register
  interface register
    module procedure remove_rotation_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine remove_rotation_init(this)
    implicit none

    type(remove_rotation_t), intent(inout)  :: this

    ! ---

    call prlog("- remove_rotation_init -")

    this%p  => NULL()

    this%d    = this%d / sqrt(dot_product(this%d, this%d))

    call prlog("     " // this%elements)
    call prlog("     d  = " // this%d)

    call prlog

  endsubroutine remove_rotation_init


  !>
  !! Apply remove_rotation
  !!
  !! Apply remove_rotation
  !<
  subroutine remove_rotation_invoke(this, dyn, nl, ierror)
    implicit none

    type(remove_rotation_t), intent(inout)  :: this
    type(dynamics_t), target                :: dyn
    type(neighbors_t), intent(in)           :: nl
    integer, intent(inout), optional        :: ierror

    ! ---

    integer   :: i, k, n, info, ipiv(2)
    real(DP)  :: help, omega, alpha
    real(DP)  :: Iij(3, 3), Ltot(3), Mtot(3), r(3)

    ! ---

    if (.not. associated(this%p, dyn%p)) then
       this%p     => dyn%p
       this%els   = filter_from_string(this%elements, dyn%p)
       dyn%p%dof  = dyn%p%dof - 2
    endif


    !
    ! Determine angular momentum (Ltot) and torque (Mtot)
    !

    Ltot  = 0.0_DP
    Mtot  = 0.0_DP
    do i = 1, dyn%p%natloc
       if (dyn%p%g(i) > 0 .and. IS_EL(this%els, dyn%p, i)) then
          r     = POS3(dyn%p, i) - this%r0

          Ltot  = Ltot + dyn%p%m(i)*cross_product(r, VEC3(dyn%v, i))
          Mtot  = Mtot + cross_product(r, VEC3(dyn%f, i))
       endif
    enddo


    !
    ! Calculate inertia tensor
    !

    Iij = 0.0_DP
    do k = 1, dyn%p%natloc
       if (dyn%p%g(k) > 0 .and. IS_EL(this%els, dyn%p, k)) then
          r     = POS3(dyn%p, k) - this%r0

          Iij   = Iij - dyn%p%m(k)*outer_product(r, r)

          help  = dyn%p%m(k)*dot_product(r, r)
          do i = 1, 3
             Iij(i, i) = Iij(i, i) + help
          enddo
       endif
    enddo

    !
    ! Solve the equation I*omega = Ltot and I*alpha = Mtot
    !

    omega  = dot_product(Ltot, this%d) / (dot_product(this%d, matmul(Iij, this%d)))
    alpha  = dot_product(Mtot, this%d) / (dot_product(this%d, matmul(Iij, this%d)))

    do i = 1, dyn%p%natloc
       if (dyn%p%g(i) > 0 .and. IS_EL(this%els, dyn%p, i)) then
          r  = POS3(dyn%p, i) - this%r0

          VEC3(dyn%v, i)  = VEC3(dyn%v, i) - omega*cross_product(this%d, r)
          VEC3(dyn%f, i)  = VEC3(dyn%f, i) - &
               dyn%p%m(i) * alpha*cross_product(this%d, r)
       endif
    enddo

  endsubroutine remove_rotation_invoke


  subroutine remove_rotation_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(remove_rotation_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                 :: cfg
    type(c_ptr), intent(out)                :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("RemoveRotation"), &
         CSTR("Remove a rotation along a certain axis."))

    call ptrdict_register_string_property(m, c_loc(this%elements), MAX_EL_STR, CSTR("elements"), &
         CSTR("Elements for which to enable this module."))

    call ptrdict_register_point_property(m, c_loc(this%r0(1)), CSTR("r0"), &
         CSTR("Anchor."))
    call ptrdict_register_point_property(m, c_loc(this%d(1)), CSTR("d"), &
         CSTR("Rotation axis."))

  endsubroutine remove_rotation_register

endmodule remove_rotation
