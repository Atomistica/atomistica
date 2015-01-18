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
!   classtype:constraints_t classname:Constraints interface:callables
! @endmeta

!>
!! Constrain coordinates according to certain rules. Right now,
!! we can remove:
!!  - Linear velocities and forces
!!  - Angular momentum and torque
!<

#include "macros.inc"

module constraints
  use libAtoms_module

  use particles
  use neighbors
  use dynamics

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  character(*), parameter, private  :: MODULE_STR = "Constraints"

  public :: constraints_t
  type constraints_t
     
     integer        :: interval  = 1   !< Interval in which to remove the momenta
     integer        :: group     = -1  !< Group of atoms

     logical(BOOL)  :: linear    = .true.   !< Remove linear velocities and forces
     logical(BOOL)  :: angular   = .false.  !< Remove angular velocities and forces

  endtype constraints_t


  public :: invoke
  interface invoke
     module procedure constraints_invoke
  endinterface

  public :: register
  interface register
    module procedure constraints_register
  endinterface

contains

  !>
  !! Apply constraints
  !!
  !! Apply constraints
  !<
  subroutine constraints_invoke(this, dyn, nl, ierror)
    implicit none

    type(constraints_t), intent(inout)  :: this
    type(dynamics_t), intent(inout)     :: dyn
    type(neighbors_t), intent(in)       :: nl
    integer, intent(inout), optional    :: ierror

    ! ---

    integer   :: i, j, k, n, info, ipiv(3)
    real(DP)  :: mass, momentum(3), force(3), Iij(3, 3), Ltot(3), Mtot(3)

    ! ---

    call timer_start("constraints_invoke")

    if ( (this%interval == -1 .and. dyn%it == 1) .or. (this%interval > 0 .and. mod(dyn%it, this%interval) == 0) ) then

       if (this%linear) then

          n         = 0
          mass      = 0.0_DP
          momentum  = 0.0_DP
          force     = 0.0_DP
          do i = 1, dyn%p%natloc
             if (this%group < 0 .or. dyn%p%g(i) == this%group) then
                n         = n + 1
                mass      = mass + dyn%p%m(i)
                momentum  = momentum + dyn%p%m(i)*VEC3(dyn%v, i)
                force     = force + VEC3(dyn%f, i)
             endif
          enddo

#ifdef _MP
          call sum_in_place(mod_communicator%mpi, mass)
          call sum_in_place(mod_communicator%mpi, momentum)
          call sum_in_place(mod_communicator%mpi, force)
          call sum_in_place(mod_communicator%mpi, n)
#endif

          if (n > 0) then

             momentum  = momentum/mass
             force     = force/n

             do i = 1, dyn%p%natloc
                if (this%group < 0 .or. dyn%p%g(i) == this%group) then
                   VEC3(dyn%v, i)  = VEC3(dyn%v, i) - momentum
                   VEC3(dyn%f, i)  = VEC3(dyn%f, i) - force
                endif
             enddo

          else

             WARN("No matching atoms found.")

          endif
       endif

       if (this%angular) then
       
          !
          ! Determine angular momentum (Ltot) and torque (Mtot)
          !
          
          Ltot  = 0.0_DP
          Mtot  = 0.0_DP
          do i = 1, dyn%p%natloc
             if (this%group < 0 .or. dyn%p%g(i) == this%group) then
                Ltot  = Ltot + &
                     dyn%p%m(i)*cross_product(POS3(dyn%p, i), VEC3(dyn%v, i))
                Mtot  = Mtot + &
                     cross_product(POS3(dyn%p, i), VEC3(dyn%f, i))
             endif
          enddo

          !
          ! Calculate inertia tensor
          !

          Iij = 0
          do i = 1, 3
             do j = 1, 3
                do k = 1, dyn%p%natloc
                   if (this%group < 0 .or. dyn%p%g(k) == this%group) then
                      Iij(i, j) = Iij(i, j) - POS(dyn%p, k, i)*POS(dyn%p, k, j)*dyn%p%m(k)
                      if (i == j)  Iij(i, j) = Iij(i, j) + dyn%p%m(k)*dot_product(POS3(dyn%p, k), POS3(dyn%p, k))
                   endif
                enddo
             enddo
          enddo

#ifdef _MP
          call sum_in_place(mod_communicator%mpi, Ltot)
          call sum_in_place(mod_communicator%mpi, Mtot)
          call sum_in_place(mod_communicator%mpi, Iij)
#endif

          !
          ! Solve the equation I*omega = Ltot and I*alpha = Mtot
          !

          call dgesv(3, 1, Iij, 3, ipiv, Ltot, 3, info)
          if (info /= 0) then
             RAISE_ERROR("dgesv failed.", ierror)
          endif
          
          call dgesv(3, 1, Iij, 3, ipiv, Mtot, 3, info)
          if (info /= 0) then
             RAISE_ERROR("dgesv failed.", ierror)
          endif

          do i = 1, dyn%p%natloc
             if (this%group < 0 .or. dyn%p%g(i) == this%group) then
                VEC3(dyn%v, i)  = VEC3(dyn%v, i) - &
                     cross_product(Ltot, POS3(dyn%p, i))
                VEC3(dyn%f, i)  = VEC3(dyn%f, i) - &
                     dyn%p%m(i) * cross_product(Mtot, POS3(dyn%p, i))
             endif
          enddo

       endif

    endif

    call timer_stop("constraints_invoke")

  endsubroutine constraints_invoke


  subroutine constraints_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(constraints_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                    :: cfg
    type(c_ptr), intent(out)                   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Constraints"), &
         CSTR("Constraints module."))

    call ptrdict_register_integer_property(m, c_loc(this%interval), CSTR("interval"), &
         CSTR("Interval in which to remove linear and/or angular momentum (-1 = only at the beginning)."))

    call ptrdict_register_integer_property(m, c_loc(this%group), CSTR("group"), &
         CSTR("Group of atoms for which to apply constrains (all atoms if < 0)."))

    call ptrdict_register_boolean_property(m, c_loc(this%linear), CSTR("linear"), &
         CSTR("Remove linear velocities and forces."))

    call ptrdict_register_boolean_property(m, c_loc(this%angular), CSTR("angular"), &
         CSTR("Remove angular momentums and torques."))

  endsubroutine constraints_register

endmodule constraints
