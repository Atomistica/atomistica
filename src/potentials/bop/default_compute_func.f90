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

  !>
  !! Compute energy, force, virial etc.
  !<
  subroutine COMPUTE_FUNC(this, p, nl, epot, f, wpot, mask, epot_per_at, &
       epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)
    implicit none

    type(BOP_TYPE),     intent(inout) :: this
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(inout) :: epot
    real(DP),           intent(inout) :: f(3, p%maxnatloc)  !< forces
    real(DP),           intent(inout) :: wpot(3, 3)
    integer,  optional, intent(in)    :: mask(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_at(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_bond(nl%neighbors_size)
    real(DP), optional, intent(inout) :: f_per_bond(3, nl%neighbors_size)
#ifdef LAMMPS
    real(DP), optional, intent(inout) :: wpot_per_at(6, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot_per_bond(6, nl%neighbors_size)
#else
    real(DP), optional, intent(inout) :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot_per_bond(3, 3, nl%neighbors_size)
#endif
    integer,  optional, intent(out)   :: ierror

    ! ---

    integer  :: i, d, el(p%maxnatloc), nebmax, nebavg

    ! ---

    INIT_ERROR(ierror)

    call timer_start(BOP_NAME_STR // "_force")

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    ! Internal element numbers
    el = -1
    nebmax = 0
    nebavg = 0
    do i = 1, p%nat
       if (p%el2Z(p%el(i)) > 0) then
         el(i) = this%Z2db(p%el2Z(p%el(i)))
       endif
       d = nl%last(i)-nl%seed(i)+1
       nebmax = max(nebmax, d)
       nebavg = nebavg + d
    enddo
    nebavg = (nebavg+1)/max(p%nat, 1)+1

#ifdef LAMMPS
    call BOP_KERNEL( &
         this, &
         p%maxnatloc, p%natloc, p%nat, p%r_non_cyc, el, &
         nebmax, nebavg, nl%seed, nl%last, nl%neighbors, nl%neighbors_size, &
         epot, f, wpot, mask, &
         epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
         ierror)
#else
    call BOP_KERNEL( &
         this, p%Abox, &
         p%maxnatloc, p%natloc, p%nat, p%r_non_cyc, el, &
         nebmax, nebavg, nl%seed, nl%last, nl%neighbors, nl%neighbors_size, &
         nl%dc, p%shear_dx, &
         epot, f, wpot, mask, &
         epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
         ierror)
#endif
    PASS_ERROR(ierror)

    call timer_stop(BOP_NAME_STR // "_force")

  endsubroutine COMPUTE_FUNC
