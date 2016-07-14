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

#include "macros.inc"
#include "filter.inc"

module dense_repulsion
  use, intrinsic :: iso_c_binding

  use supplib

  use nonuniform_spline

  use particles
  use neighbors

  use materials

  use dense_hamiltonian_type
  use dense_hamiltonian

  implicit none

  private

  public  :: E_rep

contains

  !**********************************************************************
  ! Returns the energy related to the short-range repulsive potential
  !**********************************************************************
  function E_rep(tb, db, p, nl) result(res)
    implicit none

    type(dense_hamiltonian_t), intent(in)  :: tb
    type(materials_t), intent(in)          :: db
    type(particles_t), intent(in)          :: p
    type(neighbors_t), intent(in)          :: nl

    real(DP)                               :: res

    ! ---

    integer   :: i, j, ni
    real(DP)  :: dr, erep, erep0
    integer   :: a

    type(notb_element_t), pointer  :: tb_at(:)

    ! --

    call c_f_pointer(tb%at, tb_at, [tb%nat])

    erep = 0
    i_loop: do i = 1, p%natloc

       if (IS_EL(tb%f, p, i)) then

          ni_loop: do ni = nl%seed(i), nl%last(i)
             j = GET_NEIGHBOR(nl, ni)

             if (IS_EL(tb%f, p, j)) then

                if (i <= j) then

                   dr  = GET_ABS_DRJ(p, nl, i, j, ni)
                   if (dr < db%R(tb_at(i)%enr, tb_at(j)%enr)%cut) then
                      a = interval(db%R(tb_at(i)%enr, tb_at(j)%enr), dr)
                      erep0 = f(db%R(tb_at(i)%enr, tb_at(j)%enr), 1, dr, a)
                      if (i == j .or. j > p%natloc) then
                         ! Only count half of the energy if one atom is ghost
                         erep = erep + 0.5_DP*erep0
                      else
                         erep = erep + erep0
                      endif
                   endif

                endif

             endif

          enddo ni_loop

       endif

    enddo i_loop
    
    res = erep

  endfunction E_rep

endmodule dense_repulsion
