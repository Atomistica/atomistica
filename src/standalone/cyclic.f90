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
!**********************************************************************
! Box size information
!
! FIXME!!! The whole code needs to be updated for nonorthogonal boxes
!**********************************************************************
module cyclic
  use supplib

  use particles

  implicit none

  private

  public :: cyclic_get_repetition, cyclic_in_reciprocal_bounds, cyc, cyc2

contains

  !**********************************************************************
  ! Get the number of repetitions of the unit cell such that
  ! a maximum distance of cutoff is covered.
  !**********************************************************************
  subroutine cyclic_get_repetition(p, cutoff, x_max, y_max, z_max)
    implicit none

    type(particles_t), intent(in)  :: p
    real(8), intent(in)   :: cutoff
    integer, intent(out)  :: x_max, y_max, z_max

    ! ---

    x_max = int(cutoff / sqrt(dot_product(p%Abox(1, :), p%Abox(1, :))))+1
    y_max = int(cutoff / sqrt(dot_product(p%Abox(2, :), p%Abox(2, :))))+1
    z_max = int(cutoff / sqrt(dot_product(p%Abox(3, :), p%Abox(3, :))))+1

  endsubroutine cyclic_get_repetition


  !**********************************************************************
  ! Project r into the reciprocal box
  !**********************************************************************
  function cyclic_in_reciprocal_bounds(p, r) result(cyc)
    implicit none

    type(particles_t), intent(in)  :: p
    real(DP), intent(in)  :: r(3)
    
    real(DP)              :: cyc(3)

    ! ---

    real(DP)  :: s(3)

    s    = matmul(transpose(p%Abox), r)
    s    = s - nint(s)
    cyc  = matmul(transpose(p%Bbox), s)

  endfunction cyclic_in_reciprocal_bounds


  !**********************************************************************
  ! Project x, y, z into the box
  !**********************************************************************
  subroutine cyc(p,x,y,z)
    implicit none

    type(particles_t), intent(in)  :: p
    real(DP), intent(inout)  :: x, y, z

    ! ---

    real(DP)  :: s(3)

    s = p%Bbox(:,1)*x+p%Bbox(:,2)*y+p%Bbox(:,3)*z
    s = s-nint(s)

    x = p%Abox(1,1)*s(1)+p%Abox(1,2)*s(2)+p%Abox(1,3)*s(3)
    y = p%Abox(2,1)*s(1)+p%Abox(2,2)*s(2)+p%Abox(2,3)*s(3)	
    z = p%Abox(3,1)*s(1)+p%Abox(3,2)*s(2)+p%Abox(3,3)*s(3)

  endsubroutine cyc


  !**********************************************************************
  ! Project x, y, z into the box
  !**********************************************************************
  subroutine cyc2(p, x,y,z)
    implicit none

    type(particles_t), intent(in)  :: p
    real(DP), intent(inout)  :: x, y, z

    ! ---

    real(DP)  :: s(3)

    s = p%Bbox(:,1)*x+p%Bbox(:,2)*y+p%Bbox(:,3)*z
    s = s-int(s)

    x = p%Abox(1,1)*s(1)+p%Abox(1,2)*s(2)+p%Abox(1,3)*s(3)
    y = p%Abox(2,1)*s(1)+p%Abox(2,2)*s(2)+p%Abox(2,3)*s(3)	
    z = p%Abox(3,1)*s(1)+p%Abox(3,2)*s(2)+p%Abox(3,3)*s(3)

  endsubroutine cyc2

endmodule cyclic
          

