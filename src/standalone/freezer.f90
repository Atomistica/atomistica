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
!   private broken
! @endmeta

!>
!! Freeze particles according to some criterium
!<

#include "macros.inc"

module freezer
  use supplib

  use particles
  use neighbors
  use dynamics

  implicit none

  type freezer_t
     
     integer   :: enabled  = 0

     real(DP)  :: z        = 10.0  ! Freeze particles above this coordinate

  endtype freezer_t

contains

  !**********************************************************************
  ! Freeze particles
  !**********************************************************************
  subroutine freezer_invoke(this, dyn, nl)
    implicit none

    type(freezer_t), intent(inout)   :: this
    type(dynamics_t), intent(inout)  :: dyn
    type(neighbors_t), intent(in)    :: nl

    ! ---

    integer   :: i

    ! ---

    do i = 1, dyn%p%nat

       if ( POS(dyn%p, i, 3) > this%z ) then
          dyn%p%g(i)           = 0

          VEC3(dyn%v, i)  = 0.0_DP
          VEC3(dyn%f, i)  = 0.0_DP
       endif

    enddo

  endsubroutine freezer_invoke


  !****************************************************************
  ! Initialize the property list
  !****************************************************************
  subroutine freezer_register(this, cfg, m)
    implicit none

    type(freezer_t), intent(inout)  :: this
    CPOINTER, intent(in)            :: cfg
    CPOINTER, intent(out)           :: m

    ! ---

    m = ptrdict_register_module(cfg, this%enabled, "Freezer" // char(0), &
         "Freeze particles with large z-coordinates." // char(0))

    call ptrdict_register_real_property(m, this%z, "z" // char(0), &
         "Freeze particles above this z coordinate." // char(0))

  endsubroutine freezer_register
  
endmodule freezer
