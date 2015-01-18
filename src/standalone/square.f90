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
!! The Square interpolation kernel
!!
!! The Square interpolation kernel
!<

#include "macros.inc"
#include "filter.inc"

module square
  use libAtoms_module

  type square_t

     real(DP)  :: cutoff  = 1.0_DP

     real(DP)  :: factor

  endtype square_t


  interface init
     module procedure square_init
  endinterface

  interface del
     module procedure square_del
  endinterface

  interface get_cutoff
     module procedure square_get_cutoff
  endinterface

  interface value_and_derivative
     module procedure square_value_and_derivative
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine square_init(this, cutoff)
    implicit none

    type(square_t), intent(inout)     :: this
    real(DP), intent(in), optional  :: cutoff

    ! ---

    if (present(cutoff)) then
       this%cutoff  = cutoff
    endif

    this%factor  = 15.0_DP/(2*PI*this%cutoff**5);

  endsubroutine square_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine square_del(this)
    implicit none

    type(square_t), intent(inout)  :: this

    ! ---

  endsubroutine square_del


  !>
  !! Return the absolute cutoff
  !!
  !! Return the absolute cutoff
  !<
  real(DP) function square_get_cutoff(this)
    implicit none

    type(square_t), intent(in)  :: this

    ! ---

    square_get_cutoff  = this%cutoff

  endfunction square_get_cutoff


  !>
  !! Compute value and derivative
  !!
  !! Compute value and derivative
  !<
  subroutine square_value_and_derivative(this, r, v, w)
    implicit none

    type(square_t), intent(inout)  :: this
    real(DP), intent(in)         :: r
    real(DP), intent(out)        :: v
    real(DP), intent(out)        :: w

    ! ---

    real(DP)  :: h

    ! ---

    w  = (this%cutoff - r)**2;

    v  = this%factor * (this%cutoff - r)**2;
    w  = 2*this%factor * (this%cutoff/r - 1.0_DP);

  endsubroutine square_value_and_derivative

endmodule square

