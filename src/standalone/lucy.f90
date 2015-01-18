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
!! The Lucy interpolation kernel
!!
!! The Lucy interpolation kernel
!<

#include "macros.inc"
#include "filter.inc"

module lucy
  use libAtoms_module

  type lucy_t

     real(DP)  :: cutoff  = 1.0_DP

     real(DP)  :: factor_i
     real(DP)  :: factor_w

  endtype lucy_t


  interface init
     module procedure lucy_init
  endinterface

  interface del
     module procedure lucy_del
  endinterface

  interface get_cutoff
     module procedure lucy_get_cutoff
  endinterface

  interface value_and_derivative
     module procedure lucy_value_and_derivative
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine lucy_init(this, cutoff)
    implicit none

    type(lucy_t), intent(inout)     :: this
    real(DP), intent(in), optional  :: cutoff

    ! ---

    if (present(cutoff)) then
       this%cutoff  = cutoff
    endif

    this%factor_i  = 105.0_DP/(16*PI*this%cutoff**7);
    this%factor_w  = 315.0_DP/(4*PI*this%cutoff**7);

  endsubroutine lucy_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine lucy_del(this)
    implicit none

    type(lucy_t), intent(inout)  :: this

    ! ---

  endsubroutine lucy_del


  !>
  !! Return the absolute cutoff
  !!
  !! Return the absolute cutoff
  !<
  real(DP) function lucy_get_cutoff(this)
    implicit none

    type(lucy_t), intent(in)  :: this

    ! ---

    lucy_get_cutoff  = this%cutoff

  endfunction lucy_get_cutoff


  !>
  !! Compute value and derivative
  !!
  !! Compute value and derivative
  !<
  subroutine lucy_value_and_derivative(this, r, v, w)
    implicit none

    type(lucy_t), intent(inout)  :: this
    real(DP), intent(in)         :: r
    real(DP), intent(out)        :: v
    real(DP), intent(out)        :: w

    ! ---

    real(DP)  :: h

    ! ---

    w  = (this%cutoff - r)**2;

    v  = this%factor_i * (this%cutoff + 3*r) * (this%cutoff - r) * w;
    w  = this%factor_w * w;

  endsubroutine lucy_value_and_derivative

endmodule lucy

