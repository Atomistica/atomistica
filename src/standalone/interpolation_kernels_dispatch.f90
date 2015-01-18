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
!! Interpolation kernel dispatch module
!!
!! Interpolation kernel dispatch module
!<

#include "macros.inc"

module interpolation_kernels
  use lucy
  use square

  implicit none

  !
  ! Interpolation_kernels
  !

  type interpolation_kernels_t

     type(lucy_t), pointer    :: lucy    => NULL()
     type(square_t), pointer  :: square  => NULL()

  endtype interpolation_kernels_t


  interface del
     module procedure interpolation_kernels_del
  endinterface

  interface get_cutoff
     module procedure interpolation_kernels_get_cutoff
  endinterface

  interface value_and_derivative
     module procedure interpolation_kernels_value_and_derivative
  endinterface

contains

  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine interpolation_kernels_del(this)
    implicit none

    type(interpolation_kernels_t), intent(inout)  :: this

    ! ---

    if (associated(this%lucy)) &
         call del(this%lucy)
    if (associated(this%square)) &
         call del(this%square)

  endsubroutine interpolation_kernels_del


  !>
  !! Return the absolute cutoff
  !!
  !! Return the absolute cutoff
  !<
  real(DP) function interpolation_kernels_get_cutoff(this)
    implicit none

    type(interpolation_kernels_t), intent(in)  :: this

    ! ---

    if (associated(this%lucy)) &
         interpolation_kernels_get_cutoff  = get_cutoff(this%lucy)
    if (associated(this%square)) &
         interpolation_kernels_get_cutoff  = get_cutoff(this%square)

  endfunction interpolation_kernels_get_cutoff


  !>
  !! Value and derivative
  !!
  !! Returns the value and the derivative of the interpolation kernel for
  !! particles \param r apart. Note that \param v gives the interpolation value,
  !! while the derivative \param w is defined as
  !!
  !! \f[
  !!     \nabla W = -\vec{r} w
  !! \f]
  !<
  subroutine interpolation_kernels_value_and_derivative(this, r, v, w)
    implicit none

    type(interpolation_kernels_t), intent(inout)  :: this
    real(DP), intent(in)                          :: r
    real(DP), intent(out)                         :: v
    real(DP), intent(out)                         :: w


    ! ---

    if (associated(this%lucy)) &
         call value_and_derivative(this%lucy, r, v, w)
    if (associated(this%square)) &
         call value_and_derivative(this%square, r, v, w)

  endsubroutine interpolation_kernels_value_and_derivative

endmodule interpolation_kernels
