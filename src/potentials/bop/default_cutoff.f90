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
!! Cut-off function: fCx(r), dfCccx(r)
!!
!! Cut-off function: fCx(r), dfCccx(r)
!<
subroutine fCin(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in) :: this
  integer, intent(in)        :: ijpot
  real(DP), intent(in)       :: dr
  real(DP), intent(out)      :: val
  real(DP), intent(out)      :: dval

  ! ---

  call fc(this%cut_in(ijpot), dr, val, dval)

endsubroutine fCin


#ifdef SCREENING

!>
!! Cut-off function: fCx(r), dfCccx(r)
!!
!! Cut-off function: fCx(r), dfCccx(r)
!<
subroutine fCar(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in) :: this
  integer, intent(in)        :: ijpot
  real(DP), intent(in)       :: dr
  real(DP), intent(out)      :: val
  real(DP), intent(out)      :: dval

  ! ---

  call fc(this%cut_out(ijpot), dr, val, dval)

endsubroutine fCar


!>
!! Cut-off function: fCx(r), dfCccx(r)
!!
!! Cut-off function: fCx(r), dfCccx(r)
!<
subroutine fCbo(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in) :: this
  integer, intent(in)        :: ijpot
  real(DP), intent(in)       :: dr
  real(DP), intent(out)      :: val
  real(DP), intent(out)      :: dval

  ! ---

  call fc(this%cut_bo(ijpot), dr, val, dval)

endsubroutine fCbo

#endif
