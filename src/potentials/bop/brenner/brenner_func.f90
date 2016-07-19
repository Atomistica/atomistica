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

#include "../default_cutoff.f90"

!>
!! Attractive potential: VA(r), dVA(r)
!<
elemental subroutine VA(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  real(DP)  :: expval

  ! ---

  expval  = exp(-this%expA(ijpot)*(dr-this%db%r0(ijpot)))
  val     = -this%VA_f(ijpot)*expval
  dval    = this%VA_f(ijpot)*this%expA(ijpot)*expval

endsubroutine VA


!>
!! Repulsive potential: VA(r), dVA(r)
!<
elemental subroutine VR(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  real(DP)  :: expval

  ! ---

  expval  = exp(-this%expR(ijpot)*(dr-this%db%r0(ijpot)))
  val     = this%VR_f(ijpot)*expval
  dval    = -this%VR_f(ijpot)*this%expR(ijpot)*expval

endsubroutine VR


!>
!! Angular contribution to the bond order: g(cos(theta)), dg(cos(theta))
!<
elemental subroutine g(this, ktypj, ktypi, ktypk, ijpot, ikpot, costh, val, dval_dcosth)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ktypj
  integer, intent(in)         :: ktypi
  integer, intent(in)         :: ktypk
  integer, intent(in)         :: ijpot
  integer, intent(in)         :: ikpot
  real(DP), intent(in)        :: costh
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval_dcosth

  ! ---

  real(DP)  :: h

  ! ---


  h           = this%d_sq(ikpot)+(this%db%h(ikpot)+costh)**2
  val         = this%db%gamma(ikpot)*(1+this%c_d(ikpot)-this%c_sq(ikpot)/h)
  dval_dcosth = 2*this%db%gamma(ikpot)*this%c_sq(ikpot)*(this%db%h(ikpot)+costh)/(h**2)

endsubroutine g


!>
!! Bond order function
!<
subroutine bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer,        intent(in)  :: ktypi
  integer,        intent(in)  :: ijpot
  real(DP),       intent(in)  :: zij
  real(DP),       intent(in)  :: fcij
  real(DP),       intent(in)  :: faij
  real(DP),       intent(out) :: bij
  real(DP),       intent(out) :: dfbij

  ! ---

  real(DP) :: arg

  ! ---

  if (this%db%n(ijpot) == 1.0_DP) then

     arg   = 1.0_DP + zij
     bij   = arg ** this%bo_exp(ijpot)
     dfbij = this%bo_fac(ijpot) * fcij * faij * arg ** this%bo_exp1(ijpot)

  else

     if (zij > 0.0_DP) then
        arg   = 1.0_DP + zij ** this%db%n(ijpot)
        bij   = arg ** this%bo_exp(ijpot)
        dfbij = &
             this%bo_fac(ijpot) * fcij * faij &
             * zij ** ( this%db%n(ijpot) - 1.0_DP ) &
             * arg ** this%bo_exp1(ijpot)
     else
        bij   = 0.0_DP
        dfbij = 0.0_DP
     endif

  endif

endsubroutine bo


!>
!! Length dependent contribution to the bond order: h(dr), dh(dr)
!<
elemental subroutine h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ktypj
  integer, intent(in)         :: ktypi
  integer, intent(in)         :: ktypk
  integer, intent(in)         :: ijpot
  integer, intent(in)         :: ikpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  real(DP)  :: mu, arg
  integer   :: m

  ! ---

  mu = this%db%mu(ikpot)

  if (mu == 0.0_DP) then
     val  = 1.0_DP
     dval = 0.0_DP
  else
     m = this%db%m(ikpot)

     if (m == 1) then
        val  = exp(2*mu*dr)
        dval = 2*mu*val
     else
        if (m == 3) then
           arg  = 2*mu*dr
           val  = exp(arg*arg*arg)
           dval = 6*mu * arg*arg * val 
        else
           val  = exp((2*mu*dr)**m)
           dval = 2*mu*m * (2*mu*dr)**(m-1) * val
        endif
     endif
  endif

endsubroutine h


!>
!! Generate an index for this *pair* if elements
!<
elemental function Z2pair(this, i, j)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: i
  integer, intent(in)         :: j
  integer                     :: Z2pair

  ! ---

  Z2pair = PAIR_INDEX(i, j, this%db%nel)

endfunction Z2pair
  
