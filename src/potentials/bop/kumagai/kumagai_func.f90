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
!! All functions specific to this potential
!<

#include "../default_cutoff.f90"

!>
!! Attractive potential: VA(r), dVA(r)
!!
!! Attractive potential: VA(r), dVA(r)
!<
subroutine VA(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ijpot
  real(DP), intent(in)         :: dr
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  real(DP)  :: expval

  ! ---

  expval  = exp(-this%db%lambda2(ijpot)*dr)
  val     = -this%db%B(ijpot)*expval
  dval    = this%db%B(ijpot)*this%db%lambda2(ijpot)*expval

endsubroutine VA


!>
!! Repulsive potential: VA(r), dVA(r)
!!
!! Repulsive potential: VA(r), dVA(r)
!<
subroutine VR(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ijpot
  real(DP), intent(in)         :: dr
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  real(DP)  :: expval

  ! ---

  expval  = exp(-this%db%lambda1(ijpot)*dr)
  val     = this%db%A(ijpot)*expval
  dval    = -this%db%A(ijpot)*this%db%lambda1(ijpot)*expval

endsubroutine VR


!>
!! Angular contribution to the bond order
!!
!! Angular contribution to the bond order: g(cos(theta)), dg(cos(theta))
!<
subroutine g(this, ktypj, ktypi, ktypk, ijpot, ikpot, costh, val, dval_dcosth)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ktypj
  integer, intent(in)          :: ktypi
  integer, intent(in)          :: ktypk
  integer, intent(in)          :: ijpot
  integer, intent(in)          :: ikpot
  real(DP), intent(in)         :: costh
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval_dcosth

  ! ---

  real(DP)  :: c1, c2, c3, c4, c5, h, h_cos, h_cos_sq, go, ga1, tmp

  ! ---

  c1          = this%db%c1(ktypi)
  c2          = this%db%c2(ktypi)
  c3          = this%db%c3(ktypi)
  c4          = this%db%c4(ktypi)
  c5          = this%db%c5(ktypi)
  h           = this%db%h(ktypi)

  h_cos       = h-costh
  h_cos_sq    = h_cos*h_cos

  tmp         = h_cos/(c3+h_cos_sq)
  go          = c2*tmp
  ga1         = c4*exp(-c5*h_cos_sq)

  val         = go*(1.0_DP+ga1)
  dval_dcosth = -2*(1.0_DP-h_cos*tmp)*val+2*c5*h_cos_sq*go*ga1
  val         = c1+h_cos*val

! ok
!  val  = c1+h_cos*go
!  dval_dcosth = -2*(1.0_DP-h_cos*tmp)*go

! ok
!  val = c1+1.0_DP+ga1
!  dval_dcosth = 2*c5*h_cos*ga1

endsubroutine g


!>
!! Bond order function
!!
!! Determines how the bond-order is computed from zij
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

  real(DP) :: arg, eta, delta

  ! ---

  if (zij > 0.0_DP) then
     eta   = this%db%eta(ktypi)
     delta = -this%db%delta(ktypi)

     arg   = 1.0_DP + zij ** eta
     bij   = arg ** delta
     dfbij = 0.5_DP*fcij*faij* &
          eta*zij**(eta-1.0_DP)* &
          delta*arg**(delta-1.0_DP)
  else
     bij   = 0.0_DP
     dfbij = 0.0_DP
  endif

endsubroutine bo


!>
!! Length dependent contribution to the bond order: h(dr), dh(dr)
!!
!! Length dependent contribution to the bond order: h(dr), dh(dr)
!<
subroutine h(this, ktypj, ktypi, ktypk, ijpot, ikpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ktypj
  integer, intent(in)          :: ktypi
  integer, intent(in)          :: ktypk
  integer, intent(in)          :: ijpot
  integer, intent(in)          :: ikpot
  real(DP), intent(in)         :: dr
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  real(DP)  :: alpha
  integer   :: beta

  ! ---

  alpha = this%db%alpha(ikpot)

  if (alpha == 0.0_DP) then
     val  = 1.0_DP
     dval = 0.0_DP
  else
     beta = this%db%beta(ikpot)

     if (beta == 1) then
        val  = exp(alpha*dr)
        dval = alpha*val
     else
        if (beta == 3) then
           val  = exp(dr*dr*dr)
           dval = 3*alpha * dr*dr * val 
        else
           val  = exp(alpha*dr**beta)
           dval = beta*alpha * dr**(beta-1) * val
        endif
     endif
  endif

endsubroutine h


!>
!! Translation of pairs to pair indices
!!
!! Generate a unique index for the pair \p ktypi \p ktypj
!! of elements
!<
function Z2pair(this, ktypi, ktypj)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ktypi
  integer, intent(in)          :: ktypj
  integer                      :: Z2pair

  ! ---

  Z2pair = PAIR_INDEX(ktypi, ktypj, this%db%nel)

endfunction Z2pair
