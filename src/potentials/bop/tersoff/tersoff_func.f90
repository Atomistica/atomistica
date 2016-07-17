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

  expval  = exp(-this%db%mu(ijpot)*dr)
  val     = -this%db%B(ijpot)*expval
  dval    = this%db%B(ijpot)*this%db%mu(ijpot)*expval

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

  expval  = exp(-this%db%lambda(ijpot)*dr)
  val     = this%db%A(ijpot)*expval
  dval    = -this%db%A(ijpot)*this%db%lambda(ijpot)*expval

endsubroutine VR


!>
!! Angular contribution to the bond order: g(cos(theta)), dg(cos(theta))
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

  real(DP)  :: omega, h, c_sq, d_sq, h_c

  ! ---

  omega        = this%db%omega(ikpot)
  h_c          = this%db%h(ktypi) - costh
  c_sq         = this%db%c(ktypi)**2
  d_sq         = this%db%d(ktypi)**2

  h            = d_sq + h_c**2
  val          = omega*(1.0_DP + c_sq/d_sq - c_sq/h)
  dval_dcosth  = -2*omega*c_sq*h_c/(h**2)

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

  real(DP) :: arg, e, b

  ! ---

  if (zij > 0.0_DP) then

     e     = -0.5_DP/this%db%n(ktypi)
     b     = this%db%beta(ktypi) ** this%db%n(ktypi)

     arg   = 1.0_DP + b * zij ** this%db%n(ktypi)
     bij   = this%db%xi(ijpot) * arg ** e

     dfbij = &
          -0.25_DP * fcij * faij * this%db%xi(ijpot) * b &
          * zij ** ( this%db%n(ktypi) - 1.0_DP ) &
          * arg ** ( e - 1.0_DP )

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

  real(DP)  :: mu, arg
  integer   :: m

  ! ---

  mu = this%db%mubo(ikpot)

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
           dval = 2*mu*m * arg*arg * val 
        else
           val  = exp((2*mu*dr)**m)
           dval = 2*mu*m * (2*mu*dr)**(m-1) * val
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
