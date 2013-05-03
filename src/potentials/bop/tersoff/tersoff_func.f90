!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================
!>
!! All functions specific to this potential
!<


!>
!! Cut-off function: fCx(r), dfCccx(r)
!!
!! Cut-off function: fCx(r), dfCccx(r)
!<
subroutine fCin(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ijpot
  real(DP), intent(in)         :: dr
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  real(DP)  :: fac, arg
#ifdef EXP_CUTOFF
  real(DP)  :: arg2
#endif

  ! ---

  if (dr < this%db%r1(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
#ifdef EXP_CUTOFF
  else
     fac  = 2.0_DP/( this%db%r2(ijpot) - this%db%r1(ijpot) )
     arg  = fac*( dr-this%db%r1(ijpot) )
     arg2 = arg*arg
     val  = exp(-arg*arg2)
     dval = -3*fac*arg2*val
  endif
#else
  else if (dr > this%db%r2(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else
     fac   = PI/( this%db%r2(ijpot) - this%db%r1(ijpot) )
     arg   = fac*( dr-this%db%r1(ijpot) )
     val   = 0.5_DP *  ( 1.0_DP + cos( arg ) )
     dval  = -0.5_DP * fac * sin( arg )
  endif
#endif

endsubroutine fCin


#ifdef SCREENING

!>
!! Cut-off function: fCx(r), dfCccx(r)
!!
!! Cut-off function: fCx(r), dfCccx(r)
!<
subroutine fCar(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ijpot
  real(DP), intent(in)         :: dr
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  real(DP)  :: fac, arg
#ifdef EXP_CUTOFF
  real(DP)  :: arg2
#endif

  ! ---

  if (dr < this%db%or1(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
#ifdef EXP_CUTOFF
  else
     fac  = 2.0_DP/( this%db%or2(ijpot) - this%db%or1(ijpot) )
     arg  = fac*( dr-this%db%or1(ijpot) )
     arg2 = arg*arg
     val  = exp(-arg*arg2)
     dval = -3*fac*arg2*val
  endif
#else
  else if (dr > this%db%or2(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else
     fac   = PI/( this%db%or2(ijpot) - this%db%or1(ijpot) )
     arg   = fac*( dr-this%db%or1(ijpot) )
     val   = 0.5_DP *  ( 1.0_DP + cos( arg ) )
     dval  = -0.5_DP * fac * sin( arg )
  endif
#endif

endsubroutine fCar


!>
!! Cut-off function: fCx(r), dfCccx(r)
!!
!! Cut-off function: fCx(r), dfCccx(r)
!<
subroutine fCbo(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ijpot
  real(DP), intent(in)         :: dr
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  real(DP)  :: fac, arg
#ifdef EXP_CUTOFF
  real(DP)  :: arg2
#endif

  ! ---

  if (dr < this%db%bor1(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
#ifdef EXP_CUTOFF
  else
     fac  = 2.0_DP/( this%db%bor2(ijpot) - this%db%bor1(ijpot) )
     arg  = fac*( dr-this%db%bor1(ijpot) )
     arg2 = arg*arg
     val  = exp(-arg*arg2)
     dval = -3*fac*arg2*val
  endif
#else
  else if (dr > this%db%bor2(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else
     fac   = PI/( this%db%bor2(ijpot) - this%db%bor1(ijpot) )
     arg   = fac*( dr-this%db%bor1(ijpot) )
     val   = 0.5_DP *  ( 1.0_DP + cos( arg ) )
     dval  = -0.5_DP * fac * sin( arg )
  endif
#endif

endsubroutine fCbo

#endif


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

  real(DP)  :: h, c_sq, d_sq, h_c

  ! ---

  h_c          = this%db%h(ktypi) - costh
  c_sq         = this%db%c(ktypi)**2
  d_sq         = this%db%d(ktypi)**2

  h            = d_sq + h_c**2
  val          = 1.0_DP + c_sq/d_sq - c_sq/h
  dval_dcosth  = -2*c_sq*h_c/(h**2)

endsubroutine g


!>
!! Bond order function
!!
!! Determines how the bond-order is computed from zij
!<
subroutine bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)          :: ktypi
  integer, intent(in)          :: ijpot
  real(DP), intent(in)         :: zij
  real(DP), intent(in)         :: fcij
  real(DP), intent(in)         :: faij
  real(DP), intent(out)        :: bij
  real(DP), intent(out)        :: dfbij

  ! ---

  real(DP) :: arg, e, b

  ! ---

  e      = -0.5_DP/this%db%n(ktypi)
  b      = this%db%beta(ktypi) ** this%db%n(ktypi)

  arg    = 1.0_DP + b * zij ** this%db%n(ktypi)
  bij    = this%db%xi(ijpot) * arg ** e

  if (zij /= 0.0_DP) then

     dfbij  = &
          -0.25_DP * fcij * faij * this%db%xi(ijpot) * b &
          * zij ** ( this%db%n(ktypi) - 1.0_DP ) &
          * arg ** ( e - 1.0_DP )

  else

     dfbij  = 0.0_DP

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
