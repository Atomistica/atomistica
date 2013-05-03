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
