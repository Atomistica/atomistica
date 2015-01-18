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

!**********************************************************************
! Functions for the bond-order potential, i.e. attractive, repulsive
! parts, etc.
!**********************************************************************


!**********************************************************************
! Conjugation counter
!**********************************************************************
elemental subroutine fconj(this, x, fx, dfx)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  real(DP), intent(in)        :: x
  real(DP), intent(out)       :: fx
  real(DP), intent(out)       :: dfx

  ! ---

  real(DP)  :: arg

  ! ---

  if ( x .le. 2.0_DP ) then
     fx   = 1.0_DP
     dfx  = 0.0_DP
  else if( x .ge. 3.0_DP ) then
     fx   = 0.0_DP
     dfx  = 0.0_DP
  else
     arg  = pi  * ( x - 2.0_DP )
     fx   = 0.5_DP * ( 1.0_DP + cos( arg ) )
     dfx  =-0.5_DP * pi * sin( arg )
  endif

endsubroutine fconj


!**********************************************************************
! Cut-off function: fCin(r), dfCin(r)
!**********************************************************************
subroutine fCin(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  if (dr > this%cut_in_h(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else if (dr < this%cut_in_l(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
  else
     call f_and_df(this%spl_fCin(ijpot), dr, val, dval)
  endif

endsubroutine fCin

#ifdef SCREENING

!**********************************************************************
! Cut-off function: fCar(r), dfCar(r)
!**********************************************************************
subroutine fCar(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  if (dr > this%cut_ar_h(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else if (dr < this%cut_ar_l(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
  else
     call f_and_df(this%spl_fCar(ijpot), dr, val, dval)
  endif

endsubroutine fCar


!**********************************************************************
! Cut-off function: fCbo(r), dfCbo(r)
!**********************************************************************
subroutine fCbo(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  if (dr > this%cut_bo_h(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else if (dr < this%cut_bo_l(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
  else
     call f_and_df(this%spl_fCbo(ijpot), dr, val, dval)
  endif

endsubroutine fCbo


!**********************************************************************
! Cut-off function: fCnc(r), dfCnc(r)
!**********************************************************************
subroutine fCnc(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  if (dr > this%cut_nc_h(ijpot)) then
     val   = 0.0_DP
     dval  = 0.0_DP
  else if (dr < this%cut_nc_l(ijpot)) then
     val   = 1.0_DP
     dval  = 0.0_DP
  else
     call f_and_df(this%spl_fCnc(ijpot), dr, val, dval)
  endif

endsubroutine fCnc

#endif

!**********************************************************************
! Attractive potential: VA(r), dVA(r)
!**********************************************************************
subroutine VA(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  real(DP)  :: exp1, exp2, exp3

  ! --- 

  if (dr > this%spl_VA(ijpot)%x0 .and. dr < this%spl_VA(ijpot)%cut) then

     call f_and_df(this%spl_VA(ijpot), dr, val, dval)

  else

     if (ijpot == C_C) then

        exp1 = this%cc_B1*exp(-this%cc_beta1*dr)
        exp2 = this%cc_B2*exp(-this%cc_beta2*dr)
        exp3 = this%cc_B3*exp(-this%cc_beta3*dr)

        val  = - ( exp1 + exp2 + exp3 )
        dval = - ( -this%cc_beta1*exp1 - this%cc_beta2*exp2 - this%cc_beta3*exp3 )

     else if (ijpot == C_H) then

        exp1 = this%ch_B1*exp(-this%ch_beta1*dr)

        val  = - exp1
        dval =   this%ch_beta1*exp1

     else ! if (ijpot == H_H) then                                                                                                                                                    
        exp1 = this%hh_B1*exp(-this%hh_beta1*dr)

        val  = - exp1
        dval =   this%hh_beta1*exp1

     endif

  endif

endsubroutine VA


!**********************************************************************
! Repulsive potential: VA(r), dVA(r)
!**********************************************************************
subroutine VR(this, ijpot, dr, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: dr
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval

  ! ---

  real(DP)  :: exp1, hlp1

  ! ---

  if (dr > this%spl_VR(ijpot)%x0 .and. dr < this%spl_VR(ijpot)%cut) then

     call f_and_df(this%spl_VR(ijpot), dr, val, dval)

  else

     if (ijpot == C_C) then

        exp1 = this%cc_A*exp(-this%cc_alpha*dr)
        hlp1 = 1+this%cc_Q/dr

        val  = hlp1*exp1
        dval = ( -this%cc_Q/(dr**2) - hlp1*this%cc_alpha ) * exp1

     else if (ijpot == C_H) then

        exp1 = this%ch_A*exp(-this%ch_alpha*dr)
        hlp1 = 1+this%ch_Q/dr

        val  = hlp1*exp1
        dval = ( -this%ch_Q/(dr**2) - hlp1*this%ch_alpha ) * exp1

     else ! if (ijpot == H_H) then

        exp1 = this%hh_A*exp(-this%hh_alpha*dr)
        hlp1 = 1+this%hh_Q/dr

        val  = hlp1*exp1
        dval = ( -this%hh_Q/(dr**2) - hlp1*this%hh_alpha ) * exp1

     endif

  endif

endsubroutine VR


!**********************************************************************
! Angular contribution to the bond order: g(cos(theta)), dg(cos(theta))
!**********************************************************************
elemental subroutine g(this, ktyp, costh, n, val, dval_dcosth, dval_dN)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ktyp
  real(DP), intent(in)        :: costh
  real(DP), intent(in)        :: n
  real(DP), intent(out)       :: val
  real(DP), intent(out)       :: dval_dcosth
  real(DP), intent(out)       :: dval_dN

  ! ---

  real(DP)  :: v1, v2, dv1, dv2, s, ds, arg
  integer   :: i, ig

  ! ---

  if (ktyp == rebo2_C_) then

     if (n < 3.2_DP) then

        call cc_g_from_spline(this, this%cc_g2_coeff, costh, val, dval_dcosth)
        dval_dN = 0.0_DP

     else if (n > 3.7) then

        call cc_g_from_spline(this, this%cc_g1_coeff, costh, val, dval_dcosth)
        dval_dN = 0.0_DP

     else

        call cc_g_from_spline(this, this%cc_g1_coeff, costh, v1, dv1)
        call cc_g_from_spline(this, this%cc_g2_coeff, costh, v2, dv2)

        arg = 2*PI*(n-3.2_DP)
        s  = (1+cos(arg))/2
        ds = -PI*sin(arg)

        val          =  v1*(1-s) +  v2*s
        dval_dcosth  = dv1*(1-s) + dv2*s
        dval_dN      = (v2-v1)*ds

     endif

  else

     ig = this%igh(int(-costh*12.0D0)+13)

     val          = this%spgh(1, ig) + this%spgh(2, ig)*costh
     dval_dcosth  = this%spgh(2, ig)
     do i = 3, 6
        val          = val  + this%spgh(i, ig)*costh**(i-1)
        dval_dcosth  = dval_dcosth + (i-1)*this%spgh(i, ig)*costh**(i-2)
     enddo

  endif

endsubroutine g


!**********************************************************************
! Angular contribution to the bond order: g(cos(theta)), dg(cos(theta))
!**********************************************************************
elemental subroutine cc_g_from_spline(this, g_coeff, costh, val, dval)
  implicit none

  type(BOP_TYPE), intent(in)   :: this
  type(g_coeff_t), intent(in)  :: g_coeff
  real(DP), intent(in)         :: costh
  real(DP), intent(out)        :: val
  real(DP), intent(out)        :: dval

  ! ---

  integer   :: i, j
  real(DP)  :: h, dh

  ! ---

  !    if (costh >= cc_g_theta(1)) then

  if (costh < this%cc_g_theta(2)) then
     j = 1
  else if (costh < this%cc_g_theta(3)) then
     j = 2
  else ! if (costh <= cc_g_theta(6)) then
     j = 3
     !       else
     !          write (*, '(A,ES20.10,A)')  "[g] value ", costh, " outside the region for which the spline is defined."
     !          stop
  endif

  h   = g_coeff%c(1, j) + g_coeff%c(2, j)*costh
  dh  = g_coeff%c(2, j)
  do i = 3, 6
     h   = h  + g_coeff%c(i, j)*costh**(i-1)
     dh  = dh + (i-1)*g_coeff%c(i, j)*costh**(i-2)
  enddo

  val   = h
  dval  = dh

  !    else
  !       write (*, '(A,ES20.10,A)')  "[g] value ", costh, " outside the region for which the spline is defined."
  !       stop
  !    endif

endsubroutine cc_g_from_spline


!**********************************************************************
! Bond order function
!**********************************************************************
elemental subroutine bo(this, ktypi, ijpot, zij, fcij, faij, bij, dfbij)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ktypi
  integer, intent(in)         :: ijpot
  real(DP), intent(in)        :: zij
  real(DP), intent(in)        :: fcij
  real(DP), intent(in)        :: faij
  real(DP), intent(out)       :: bij
  real(DP), intent(out)       :: dfbij

  ! ---

  real(DP) :: arg

  ! ---

  arg    = 1.0 + zij
  bij    = arg ** this%conpe(ktypi)
  dfbij  = this%conan(ktypi) * fcij * faij * arg ** this%conpf(ktypi)

endsubroutine bo


!**********************************************************************
! Length dependent contribution to the bond order: h(dr), dh(dr)
!**********************************************************************
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

  if ( (ijpot+ikpot) <= 4 ) then
     val  = 1.0_DP
     dval = 0.0_DP
  else

     !
     ! exp( alpha * ( rij - rik ) ) & dexp / d(rij-rik).
     !

     val  = this%conear(ijpot, ikpot) * exp( this%conalp * dr )
     dval = this%conalp * val

     !       write (*, '(2I5,5F20.10)')  ijpot, ikpot, dr, val, dval, conear(ijpot, ikpot), conalp
  endif

endsubroutine h


!**********************************************************************
! Generate an index for this *pair* if elements
!**********************************************************************
elemental function Z2pair(this, ktypi, ktypj)
  implicit none

  type(BOP_TYPE), intent(in)  :: this
  integer, intent(in)         :: ktypi
  integer, intent(in)         :: ktypj
  integer                     :: Z2pair

  ! ---

  if (ktypi == rebo2_C_) then
     Z2pair = ktypj
  else if (ktypj == rebo2_C_) then
     Z2pair = ktypi
  else
     Z2pair = ktypi + ktypj
  endif

endfunction Z2pair
