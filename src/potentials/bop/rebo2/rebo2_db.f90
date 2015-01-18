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
! Database containing materials parameters for 2nd generation REBO
!**********************************************************************

  !**********************************************************************
  ! Initialize the REBO with default parameters
  !********************************************************************** 
  subroutine rebo2_db_init(this)
    implicit none

    type(BOP_TYPE)  :: this

    ! ---

    real(DP)  :: in_Fcc(0:4, 0:4, 0:9)
    real(DP)  :: in_dFdi(0:4, 0:4, 0:9)
    real(DP)  :: in_dFdj(0:4, 0:4, 0:9)
    real(DP)  :: in_dFdk(0:4, 0:4, 0:9)

    real(DP)  :: in_Fch(0:4, 0:4, 0:9)
    real(DP)  :: in_Fhh(0:4, 0:4, 0:9)

    real(DP)  :: in_Pcc(0:5, 0:5)
    real(DP)  :: in_Pch(0:5, 0:5)
    real(DP)  :: in_Tcc(0:4, 0:4, 0:9)

    ! ---

#ifdef ZERO_TABLES
    in_Fcc   = 0.0_DP
    in_dFdi  = 0.0_DP
    in_dFdj  = 0.0_DP
    in_dFdk  = 0.0_DP
    in_Fch   = 0.0_DP
    in_Fhh   = 0.0_DP
    in_Pcc   = 0.0_DP
    in_Pch   = 0.0_DP
    in_Tcc   = 0.0_DP
#else
    call rebo2_default_Fcc_table(in_Fcc, in_dFdi, in_dFdj, in_dFdk)
    call rebo2_default_Fch_table(in_Fch)
    call rebo2_default_Fhh_table(in_Fhh)
    call rebo2_default_Pcc_table(in_Pcc)
    call rebo2_default_Pch_table(in_Pch)
    call rebo2_default_Tcc_table(in_Tcc)
#endif

    call rebo2_db_init_with_parameters( &
         this, in_Fcc, in_dFdi, in_dFdj, in_dFdk, in_Fch, in_Fhh, in_Pcc, &
         in_Pch, in_Tcc)

  endsubroutine rebo2_db_init


  !**********************************************************************
  ! Initialize the REBO with a set of chosen parameters.
  ! Note: General parameters are taken from this and need to be set
  ! before call to this init routine.
  !********************************************************************** 
  subroutine rebo2_db_init_with_parameters( &
       this, in_Fcc, in_dFdi, in_dFdj, in_dFdk, in_Fch, in_Fhh, in_Pcc, in_Pch, in_Tcc)
    implicit none

    type(BOP_TYPE), intent(inout)  :: this

    real(DP), intent(in)          :: in_Fcc(0:4, 0:4, 0:9)
    real(DP), intent(in)          :: in_dFdi(0:4, 0:4, 0:9)
    real(DP), intent(in)          :: in_dFdj(0:4, 0:4, 0:9)
    real(DP), intent(in)          :: in_dFdk(0:4, 0:4, 0:9)

    real(DP), intent(in)          :: in_Fch(0:4, 0:4, 0:9)
    real(DP), intent(in)          :: in_Fhh(0:4, 0:4, 0:9)

    real(DP), intent(in)          :: in_Pcc(0:5, 0:5)
    real(DP), intent(in)          :: in_Pch(0:5, 0:5)
    real(DP), intent(in)          :: in_Tcc(0:4, 0:4, 0:9)

    ! ---

    integer  :: i

    ! ---

#ifdef SCREENING
    this%dC        = this%Cmax-this%Cmin
    this%C_dr_cut  = this%Cmax**2/(4*(this%Cmax-1))
#endif

    if (ilog /= -1) then
       write (ilog, '(A)')   "- rebo2_db_init -"

#ifdef SCREENING
       write (ilog, '(5X,A,F20.10)')  "C_min     = ", this%Cmin
       write (ilog, '(5X,A,F20.10)')  "C_max     = ", this%Cmax
       write (ilog, '(5X,A,F20.10)')  "dC        = ", this%dC
       write (ilog, '(5X,A,F20.10)')  "C_dr_cut  = ", this%C_dr_cut
       write (ilog, '(5X,A,F20.10)')  "cc_in_r1  = ", this%cc_in_r1
       write (ilog, '(5X,A,F20.10)')  "cc_in_r2  = ", this%cc_in_r2
       write (ilog, '(5X,A,F20.10)')  "cc_ar_r1  = ", this%cc_ar_r1
       write (ilog, '(5X,A,F20.10)')  "cc_ar_r2  = ", this%cc_ar_r2
       write (ilog, '(5X,A,F20.10)')  "cc_bo_r1  = ", this%cc_bo_r1
       write (ilog, '(5X,A,F20.10)')  "cc_bo_r2  = ", this%cc_bo_r2
       write (ilog, '(5X,A,F20.10)')  "cc_nc_r1  = ", this%cc_nc_r1
       write (ilog, '(5X,A,F20.10)')  "cc_nc_r2  = ", this%cc_nc_r2
#else
       write (ilog, '(5X,A,F20.10)')  "cc_r1     = ", this%cc_in_r1
       write (ilog, '(5X,A,F20.10)')  "cc_r2     = ", this%cc_in_r2
#endif
       write (ilog, '(5X,A,L1)')      "dihedral  = ", this%with_dihedral

#ifdef ZERO_TABLES
       write (ilog, '(5X,A)')  "Warning: All tables are set to zero!"
#endif

#ifdef EXP_CUT
       write (ilog, '(5X,A)')  "Using exponential cut-off function."
#endif

       write (ilog, *)
    endif

    !
    ! bond order constants.
    !

    this%conpe(1) = - 0.5
    this%conpe(3) = - 0.5
    this%conan(1) = 0.5 * this%conpe(1)
    this%conan(3) = 0.5 * this%conpe(3)
    this%conpf(1) = this%conpe(1) - 1.0
    this%conpf(3) = this%conpe(3) - 1.0

    !
    ! bond order penalty function constants
    !

    this%conalp          = this%hhh_lambda
    this%conear    = 0.0
    this%conear(C_C,C_C) = 1.0
    this%conear(C_C,C_H) = exp( this%conalp * ( this%ch_re - this%cc_re ) )
    this%conear(C_C,H_H) = exp( this%conalp * ( this%hh_re - this%cc_re ) )
    this%conear(C_H,C_C) = 1.0 / this%conear(C_C,C_H)
    this%conear(C_H,C_H) = 1.0
    this%conear(C_H,H_H) = exp( this%conalp * ( this%hh_re - this%ch_re ) )
    this%conear(H_H,C_C) = 1.0 / this%conear(C_C,H_H)
    this%conear(H_H,C_H) = 1.0 / this%conear(C_H,H_H)
    this%conear(H_H,H_H) = 1.0

    !
    ! cutoff constants.
    !

    this%cut_in_l(:)     = 0.0
    this%cut_in_l(C_C)   = this%cc_in_r1
    this%cut_in_l(C_H)   = this%ch_r1
    this%cut_in_l(H_H)   = this%hh_r1

    this%cut_in_h(:)     = 0.0
    this%cut_in_h(C_C)   = this%cc_in_r2
    this%cut_in_h(C_H)   = this%ch_r2
    this%cut_in_h(H_H)   = this%hh_r2

    this%cut_in_h2(:)    = 0.0
    this%cut_in_h2(C_C)  = this%cc_in_r2 ** 2
    this%cut_in_h2(C_H)  = this%ch_r2 ** 2
    this%cut_in_h2(H_H)  = this%hh_r2 ** 2

    this%cut_in_m(:)     = 0.0
    this%cut_in_m(C_C)   = (this%cc_in_r1+this%cc_in_r2)/2
    this%cut_in_m(C_H)   = (this%ch_r1+this%ch_r2)/2
    this%cut_in_m(H_H)   = (this%hh_r1+this%hh_r2)/2

#ifdef SCREENING
    this%cut_ar_l(:)     = 0.0
    this%cut_bo_l(:)     = 0.0
    this%cut_nc_l(:)     = 0.0
    this%cut_ar_l(C_C)   = this%cc_ar_r1
    this%cut_ar_l(C_H)   = this%ch_r1
    this%cut_ar_l(H_H)   = this%hh_r1
    this%cut_bo_l(C_C)   = this%cc_bo_r1
    this%cut_bo_l(C_H)   = this%ch_r1
    this%cut_bo_l(H_H)   = this%hh_r1
    this%cut_nc_l(C_C)   = this%cc_nc_r1
    this%cut_nc_l(C_H)   = this%ch_r1
    this%cut_nc_l(H_H)   = this%hh_r1

    this%cut_ar_h(:)     = 0.0
    this%cut_bo_h(:)     = 0.0
    this%cut_nc_h(:)     = 0.0
    this%cut_ar_h(C_C)   = this%cc_ar_r2
    this%cut_ar_h(C_H)   = this%ch_r2
    this%cut_ar_h(H_H)   = this%hh_r2
    this%cut_bo_h(C_C)   = this%cc_bo_r2
    this%cut_bo_h(C_H)   = this%ch_r2
    this%cut_bo_h(H_H)   = this%hh_r2
    this%cut_nc_h(C_C)   = this%cc_nc_r2
    this%cut_nc_h(C_H)   = this%ch_r2
    this%cut_nc_h(H_H)   = this%hh_r2

    this%cut_ar_h2(:)    = 0.0
    this%cut_bo_h2(:)    = 0.0
    this%cut_nc_h2(:)    = 0.0
    this%cut_ar_h2(C_C)  = this%cc_ar_r2 ** 2
    this%cut_ar_h2(C_H)  = this%ch_r2 ** 2
    this%cut_ar_h2(H_H)  = this%hh_r2 ** 2
    this%cut_bo_h2(C_C)  = this%cc_bo_r2 ** 2
    this%cut_bo_h2(C_H)  = this%ch_r2 ** 2
    this%cut_bo_h2(H_H)  = this%hh_r2 ** 2
    this%cut_nc_h2(C_C)  = this%cc_nc_r2 ** 2
    this%cut_nc_h2(C_H)  = this%ch_r2 ** 2
    this%cut_nc_h2(H_H)  = this%hh_r2 ** 2

    this%cut_ar_m(:)     = 0.0
    this%cut_bo_m(:)     = 0.0
    this%cut_nc_m(:)     = 0.0
    this%cut_ar_m(C_C)   = (this%cc_ar_r1+this%cc_ar_r2)/2
    this%cut_ar_m(C_H)   = (this%ch_r1+this%ch_r2)/2
    this%cut_ar_m(H_H)   = (this%hh_r1+this%hh_r2)/2
    this%cut_bo_m(C_C)   = (this%cc_bo_r1+this%cc_bo_r2)/2
    this%cut_bo_m(C_H)   = (this%ch_r1+this%ch_r2)/2
    this%cut_bo_m(H_H)   = (this%hh_r1+this%hh_r2)/2
    this%cut_nc_m(C_C)   = (this%cc_nc_r1+this%cc_nc_r2)/2
    this%cut_nc_m(C_H)   = (this%ch_r1+this%ch_r2)/2
    this%cut_nc_m(H_H)   = (this%hh_r2+this%hh_r1)/2

    do i = 1, 10
       this%max_cut_sq(i) = maxval( (/ this%cut_ar_h2(i), this%cut_in_h2(i), this%cut_bo_h2(i), this%cut_nc_h2(i) /) )
    enddo
#else
    do i = 1, 10
       this%max_cut_sq(i) = this%cut_in_h2(i)
    enddo
#endif

    !
    ! Generate the coefficients for
    ! the bi- and tri- cubic interpolation functions.
    !

    call rebo2_db_make_cc_g_spline(this)

    !
    ! Initialize look-up tables
    !

    call init(this%Fcc, 4, 4, 9, in_Fcc, in_dFdi, in_dFdj, in_dFdk)
    call init(this%Fch, 4, 4, 9, in_Fch)
    call init(this%Fhh, 4, 4, 9, in_Fhh)
    call init(this%Pcc, 5, 5, in_Pcc)
    call init(this%Pch, 5, 5, in_Pch)
    call init(this%Tcc, 4, 4, 9, in_Tcc)

    if (ilog /= -1) then
       write (ilog, '(5X,A,F20.10)')  "C-C cut-off = ", sqrt(this%max_cut_sq(C_C))
       write (ilog, '(5X,A,F20.10)')  "C-H cut-off = ", sqrt(this%max_cut_sq(C_H))
       write (ilog, '(5X,A,F20.10)')  "H-H cut-off = ", sqrt(this%max_cut_sq(H_H))

       call prlog("     Fcc:")
       call table3d_prlog(this%Fcc, indent=5)
       call prlog("     Fch:")
       call table3d_prlog(this%Fch, indent=5)
       call prlog("     Fhh:")
       call table3d_prlog(this%Fhh, indent=5)
       call prlog("     Pcc:")
       call table2d_prlog(this%Pcc, indent=5)
       call prlog("     Pch:")
       call table2d_prlog(this%Pch, indent=5)
       call prlog("     Tcc:")
       call table3d_prlog(this%Tcc, indent=5)
       
       write (ilog, *)
    endif
 
    !
    ! Make splines for attractive, repulsive functions
    !

    call rebo2_db_make_splines(this)

    this%tables_allocated = .true.

  endsubroutine rebo2_db_init_with_parameters


  !**********************************************************************
  ! Free all resources
  !********************************************************************** 
  subroutine rebo2_db_del(this)
    implicit none

    type(BOP_TYPE), intent(inout)  :: this

    ! ---

    if (this%neighbor_list_allocated) then
#ifdef NUM_NEIGHBORS
       deallocate(this%nn)
#endif

       deallocate(this%neb_seed)
       deallocate(this%neb_last)

       deallocate(this%neb)
       deallocate(this%nbb)
#ifndef LAMMPS
       deallocate(this%dcell)
#endif
       deallocate(this%bndtyp)
       deallocate(this%bndlen)
       deallocate(this%bndnm)
       deallocate(this%cutfcnar)
       deallocate(this%cutdrvar)

#ifdef SCREENING
       deallocate(this%cutfcnbo)
       deallocate(this%cutdrvbo)
       deallocate(this%cutfcnnc)
       deallocate(this%cutdrvnc)
       deallocate(this%sneb_seed)
       deallocate(this%sneb_last)
       deallocate(this%sneb)
       deallocate(this%sbnd)
       deallocate(this%sfacbo)
       deallocate(this%sfacnc)
       deallocate(this%cutdrarik)
       deallocate(this%cutdrarjk)
       deallocate(this%cutdrboik)
       deallocate(this%cutdrbojk)
       deallocate(this%cutdrncik)
       deallocate(this%cutdrncjk)
#endif

       this%neighbor_list_allocated = .false.
    endif

    if (this%tables_allocated) then
       call del(this%Fcc)
       call del(this%Fch)
       call del(this%Fhh)
       call del(this%Pcc)
       call del(this%Pch)
       call del(this%Tcc)

       call del(this%spl_VA(C_C))
       call del(this%spl_VA(C_H))
       call del(this%spl_VA(H_H))

       call del(this%spl_VR(C_C))
       call del(this%spl_VR(C_H))
       call del(this%spl_VR(H_H))

       call del(this%spl_fCin(C_C))
       call del(this%spl_fCin(C_H))
       call del(this%spl_fCin(H_H))

#ifdef SCREENING
       call del(this%spl_fCar(C_C))
       call del(this%spl_fCar(C_H))
       call del(this%spl_fCar(H_H))

       call del(this%spl_fCbo(C_C))
       call del(this%spl_fCbo(C_H))
       call del(this%spl_fCbo(H_H))

       call del(this%spl_fCnc(C_C))
       call del(this%spl_fCnc(C_H))
       call del(this%spl_fCnc(H_H))
#endif

       this%tables_allocated = .false.
    endif

  endsubroutine rebo2_db_del


  !**********************************************************************
  ! Compute the coefficients for the g(cos(theta)) spline
  ! for C-C interaction
  !********************************************************************** 
  subroutine rebo2_db_make_cc_g_spline(this)
    implicit none

    type(BOP_TYPE), intent(inout)  :: this

    ! ---

    real(DP)  :: A(6, 6), Asave(6, 6)
    real(DP)  :: B(6)

    real(DP)  :: z

    integer   :: i, j, k

    integer   :: ipiv(6)

!    real(DP)  :: h, dh

    ! ---

    !
    ! Third interval
    !

    do i = 3, 6
!       z = (g_theta(i)-g_theta(3))/(g_theta(6)-g_theta(3))
       z = this%cc_g_theta(i)

       do j = 1, 6
          A(i-2, j) = z**(j-1)
       enddo
    enddo

    z = this%cc_g_theta(3)
    A(5, :) = 0.0
    A(6, :) = 0.0
    A(5, 2) = 1.0
    A(6, 3) = 2.0
    do j = 3, 6
                    A(5, j) = (j-1)*z**(j-2)        ! First derivative on left boundary
       if (j >= 4)  A(6, j) = (j-2)*(j-1)*z**(j-3)  ! Second derivative on left boundary
    enddo

    B(1:4) = this%cc_g_g1(3:6)
    B(5)   = this%cc_g_dg1(3)
    B(6)   = this%cc_g_d2g1(3)

!    write (*, *)  1

    Asave = A
    call dgesv(6, 1, A, 6, ipiv, B, 6, i)

    if (i /= 0) then
       write (*, '(A,I5)')  "[rebo2_make_cc_g_spline] dgesv failed. info = ", i
       stop
    endif

    this%cc_g1_coeff%c(:, 3) = B(:)

    B(1:4) = this%cc_g_g2(3:6)
    B(5)   = this%cc_g_dg1(3)
    B(6)   = this%cc_g_d2g1(3)

!    write (*, *)  2

    A = Asave
    call dgesv(6, 1, A, 6, ipiv, B, 6, i)

    if (i /= 0) then
       write (*, '(A,I5)')  "[rebo2_make_cc_g_spline] dgesv failed. info = ", i
       stop
    endif

    this%cc_g2_coeff%c(:, 3) = B(:)


    !
    ! First interval and second interval
    !

    do k = 0, 1

       A = 0.0_DP

       do i = 0, 1
          z = this%cc_g_theta(1+k)*(1-i) + this%cc_g_theta(2+k)*i

          A(3*i+1, 1) = 1.0
          A(3*i+2, 2) = 1.0
          A(3*i+3, 3) = 2.0
          do j = 2, 6
                           A(3*i+1, j) = z**(j-1)
             if (j >= 3)   A(3*i+2, j) = (j-1)*z**(j-2)
             if (j >= 4)   A(3*i+3, j) = (j-2)*(j-1)*z**(j-3)
          enddo
       enddo

       B(1)   = this%cc_g_g1(1+k)
       B(2)   = this%cc_g_dg1(1+k)
       B(3)   = this%cc_g_d2g1(1+k)
       B(4)   = this%cc_g_g1(2+k)
       B(5)   = this%cc_g_dg1(2+k)
       B(6)   = this%cc_g_d2g1(2+k)

!       write (*, *)  3

       call dgesv(6, 1, A, 6, ipiv, B, 6, i)

       if (i /= 0) then
          write (*, '(A,I5)')  "[rebo2_make_cc_g_spline] dgesv failed. info = ", i
          stop
       endif

       this%cc_g1_coeff%c(:, 1+k) = B(:)
       this%cc_g2_coeff%c(:, 1+k) = B(:)

    enddo

  endsubroutine rebo2_db_make_cc_g_spline


! --- Functions ---

  function cc_VA(dr, cc_B1, cc_B2, cc_B3, cc_beta1, cc_beta2, cc_beta3) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: cc_B1
    real(DP), intent(in)  :: cc_B2
    real(DP), intent(in)  :: cc_B3
    real(DP), intent(in)  :: cc_beta1
    real(DP), intent(in)  :: cc_beta2
    real(DP), intent(in)  :: cc_beta3
    real(DP)              :: val

    ! ---

    real(DP)  :: exp1, exp2, exp3

    ! ---

    exp1 = cc_B1*exp(-cc_beta1*dr)
    exp2 = cc_B2*exp(-cc_beta2*dr)
    exp3 = cc_B3*exp(-cc_beta3*dr)

    val  = - ( exp1 + exp2 + exp3 )

  endfunction cc_VA


  function cc_VR(dr, cc_A, cc_Q, cc_alpha) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: cc_A
    real(DP), intent(in)  :: cc_Q
    real(DP), intent(in)  :: cc_alpha
    real(DP)              :: val

    ! ---

    real(DP)  :: exp1, hlp1

    ! ---

    exp1 = cc_A*exp(-cc_alpha*dr)
    hlp1 = 1+cc_Q/dr

    val  = hlp1*exp1

  endfunction cc_VR


  function ch_VA(dr, ch_B1, ch_beta1) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: ch_B1
    real(DP), intent(in)  :: ch_beta1
    real(DP)              :: val

    ! ---

    real(DP)  :: exp1

    ! ---

    exp1 = ch_B1*exp(-ch_beta1*dr)

    val  = - exp1

  endfunction ch_VA


  function ch_VR(dr, ch_A, ch_Q, ch_alpha) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: ch_A
    real(DP), intent(in)  :: ch_Q
    real(DP), intent(in)  :: ch_alpha
    real(DP)              :: val

    ! ---

    real(DP)  :: exp1, hlp1

    ! ---

    exp1 = ch_A*exp(-ch_alpha*dr)
    hlp1 = 1+ch_Q/dr

    val  = hlp1*exp1

  endfunction ch_VR


  function hh_VA(dr, hh_B1, hh_beta1) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: hh_B1
    real(DP), intent(in)  :: hh_beta1
    real(DP)              :: val

    ! ---

    real(DP)  :: exp1

    ! ---

    exp1 = hh_B1*exp(-hh_beta1*dr)

    val  = - exp1

  endfunction hh_VA


  function hh_VR(dr, hh_A, hh_Q, hh_alpha) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: hh_A
    real(DP), intent(in)  :: hh_Q
    real(DP), intent(in)  :: hh_alpha
    real(DP)              :: val

    ! ---

    real(DP)  :: exp1, hlp1

    ! ---

    exp1 = hh_A*exp(-hh_alpha*dr)
    hlp1 = 1+hh_Q/dr

    val  = hlp1*exp1

  endfunction hh_VR


#ifdef EXP_CUT

  function cutoff_f(dr, l, h, m) result(val)
    implicit none
    
    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: l
    real(DP), intent(in)  :: h
    real(DP), intent(in)  :: m
    real(DP)              :: val

    ! ---

!    if (dr < m) then 
!       val = 2**(-(2*(dr-l)/(h-l))**cutoff_k)
!    else
!       val = 1-2**(-(2*(h-dr)/(h-l))**cutoff_k)
!    endif

    val = exp(-(2*(dr-l)/(h-l))**3)

  endfunction cutoff_f

#else

  function cutoff_f(dr, l, h, m) result(val)
    implicit none

    real(DP), intent(in)  :: dr
    real(DP), intent(in)  :: l
    real(DP), intent(in)  :: h
    real(DP), intent(in)  :: m
    real(DP)              :: val

    ! ---

    real(DP)  :: fca

    ! ---

    fca = pi / ( h - l )

    val   = 0.5 *  ( 1.0 + cos( fca*( dr-l ) ) )

  endfunction cutoff_f

#endif


  !**********************************************************************
  ! Make splines for attractive, repulsive functions
  !********************************************************************** 
  subroutine rebo2_db_make_splines(this)
    implicit none

    type(BOP_TYPE), intent(inout)  :: this

    ! ---

    real(DP)  :: cc_r2

    ! ---

    !
    ! Attractive potential
    !

#ifdef SCREENING
    cc_r2 = this%cc_ar_r2
#else
    cc_r2 = this%cc_in_r2
#endif

    call init( &
         this%spl_VA(C_C), &
         this%spl_n, this%spl_x0, cc_r2, &
         cc_VA, this%cc_B1, this%cc_B2, this%cc_B3, this%cc_beta1, this%cc_beta2, this%cc_beta3)
    call init( &
         this%spl_VA(C_H), &
         this%spl_n, this%spl_x0, this%ch_r2, &
         ch_VA, this%ch_B1, this%ch_beta1)
    call init( &
         this%spl_VA(H_H), &
         this%spl_n, this%spl_x0, this%hh_r2, &
         hh_VA, this%hh_B1, this%hh_beta1)

    !
    ! Repulsive potential
    !

    call init( &
         this%spl_VR(C_C), &
         this%spl_n, this%spl_x0, cc_r2, &
         cc_VR, this%cc_A, this%cc_Q, this%cc_alpha)
    call init( &
         this%spl_VR(C_H), &
         this%spl_n, this%spl_x0, this%ch_r2, &
         ch_VR, this%ch_A, this%ch_Q, this%ch_alpha)
    call init( &
         this%spl_VR(H_H), &
         this%spl_n, this%spl_x0, this%hh_r2, &
         hh_VR, this%hh_A, this%hh_Q, this%hh_alpha)

    !
    ! Inner cut-off
    !

    call init( &
         this%spl_fCin(C_C), &
         this%spl_n, this%cc_in_r1, this%cc_in_r2, &
         cutoff_f, this%cut_in_l(C_C), this%cut_in_h(C_C), this%cut_in_m(C_C))
    call init( &
         this%spl_fCin(C_H), &
         this%spl_n, this%ch_r1, this%ch_r2, &
         cutoff_f, this%cut_in_l(C_H), this%cut_in_h(C_H), this%cut_in_m(C_H))
    call init( &
         this%spl_fCin(H_H), &
         this%spl_n, this%hh_r1, this%hh_r2, &
         cutoff_f, this%cut_in_l(H_H), this%cut_in_h(H_H), this%cut_in_m(H_H))

#ifdef SCREENING

    !
    ! Attractive-repulsive cut-off
    !

    call init( &
         this%spl_fCar(C_C), &
         this%spl_n, this%cc_ar_r1, this%cc_ar_r2, &
         cutoff_f, this%cut_ar_l(C_C), this%cut_ar_h(C_C), this%cut_ar_m(C_C))
    call init( &
         this%spl_fCar(C_H), &
         this%spl_n, this%ch_r1, this%ch_r2, &
         cutoff_f, this%cut_ar_l(C_H), this%cut_ar_h(C_H), this%cut_ar_m(C_H))
    call init( &
         this%spl_fCar(H_H), &
         this%spl_n, this%hh_r1, this%hh_r2, &
         cutoff_f, this%cut_ar_l(H_H), this%cut_ar_h(H_H), this%cut_ar_m(H_H))

    !
    ! Bond-order cut-off
    !

    call init( &
         this%spl_fCbo(C_C), &
         this%spl_n, this%cc_bo_r1, this%cc_bo_r2, &
         cutoff_f, this%cut_bo_l(C_C), this%cut_bo_h(C_C), this%cut_bo_m(C_C))
    call init( &
         this%spl_fCbo(C_H), &
         this%spl_n, this%ch_r1, this%ch_r2, &
         cutoff_f, this%cut_bo_l(C_H), this%cut_bo_h(C_H), this%cut_bo_m(C_H))
    call init( &
         this%spl_fCbo(H_H), &
         this%spl_n, this%hh_r1, this%hh_r2, &
         cutoff_f, this%cut_bo_l(H_H), this%cut_bo_h(H_H), this%cut_bo_m(H_H))

    !
    ! Neighbor and conjugation cut-off
    !

    call init( &
         this%spl_fCnc(C_C), &
         this%spl_n, this%cc_nc_r1, this%cc_nc_r2, &
         cutoff_f, this%cut_nc_l(C_C), this%cut_nc_h(C_C), this%cut_nc_m(C_C))
    call init( &
         this%spl_fCnc(C_H), &
         this%spl_n, this%ch_r1, this%ch_r2, &
         cutoff_f, this%cut_nc_l(C_H), this%cut_nc_h(C_H), this%cut_nc_m(C_H))
    call init( &
         this%spl_fCnc(H_H), &
         this%spl_n, this%hh_r1, this%hh_r2, &
         cutoff_f, this%cut_nc_l(H_H), this%cut_nc_h(H_H), this%cut_nc_m(H_H))

#endif

  endsubroutine rebo2_db_make_splines


