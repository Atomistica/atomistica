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

  public :: TERSOFF_MAX_REF, TERSOFF_MAX_EL, TERSOFF_MAX_PAIRS

  integer, parameter  :: TERSOFF_MAX_REF    = 1000

  integer, parameter  :: TERSOFF_MAX_EL     = 3
  integer, parameter  :: TERSOFF_MAX_PAIRS  = PAIR_INDEX(TERSOFF_MAX_EL, TERSOFF_MAX_EL, TERSOFF_MAX_EL)

  !>
  !! This data-type contains the parameter set.
  !<
  public :: BOP_DB_TYPE
  type BOP_DB_TYPE

     integer  :: nel = -1
     integer  :: nA
     integer  :: nB
     integer  :: nxi
     integer  :: nlambda
     integer  :: nmu
     integer  :: nomega
     integer  :: nmubo
     integer  :: nm
     integer  :: nbeta
     integer  :: nn
     integer  :: nc
     integer  :: nd
     integer  :: nh
     integer  :: nr1
     integer  :: nr2
#ifdef SCREENING
     integer  :: nor1
     integer  :: nor2
     integer  :: nbor1
     integer  :: nbor2
     integer  :: nCmin
     integer  :: nCmax
#endif

     character(TERSOFF_MAX_REF)  :: ref

     character     :: el(2, TERSOFF_MAX_EL)

     real(DP)      :: A(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: B(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: xi(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: lambda(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: mu(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: omega(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: mubo(TERSOFF_MAX_PAIRS) = 0.0_DP
     integer       :: m(TERSOFF_MAX_PAIRS) = 1
     real(DP)      :: beta(TERSOFF_MAX_EL) = 1.0_DP
     real(DP)      :: n(TERSOFF_MAX_EL) = 1.0_DP
     real(DP)      :: c(TERSOFF_MAX_EL) = 1.0_DP
     real(DP)      :: d(TERSOFF_MAX_EL) = 1.0_DP
     real(DP)      :: h(TERSOFF_MAX_EL) = 1.0_DP

     real(DP)      :: r1(TERSOFF_MAX_PAIRS) = 1.0_DP
     real(DP)      :: r2(TERSOFF_MAX_PAIRS) = 2.0_DP

#ifdef SCREENING
     real(DP)      :: or1(TERSOFF_MAX_PAIRS)       !< Outer cut-off start
     real(DP)      :: or2(TERSOFF_MAX_PAIRS)       !< Outer cut-off end

     real(DP)      :: bor1(TERSOFF_MAX_PAIRS)      !< Bond-order cut-off start
     real(DP)      :: bor2(TERSOFF_MAX_PAIRS)      !< Bond-order cut-off end

     real(DP)      :: Cmin(TERSOFF_MAX_PAIRS)      !< Inner screening parameter
     real(DP)      :: Cmax(TERSOFF_MAX_PAIRS)      !< Outer screening parameter
#endif

  endtype BOP_DB_TYPE


#define FILL1 0.0_DP
#define FILL3 0.0_DP,0.0_DP,0.0_DP
#define FILL3i 0,0,0

  type(BOP_DB_TYPE), parameter  :: Tersoff_PRB_39_5566_SiC = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Tersoff J., Phys. Rev. B 39, 5566 (1989)", &               ! ref
       reshape( (/  "C"," ",  "S","i", " "," " /), &
                (/ 2, TERSOFF_MAX_EL /) ), &                       ! el
       (/ 1.3936d3,    sqrt(1.3936d3*1.8308d3),  1.8308d3,  FILL3  /), &   ! A
       (/ 3.4674d2,    sqrt(3.4674d2*4.7118d2),  4.7118d2,  FILL3  /), &   ! B
       (/ 1.0_DP,      0.9776d0,                 1.0_DP,    FILL3  /), &   ! xi
       (/ 3.4879d0,    (3.4879d0+2.4799d0)/2,    2.4799d0,  FILL3  /), &   ! lambda
       (/ 2.2119d0,    (2.2119d0+1.7322d0)/2,    1.7322d0,  FILL3  /), &   ! mu
       (/ 1d0,         1d0,                      1d0,       FILL3  /), &   ! omega
#ifdef SCREENING
       (/ 0.69103023078057590_DP, 0.56580821386164815_DP, 0.43569872294774004_DP, FILL3 /), & ! mubo
       (/  3,          3,                        3,         FILL3i /), &   ! m
#else
       (/ 0.0d0,       0.0d0,                    0.0d0,     FILL3  /), &   ! mubo
       (/  1,          1,                        1,         FILL3i /), &   ! m
#endif
       (/ 1.5724d-7,   1.1000d-6,  FILL1 /), &                             ! beta
       (/ 7.2751d-1,   7.8734d-1,  FILL1 /), &                             ! n
       (/ 3.8049d4,    1.0039d5,   FILL1 /), &                              ! c
       (/ 4.3484d0,    1.6217d1,   FILL1 /), &                              ! d
       (/ -5.7058d-1,  -5.9825d-1, FILL1 /), &                            ! h
#ifdef SCREENING
       (/ 2.00_DP,        sqrt(2.00_DP*2.50_DP),        2.50_DP,         FILL3 /), &     ! r1
       (/ 2.00_DP*1.2_DP, sqrt(2.00_DP*2.50_DP)*1.2_DP, 2.50_DP*1.2_DP,  FILL3 /), &     ! r2
       (/ 2.00_DP,        sqrt(2.00_DP*3.00_DP),        3.00_DP       ,  FILL3 /), &     ! or1
       (/ 2.00_DP*2.0_DP, sqrt(2.00_DP*3.00_DP)*2.0_DP, 3.00_DP*2.0_DP,  FILL3 /), &     ! or2
       (/ 2.00_DP,        sqrt(2.00_DP*3.00_DP),        3.00_DP       ,  FILL3 /), &     ! bor1
       (/ 2.00_DP*2.0_DP, sqrt(2.00_DP*3.00_DP)*2.0_DP, 3.00_DP*2.0_DP,  FILL3 /), &     ! bor2
       (/  1.0_DP,        1.0_DP,                       1.0_DP,          FILL3 /), &     ! Cmin
       (/  3.0_DP,        3.0_DP,                       3.0_DP,          FILL3 /) &      ! Cmax
#else
       (/ 1.80_DP,        sqrt(1.80_DP*2.70_DP),        2.70_DP,         FILL3 /), &     ! r1
       (/ 2.10_DP,        sqrt(2.10_DP*3.00_DP),        3.00_DP,         FILL3 /) &      ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter, private  :: tersoff_db(1) = (/ &
       Tersoff_PRB_39_5566_SiC &
       /)
