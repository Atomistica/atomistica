!! ======================================================================
!! MDCORE - Interatomic potential library
!! https://github.com/pastewka/mdcore
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================

  public :: TERSOFF_MAX_REF, TERSOFF_MAX_EL, TERSOFF_MAX_PAIRS

  integer, parameter  :: TERSOFF_MAX_REF    = 1000

  integer, parameter  :: TERSOFF_MAX_EL     = 2
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

     real(DP)      :: A(TERSOFF_MAX_PAIRS)
     real(DP)      :: B(TERSOFF_MAX_PAIRS)
     real(DP)      :: xi(TERSOFF_MAX_PAIRS)
     real(DP)      :: lambda(TERSOFF_MAX_PAIRS)
     real(DP)      :: mu(TERSOFF_MAX_PAIRS)
     real(DP)      :: mubo(TERSOFF_MAX_PAIRS)
     integer       :: m(TERSOFF_MAX_PAIRS)
     real(DP)      :: beta(TERSOFF_MAX_EL)
     real(DP)      :: n(TERSOFF_MAX_EL)
     real(DP)      :: c(TERSOFF_MAX_EL)
     real(DP)      :: d(TERSOFF_MAX_EL)
     real(DP)      :: h(TERSOFF_MAX_EL)

     real(DP)      :: r1(TERSOFF_MAX_PAIRS)
     real(DP)      :: r2(TERSOFF_MAX_PAIRS)

#ifdef SCREENING
     real(DP)      :: or1(TERSOFF_MAX_PAIRS)       !< Outer cut-off start
     real(DP)      :: or2(TERSOFF_MAX_PAIRS)       !< Outer cut-off end

     real(DP)      :: bor1(TERSOFF_MAX_PAIRS)      !< Bond-order cut-off start
     real(DP)      :: bor2(TERSOFF_MAX_PAIRS)      !< Bond-order cut-off end

     real(DP)      :: Cmin(TERSOFF_MAX_PAIRS)      !< Inner screening parameter
     real(DP)      :: Cmax(TERSOFF_MAX_PAIRS)      !< Outer screening parameter
#endif

  endtype BOP_DB_TYPE


  type(BOP_DB_TYPE), parameter  :: Tersoff_PRB_39_5566_SiC = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Tersoff J., Phys. Rev. B 39, 5566 (1989)", &               ! ref
       reshape( (/  "C"," ",  "S","i", " "," " /), &
                (/ 2, TERSOFF_MAX_EL /) ), &                       ! el
       (/ 1.3936d3,    sqrt(1.3936d3*1.8308d3),  1.8308d3  /), &   ! A
       (/ 3.4674d2,    sqrt(3.4674d2*4.7118d2),  4.7118d2  /), &   ! B
       (/ 1.0_DP,      0.9776d0,                 1.0_DP    /), &   ! xi
       (/ 3.4879d0,    (3.4879d0+2.4799d0)/2,    2.4799d0  /), &   ! lambda
       (/ 2.2119d0,    (2.2119d0+1.7322d0)/2,    1.7322d0  /), &   ! mu
#ifdef SCREENING
       (/ 0.69103023078057590_DP, 0.56580821386164815_DP, 0.43569872294774004_DP /), & ! mubo
       (/  3,          3,                        3         /), &   ! m
#else
       (/ 0.0d0,       0.0d0,                    0.0d0     /), &   ! mubo
       (/  1,          1,                        1         /), &   ! m
#endif
       (/ 1.5724d-7,   1.1000d-6 /), &                             ! beta
       (/ 7.2751d-1,   7.8734d-1 /), &                             ! n
       (/ 3.8049d4,    1.0039d5 /), &                              ! c
       (/ 4.3484d0,    1.6217d1 /), &                              ! d
       (/ -5.7058d-1,  -5.9825d-1 /), &                            ! h
#ifdef SCREENING
       (/ 2.00_DP,        sqrt(2.00_DP*2.50_DP),        2.50_DP        /), &     ! r1
       (/ 2.00_DP*1.2_DP, sqrt(2.00_DP*2.50_DP)*1.2_DP, 2.50_DP*1.2_DP /), &     ! r2
       (/ 2.00_DP,        sqrt(2.00_DP*3.00_DP),        3.00_DP        /), &     ! or1
       (/ 2.00_DP*2.0_DP, sqrt(2.00_DP*3.00_DP)*2.0_DP, 3.00_DP*2.0_DP /), &     ! or2
       (/ 2.00_DP,        sqrt(2.00_DP*3.00_DP),        3.00_DP        /), &     ! bor1
       (/ 2.00_DP*2.0_DP, sqrt(2.00_DP*3.00_DP)*2.0_DP, 3.00_DP*2.0_DP /), &     ! bor2
       (/  1.0_DP,        1.0_DP,                   1.0_DP  /), &     ! Cmin
       (/  3.0_DP,        3.0_DP,                   3.0_DP  /) &      ! Cmax
#else
       (/ 1.80_DP,     sqrt(1.80_DP*2.70_DP),    2.70_DP /), &     ! r1
       (/ 2.10_DP,     sqrt(2.10_DP*3.00_DP),    3.00_DP /) &      ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter, private  :: tersoff_db(1) = (/ &
       Tersoff_PRB_39_5566_SiC &
       /)
