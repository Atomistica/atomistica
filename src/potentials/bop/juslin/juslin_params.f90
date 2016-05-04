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

  public :: JUSLIN_MAX_REF, JUSLIN_MAX_EL, JUSLIN_MAX_PAIRS

  integer, parameter  :: JUSLIN_MAX_REF    = 1000

  integer, parameter  :: JUSLIN_MAX_EL     = 3
  integer, parameter  :: JUSLIN_MAX_PAIRS  = PAIR_INDEX_NS(JUSLIN_MAX_EL, JUSLIN_MAX_EL, JUSLIN_MAX_EL)


  !>
  !! This data-type contains the parameter set.
  !<
  public :: BOP_DB_TYPE
  type BOP_DB_TYPE

     integer   :: nel = -1
     integer   :: nD0, nr0, nS, nbeta, ngamma, nc, nd, nh, nn, nalpha, nomega, nm, nr1, nr2

#ifdef SCREENING
     integer   :: nor1, nor2, nbor1, nbor2, nCmin, nCmax
#endif

     character(JUSLIN_MAX_REF)  :: ref = "*"

     character  :: el(2, JUSLIN_MAX_EL)

     real(DP)  :: D0(JUSLIN_MAX_EL**2)        !< Binding energy of the dimer
     real(DP)  :: r0(JUSLIN_MAX_EL**2)        !< Dimer bond distance
     real(DP)  :: S(JUSLIN_MAX_EL**2)         !< Slope of Pauling plot
     real(DP)  :: beta(JUSLIN_MAX_EL**2)      !< Dimer stiffness, i.e. vibrational frequency
     real(DP)  :: gamma(JUSLIN_MAX_EL**2)     !< Scaling factor for the bond-order
     real(DP)  :: c(JUSLIN_MAX_EL**2)         !< Angular parameters
     real(DP)  :: d(JUSLIN_MAX_EL**2)         !< Angular parameters
     real(DP)  :: h(JUSLIN_MAX_EL**2)         !< Angular parameters
     real(DP)  :: n(JUSLIN_MAX_EL**2)         !< Bond-order given by ( 1 + zij ** n )
     real(DP)  :: alpha(JUSLIN_MAX_EL**3)     !< Exponential bond-order contribution
     real(DP)  :: omega(JUSLIN_MAX_EL**3)     !< Scaling of exponential bond-order contribution
     integer   :: m(JUSLIN_MAX_EL**3)     !< Scaling of exponential bond-order contribution
     real(DP)  :: r1(JUSLIN_MAX_EL**2)        !< Cut-off start (inner cut-off if screening is enabled)
     real(DP)  :: r2(JUSLIN_MAX_EL**2)        !< Cut-off end (inner cut-off if screening is enabled)
#ifdef SCREENING
     real(DP)  :: or1(JUSLIN_MAX_EL**2)       !< Outer cut-off start
     real(DP)  :: or2(JUSLIN_MAX_EL**2)       !< Outer cut-off end

     real(DP)  :: bor1(JUSLIN_MAX_EL**2)      !< Bond-order cut-off start
     real(DP)  :: bor2(JUSLIN_MAX_EL**2)      !< Bond-order cut-off end

     real(DP)  :: Cmin(JUSLIN_MAX_EL**2)      !< Inner screening parameter
     real(DP)  :: Cmax(JUSLIN_MAX_EL**2)      !< Outer screening parameter
#endif

  endtype BOP_DB_TYPE


  type(BOP_DB_TYPE), parameter  :: Juslin_J_Appl_Phys_98_123520_WCH = BOP_DB_TYPE( &
       3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 27, 27, 27, 9, 9, &
#ifdef SCREENING
       9, 9, 9, 9, 9, 9, &
#endif
       "Juslin N. et al., J. Appl. Phys. 98, 123520 (2005)", &  ! ref
       reshape( (/  "W"," ",  "C"," ", "H"," " /), &
                (/ 2, JUSLIN_MAX_EL /) ), &                   ! el
       (/ 5.41861_DP,     6.64_DP,      2.748_DP,    0.0_DP,  6.0_DP,         3.6422_DP,      0.0_DP,  3.642_DP,    4.7509_DP  /), &  ! D0
       (/ 2.34095_DP,     1.90547_DP,   1.727_DP,   -1.0_DP,  1.39_DP,        1.1199_DP,     -1.0_DP,  1.1199_DP,   0.74144_DP /), &  ! r0
       (/ 1.92708_DP,     2.96149_DP,   1.2489_DP,   0.0_DP,  1.22_DP,        1.69077_DP,     0.0_DP,  1.69077_DP,  2.3432_DP  /), &  ! S
       (/ 1.38528_DP,     1.80370_DP,   1.52328_DP,  0.0_DP,  2.1_DP,         1.9583_DP,      0.0_DP,  1.9583_DP,   1.9436_DP  /), &  ! beta
       (/ 0.00188227_DP,  0.072855_DP,  0.0054_DP,   0.0_DP,  0.00020813_DP,  0.00020813_DP,  0.0_DP,  12.33_DP,    12.33_DP   /), &  ! gamma
       (/ 2.14969_DP,     1.10304_DP,   1.788_DP,    0.0_DP,  330.0_DP,       330.0_DP,       0.0_DP,  0.0_DP,      0.0_DP     /), &  ! c
       (/ 0.17126_DP,     0.33018_DP,   0.8255_DP,   0.0_DP,  3.5_DP,         3.5_DP,         0.0_DP,  1.0_DP,      1.0_DP     /), &  ! d
       (/-0.27780_DP,     0.75107_DP,   0.38912_DP,  0.0_DP,  1.0_DP,         1.0_DP,         0.0_DP,  1.0_DP,      1.0_DP     /), &  ! h
       (/ 1.0_DP,         1.0_DP,       1.0_DP,      0.0_DP,  1.0_DP,         1.0_DP,         0.0_DP,  1.0_DP,      1.0_DP     /), &  ! n
       (/ 0.45876_DP, 0.0_DP, 0.0_DP, 0.45876_DP, 0.0_DP, 0.0_DP, 0.45876_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 4.0_DP, 0.0_DP, 4.0_DP, 4.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 4.0_DP, 4.0_DP, 0.0_DP, 4.0_DP, 4.0_DP /), &  ! alpha
       (/ 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 2.94586_DP, 4.54415_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 0.33946_DP, 0.22006_DP, 1.0_DP, 1.0_DP, 1.0_DP /), &  ! omega
       (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), &  ! m
#ifdef SCREENING
       (/ 3.20_DP,        2.60_DP,      2.68_DP,     0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! r1
       (/ 3.80_DP,        3.00_DP,      2.96_DP,     0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /), &  ! r2
       (/ 3.20_DP,        2.60_DP,      2.68_DP,     0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! or1
       (/ 3.80_DP,        3.00_DP,      2.96_DP,     0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /), &  ! or2
       (/ 3.20_DP,        2.60_DP,      2.68_DP,     0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! bor1
       (/ 3.80_DP,        3.00_DP,      2.96_DP,     0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /), &  ! bor2
       (/ 1.00_DP,        1.00_DP,      1.00_DP,     0.0_DP,  1.00_DP,        1.00_DP,        0.0_DP,  1.00_DP,     1.00_DP    /), &  ! Cmin
       (/ 3.00_DP,        3.00_DP,      3.00_DP,     0.0_DP,  3.00_DP,        3.00_DP,        0.0_DP,  3.00_DP,     3.00_DP    /) &   ! Cmax
#else
       (/ 3.20_DP,        2.60_DP,      1.95_DP,     0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! r1
       (/ 3.80_DP,        3.00_DP,      2.35_DP,     0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /) &   ! r2
#endif
       )


  type(BOP_DB_TYPE), parameter  :: Kuopanportti_P_Comp_Mat_Sci_111_525_FeCH = BOP_DB_TYPE( &
       3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 27, 27, 27, 9, 9, &
#ifdef SCREENING
       9, 9, 9, 9, 9, 9, &
#endif
       "Kuopanportti P. et al., Comp. Mat. Sci. 111, 525 (2016)", &  ! ref
       reshape( (/  "F","e",  "C"," ", "H"," " /), &
                (/ 2, JUSLIN_MAX_EL /) ), &                   ! el
       (/ 1.5_DP,         4.82645134_DP,  1.630_DP,    0.0_DP,  6.0_DP,         3.6422_DP,      0.0_DP,  3.642_DP,    4.7509_DP  /), &  ! D0
       (/ 2.29_DP,        1.47736510_DP,  1.589_DP,   -1.0_DP,  1.39_DP,        1.1199_DP,     -1.0_DP,  1.1199_DP,   0.74144_DP /), &  ! r0
       (/ 2.0693_DP,      1.43134755_DP,  4.000_DP,    0.0_DP,  1.22_DP,        1.69077_DP,     0.0_DP,  1.69077_DP,  2.3432_DP  /), &  ! S
       (/ 1.4_DP,         1.63208170_DP,  1.875_DP,    0.0_DP,  2.1_DP,         1.9583_DP,      0.0_DP,  1.9583_DP,   1.9436_DP  /), &  ! beta
       (/ 0.01158_DP,     0.00205862_DP,  0.01332_DP,  0.0_DP,  0.00020813_DP,  0.00020813_DP,  0.0_DP,  12.33_DP,    12.33_DP   /), &  ! gamma
       (/ 1.2899_DP,      8.95583221_DP,  424.5_DP,    0.0_DP,  330.0_DP,       330.0_DP,       0.0_DP,  0.0_DP,      0.0_DP     /), &  ! c
       (/ 0.3413_DP,      0.72062047_DP,  7.282_DP,    0.0_DP,  3.5_DP,         3.5_DP,         0.0_DP,  1.0_DP,      1.0_DP     /), &  ! d
       (/-0.26_DP,        0.87099874_DP, -0.1091_DP,   0.0_DP,  1.0_DP,         1.0_DP,         0.0_DP,  1.0_DP,      1.0_DP     /), &  ! h
       (/ 1.0_DP,         1.0_DP,         1.0_DP,      0.0_DP,  1.0_DP,         1.0_DP,         0.0_DP,  1.0_DP,      1.0_DP     /), &  ! n
       (/ 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 4.0_DP, 0.0_DP, 4.0_DP, 4.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 4.0_DP, 4.0_DP, 0.0_DP, 4.0_DP, 4.0_DP /), &  ! alpha
       (/ 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 2.94586_DP, 4.54415_DP, 1.0_DP, 1.0_DP, 1.0_DP, 1.0_DP, 0.33946_DP, 0.22006_DP, 1.0_DP, 1.0_DP, 1.0_DP /), &  ! omega
       (/ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 /), &  ! m
#ifdef SCREENING
       (/ 2.95_DP,        2.30_DP,        2.2974_DP,   0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! r1
       (/ 3.35_DP,        2.70_DP,        2.6966_DP,   0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /), &  ! r2
       (/ 100.0_DP,       100.0_DP,       100.0_DP,    0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! or1
       (/ 3.35_DP,        2.70_DP,        2.6966_DP,   0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /), &  ! or2
       (/ 100.0_DP,       100.0_DP,       100.0_DP,    0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! bor1
       (/ 3.35_DP,        2.70_DP,        2.6966_DP,   0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /), &  ! bor2
       (/ 1.00_DP,        1.00_DP,        1.00_DP,     0.0_DP,  1.00_DP,        1.00_DP,        0.0_DP,  1.00_DP,     1.00_DP    /), &  ! Cmin
       (/ 3.00_DP,        3.00_DP,        3.00_DP,     0.0_DP,  3.00_DP,        3.00_DP,        0.0_DP,  3.00_DP,     3.00_DP    /) &   ! Cmax
#else
       (/ 2.95_DP,        2.30_DP,        2.2974_DP,   0.0_DP,  1.70_DP,        1.30_DP,        0.0_DP,  1.30_DP,     1.10_DP    /), &  ! r1
       (/ 3.35_DP,        2.70_DP,        2.6966_DP,   0.0_DP,  2.00_DP,        1.80_DP,        0.0_DP,  1.80_DP,     1.70_DP    /) &   ! r2
#endif
       )


  type(BOP_DB_TYPE), parameter, private  :: BOP_DB(2) = (/ &
       Juslin_J_Appl_Phys_98_123520_WCH, &
       Kuopanportti_P_Comp_Mat_Sci_111_525_FeCH &
       /)

