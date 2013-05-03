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

  public :: BRENNER_MAX_REF, BRENNER_MAX_EL, BRENNER_MAX_PAIRS

  integer, parameter  :: BRENNER_MAX_REF    = 1000

  integer, parameter  :: BRENNER_MAX_EL     = 2
  integer, parameter  :: BRENNER_MAX_PAIRS  = PAIR_INDEX(BRENNER_MAX_EL, BRENNER_MAX_EL, BRENNER_MAX_EL)


  !>
  !! This data-type contains the parameter set.
  !<
  public :: BOP_DB_TYPE
  type BOP_DB_TYPE

     integer   :: nel = -1
     integer   :: nD0, nr0, nS, nbeta, ngamma, nc, nd, nh, nmu, nn, nm, nr1, nr2

#ifdef SCREENING
     integer   :: nor1, nor2, nbor1, nbor2, nCmin, nCmax
#endif

     character(BRENNER_MAX_REF)  :: ref = "*"

     character  :: el(2, BRENNER_MAX_EL)

     real(DP)  :: D0(BRENNER_MAX_PAIRS)        !< Binding energy of the dimer
     real(DP)  :: r0(BRENNER_MAX_PAIRS)        !< Dimer bond distance
     real(DP)  :: S(BRENNER_MAX_PAIRS)         !< Slope of Pauling plot
     real(DP)  :: beta(BRENNER_MAX_PAIRS)      !< Dimer stiffness, i.e. vibrational frequency
     real(DP)  :: gamma(BRENNER_MAX_PAIRS)     !< Scaling factor for the bond-order
     real(DP)  :: c(BRENNER_MAX_PAIRS)         !< Angular parameters
     real(DP)  :: d(BRENNER_MAX_PAIRS)         !< Angular parameters
     real(DP)  :: h(BRENNER_MAX_PAIRS)         !< Angular parameters
     real(DP)  :: mu(BRENNER_MAX_PAIRS)        !< Exponential bond-order contribution
     real(DP)  :: n(BRENNER_MAX_PAIRS)         !< Bond-order given by ( 1 + zij ** n )
     integer   :: m(BRENNER_MAX_PAIRS)         !< Distance dependent part of bo given by exp(2(mu*dr)**m)
     real(DP)  :: r1(BRENNER_MAX_PAIRS)        !< Cut-off start (inner cut-off if screening is enabled)
     real(DP)  :: r2(BRENNER_MAX_PAIRS)        !< Cut-off end (inner cut-off if screening is enabled)
#ifdef SCREENING
     real(DP)  :: or1(BRENNER_MAX_PAIRS)       !< Outer cut-off start
     real(DP)  :: or2(BRENNER_MAX_PAIRS)       !< Outer cut-off end

     real(DP)  :: bor1(BRENNER_MAX_PAIRS)      !< Bond-order cut-off start
     real(DP)  :: bor2(BRENNER_MAX_PAIRS)      !< Bond-order cut-off end

     real(DP)  :: Cmin(BRENNER_MAX_PAIRS)      !< Inner screening parameter
     real(DP)  :: Cmax(BRENNER_MAX_PAIRS)      !< Outer screening parameter
#endif

  endtype BOP_DB_TYPE


  type(BOP_DB_TYPE), parameter  :: Erhart_PRB_71_035211_SiC = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Erhart P., Albe K., Phys. Rev. B 71, 035211 (2005)", &     ! ref
       reshape( (/  "C"," ",  "S","i", " "," " /), &
                (/ 2, BRENNER_MAX_EL /)), &                       ! el
       (/  6.00_DP,      4.36_DP,       3.24_DP      /), &         ! D0
       (/  1.4276_DP,    1.79_DP,       2.232_DP     /), &         ! r0
       (/  2.167_DP,     1.847_DP,      1.842_DP     /), &         ! S
       (/  2.0099_DP,    1.6991_DP,     1.4761_DP    /), &         ! beta
       (/  0.11233_DP,   0.011877_DP,   0.114354_DP  /), &         ! gamma
       (/  181.910_DP,   273987.0_DP,   2.00494_DP   /), &         ! c
       (/  6.28433_DP,   180.314_DP,    0.81472_DP   /), &         ! d
       (/  0.5556_DP,    0.68_DP,       0.259_DP     /), &         ! h
#ifdef SCREENING
       (/  1.0_DP/1.4276_DP,  1.0_DP/1.79_DP,                1.0_DP/2.232_DP   /), &   ! mu
       (/  1.0_DP,            1.0_DP,                        1.0_DP            /), &   ! n
       (/  3,                 3,                             3                 /), &   ! m
       (/  2.00_DP,           sqrt(2.00_DP*2.50_DP),         2.50_DP           /), &   ! r1
       (/  2.00_DP*1.2_DP,    sqrt(2.00_DP*2.50_DP)*1.2_DP,  2.50_DP*1.2_DP    /), &   ! r2
       (/  2.00_DP,           sqrt(2.00_DP*3.00_DP),         3.00_DP           /), &   ! or1
       (/  2.00_DP*2.0_DP,    sqrt(2.00_DP*3.00_DP)*2.0_DP,  3.00_DP*2.0_DP    /), &   ! or2
       (/  2.00_DP,           sqrt(2.00_DP*3.00_DP),         3.00_DP           /), &   ! bor1
       (/  2.00_DP*2.0_DP,    sqrt(2.00_DP*3.00_DP)*2.0_DP,  3.00_DP*2.0_DP    /), &   ! bor2
       (/  1.0_DP,            1.0_DP,                        1.0_DP            /), &   ! Cmin
       (/  3.0_DP,            3.0_DP,                        3.0_DP            /) &    ! Cmax
#else
       (/  0.0_DP,       0.0_DP,        0.0_DP       /), &         ! mu
       (/  1.0_DP,       1.0_DP,        1.0_DP       /), &         ! n
       (/  1,            1,             1            /), &         ! m
       (/  1.85_DP,      2.20_DP,       2.68_DP      /), &         ! r1
       (/  2.15_DP,      2.60_DP,       2.96_DP      /) &          ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter  :: Albe_PRB_65_195124_PtC = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Albe K., Nordlund K., Averback R. S., Phys. Rev. B 65, 195124 (2002)", &     ! ref
       reshape( (/  "P","t",  "C"," ", " "," " /), &
                (/ 2, BRENNER_MAX_EL /)), &                       ! el
       (/  3.683_DP,      5.3_DP,     6.0_DP         /), &         ! D0
       (/  2.384_DP,      1.84_DP,    1.39_DP        /), &         ! r0
       (/  2.24297_DP,    1.1965_DP,  1.22_DP        /), &         ! S
       (/  1.64249_DP,    1.836_DP,   2.1_DP         /), &         ! beta
       (/  0.0008542_DP,  0.0097_DP,  0.00020813_DP  /), &         ! gamma
       (/  34.0_DP,       1.23_DP,    330.0_DP       /), &         ! c
       (/  1.1_DP,        0.36_DP,    3.5_DP         /), &         ! d
       (/  1.0_DP,        1.0_DP,     1.0_DP         /), &         ! h
       (/  1.335_DP,      0.0_DP,     0.0_DP         /), &         ! mu
       (/  1.0_DP,        1.0_DP,     1.0_DP         /), &         ! n
       (/  1,             1,          1              /), &         ! m
#ifdef SCREENING
       (/  2.9_DP,        2.5_DP,     1.7_DP         /), &         ! r1
       (/  3.3_DP,        2.8_DP,     2.0_DP         /), &         ! r2
       (/  2.9_DP,        2.5_DP,     1.7_DP         /), &         ! or1
       (/  3.3_DP,        2.8_DP,     2.0_DP         /), &         ! or2
       (/  2.9_DP,        2.5_DP,     1.7_DP         /), &         ! bor1
       (/  3.3_DP,        2.8_DP,     2.0_DP         /), &         ! bor2
       (/  1.0_DP,        1.0_DP,     1.0_DP         /), &         ! Cmin
       (/  3.0_DP,        3.0_DP,     3.0_DP         /) &          ! Cmax
#else
       (/  2.9_DP,        2.5_DP,     1.7_DP         /), &         ! r1
       (/  3.3_DP,        2.8_DP,     2.0_DP         /) &         ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter  :: Henriksson_PRB_79_144107_FeC = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Henriksson K.O.E., Nordlund K., Phys. Rev. B 79, 144107 (2009)", & ! ref
       reshape( &
       (/  "F","e",  "C"," ", " "," " /), &
       (/ 2, BRENNER_MAX_EL /)), &                       ! el
       (/  1.5_DP,         4.82645134_DP,   6.0_DP           /), &  ! D0
       (/  2.29_DP,        1.47736510_DP,   1.39_DP          /), &  ! r0
       (/  2.0693109_DP,   1.43134755_DP,   1.22_DP          /), &  ! S
       (/  1.4_DP,         1.63208170_DP,   2.1_DP           /), &  ! beta
       (/  0.0115751_DP,   0.00205862_DP,   0.00020813_DP    /), &  ! gamma
       (/  1.2898716_DP,   8.95583221_DP,   330.0_DP         /), &  ! c
       (/  0.3413219_DP,   0.72062047_DP,   3.5_DP           /), &  ! d 
       (/ -0.26_DP,        0.87099874_DP,   1.0_DP           /), &  ! h
#ifdef SCREENING
       (/  0.0_DP,         0.0_DP,          1.0_DP/1.315_DP  /), &  ! mu
       (/  1.0_DP,         1.0_DP,          1.0_DP           /), &  ! n
       (/  1,              1,               3                /), &  ! m
       (/  2.95_DP,        2.3_DP,          2.00_DP          /), &  ! r1
       (/  3.35_DP,        2.7_DP,          1.2_DP*2.00_DP   /), &  ! r2
       (/  100.0_DP,       100.0_DP,        2.00_DP          /), &  ! or1
       (/  3.35_DP,        2.7_DP,          2.0_DP*2.00_DP   /), &  ! or2
       (/  100.0_DP,       100.0_DP,        1.20_DP          /), &  ! bor1
       (/  3.35_DP,        2.7_DP,          2.0_DP*2.00_DP   /), &  ! bor2
       (/  1.00_DP,        1.00_DP,         1.00_DP          /), &  ! Cmin
       (/  3.00_DP,        3.00_DP,         3.00_DP          /) &   ! Cmax
#else
       (/  0.0_DP,         0.0_DP,          0.0_DP           /), &  ! mu
       (/  1.0_DP,         1.0_DP,          1.0_DP           /), &  ! n
       (/  1,              1,               1                /), &  ! m
       (/  2.95_DP,        2.3_DP,          1.70_DP          /), &  ! r1
       (/  3.35_DP,        2.7_DP,          2.00_DP          /) &   ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter  :: Kioseoglou_PSSb_245_1118_AlN = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Kioseoglou J., Komninou Ph., Karakostas Th., Phys. Stat. Sol. (b) 245, 1118 (2008)", &
       reshape( (/  "N"," ",  "A","l", " "," " /), &
                (/ 2, BRENNER_MAX_EL /) ), &                      ! el
       (/  9.9100_DP,     3.3407_DP,     1.5000_DP   /), &         ! D0
       (/  1.1100_DP,     1.8616_DP,     2.4660_DP   /), &         ! r0
       (/  1.4922_DP,     1.7269_DP,     2.7876_DP   /), &         ! S
       (/  2.05945_DP,    1.7219_DP,     1.0949_DP   /), &         ! beta
       (/  0.76612_DP,    1.1e-6_DP,     0.3168_DP   /), &         ! gamma
       (/  0.178493_DP,   100390.0_DP,   0.0748_DP   /), &         ! c
       (/  0.20172_DP,    16.2170_DP,    19.5691_DP  /), &         ! d
       (/  0.045238_DP,   0.5980_DP,     0.6593_DP   /), &         ! h
       (/  0.0_DP,        0.0_DP,        0.0_DP      /), &         ! mu
       (/  1.0_DP,        0.7200_DP,     6.0865_DP   /), &         ! n
       (/  1,             1,             1           /), &         ! m
#ifdef SCREENING
       (/  2.00_DP,       2.19_DP,       2.60_DP     /), &         ! r1
       (/  2.40_DP,       2.49_DP,       2.80_DP     /), &         ! r2
       (/  2.00_DP,       2.19_DP,       3.40_DP     /), &         ! or1
       (/  2.40_DP,       2.49_DP,       3.60_DP     /), &         ! or2
       (/  2.00_DP,       2.19_DP,       3.40_DP     /), &         ! bor1
       (/  2.40_DP,       2.49_DP,       3.60_DP     /), &         ! bor2
       (/  1.00_DP,       1.00_DP,       1.00_DP     /), &         ! Cmin
       (/  3.00_DP,       3.00_DP,       3.00_DP     /) &          ! Cmax
#else
       (/  2.00_DP,       2.19_DP,       3.40_DP     /), &         ! r1 
       (/  2.40_DP,       2.49_DP,       3.60_DP     /) &          ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter  :: Klemenz_unpublished_SiC = BOP_DB_TYPE( &
       2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
#ifdef SCREENING
       3, 3, 3, 3, 3, 3, &
#endif
       "Klemenz A. (2011)", &     ! ref
       reshape( (/  "C"," ",  "S","i", " "," " /), &
                (/ 2, BRENNER_MAX_EL /)), &                                        ! el
       (/  9.18655066_DP,      3.74012106_DP,      2.44699407_DP      /), &         ! D0
       (/  1.31756335_DP,      1.83432402_DP,      2.33202404_DP      /), &         ! r0
       (/  2.00832275_DP,      1.65091787_DP,      2.03923608_DP      /), &         ! S
       (/  2.00674463_DP,      1.70321937_DP,      1.46034655_DP      /), &         ! beta
       (/  4.62525610e-01_DP,  4.34885606e-03_DP,  1.60279303e-02_DP  /), &         ! gamma
       (/  1.00506717e+02_DP,  3.60235187e+05_DP,  1.76882660e+01_DP  /), &         ! c
       (/  6.16010053_DP,      175.38330641_DP,    1.57043362_DP      /), &         ! d
       (/  0.44131887_DP,      0.59895772_DP,      0.39716596_DP      /), &         ! h
       (/  2.00258220_DP,      1.87466146_DP,      1.44622936_DP      /), &         ! mu
       (/  1.0_DP,             1.0_DP,             1.0_DP             /), &         ! n
       (/  1,                  1,                  1                  /), &         ! m
#ifdef SCREENING
       (/  2.00_DP,           sqrt(2.00_DP*2.50_DP),         2.50_DP           /), &   ! r1
       (/  2.00_DP*1.2_DP,    sqrt(2.00_DP*2.50_DP)*1.2_DP,  2.50_DP*1.2_DP    /), &   ! r2
       (/  2.00_DP,           sqrt(2.00_DP*3.00_DP),         3.00_DP           /), &   ! or1
       (/  2.00_DP*2.0_DP,    sqrt(2.00_DP*3.00_DP)*2.0_DP,  3.00_DP*2.0_DP    /), &   ! or2
       (/  2.00_DP,           sqrt(2.00_DP*3.00_DP),         3.00_DP           /), &   ! bor1
       (/  2.00_DP*2.0_DP,    sqrt(2.00_DP*3.00_DP)*2.0_DP,  3.00_DP*2.0_DP    /), &   ! bor2
       (/  1.0_DP,            1.0_DP,                        1.0_DP            /), &   ! Cmin
       (/  3.0_DP,            3.0_DP,                        3.0_DP            /) &    ! Cmax
#else
       (/  1.85_DP,            2.20_DP,            2.68_DP            /), &         ! r1
       (/  2.15_DP,            2.60_DP,            2.96_DP            /) &          ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter, private  :: BOP_DB(5) = (/ &
       Erhart_PRB_71_035211_SiC, &
       Henriksson_PRB_79_144107_FeC, &
       Albe_PRB_65_195124_PtC, &
       Kioseoglou_PSSb_245_1118_AlN, &
       Klemenz_unpublished_SiC &
       /)

