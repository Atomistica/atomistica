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

  public :: KUMAGAI_MAX_REF, KUMAGAI_MAX_EL, KUMAGAI_MAX_PAIRS

  integer, parameter  :: KUMAGAI_MAX_REF    = 1000

  integer, parameter  :: KUMAGAI_MAX_EL     = 1
  integer, parameter  :: KUMAGAI_MAX_PAIRS  = &
       PAIR_INDEX(KUMAGAI_MAX_EL, KUMAGAI_MAX_EL, KUMAGAI_MAX_EL)

  !>
  !! This data-type contains the parameter set.
  !<
  public :: BOP_DB_TYPE
  type BOP_DB_TYPE

     integer  :: nel = -1
     integer  :: nA
     integer  :: nB
     integer  :: nlambda1
     integer  :: nlambda2
     integer  :: neta
     integer  :: ndelta
     integer  :: nalpha
     integer  :: nbeta
     integer  :: nc1
     integer  :: nc2
     integer  :: nc3
     integer  :: nc4
     integer  :: nc5
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

     character(KUMAGAI_MAX_REF)  :: ref

     character     :: el(2, KUMAGAI_MAX_EL)

     real(DP)      :: A(KUMAGAI_MAX_PAIRS)
     real(DP)      :: B(KUMAGAI_MAX_PAIRS)
     real(DP)      :: lambda1(KUMAGAI_MAX_PAIRS)
     real(DP)      :: lambda2(KUMAGAI_MAX_PAIRS)
     real(DP)      :: eta(KUMAGAI_MAX_EL)
     real(DP)      :: delta(KUMAGAI_MAX_EL)
     real(DP)      :: alpha(KUMAGAI_MAX_PAIRS)
     integer       :: beta(KUMAGAI_MAX_PAIRS)
     real(DP)      :: c1(KUMAGAI_MAX_EL)
     real(DP)      :: c2(KUMAGAI_MAX_EL)
     real(DP)      :: c3(KUMAGAI_MAX_EL)
     real(DP)      :: c4(KUMAGAI_MAX_EL)
     real(DP)      :: c5(KUMAGAI_MAX_EL)
     real(DP)      :: h(KUMAGAI_MAX_EL)

     real(DP)      :: r1(KUMAGAI_MAX_PAIRS)
     real(DP)      :: r2(KUMAGAI_MAX_PAIRS)

#ifdef SCREENING
     real(DP)      :: or1(KUMAGAI_MAX_PAIRS)       !< Outer cut-off start
     real(DP)      :: or2(KUMAGAI_MAX_PAIRS)       !< Outer cut-off end

     real(DP)      :: bor1(KUMAGAI_MAX_PAIRS)      !< Bond-order cut-off start
     real(DP)      :: bor2(KUMAGAI_MAX_PAIRS)      !< Bond-order cut-off end

     real(DP)      :: Cmin(KUMAGAI_MAX_PAIRS)      !< Inner screening parameter
     real(DP)      :: Cmax(KUMAGAI_MAX_PAIRS)      !< Outer screening parameter
#endif
  endtype BOP_DB_TYPE


  type(BOP_DB_TYPE), parameter  :: Kumagai_CompMaterSci_39_457_Si = BOP_DB_TYPE( &
       1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
#ifdef SCREENING
       1, 1, 1, 1, 1, 1, &
#endif
       "Kumagai, Izumi, Hara, Sakai, Comp. Mater. Sci. 39, 457 (2007)", & ! ref
       reshape( (/  "S","i",  " "," ", " "," " /), &
                (/ 2, KUMAGAI_MAX_EL /) ), &                       ! el
       (/ 3281.5905_DP   /), &   ! A
       (/ 121.00047_DP   /), &   ! B
       (/ 3.2300135_DP   /), &   ! lambda1
       (/ 1.3457970_DP   /), &   ! lambda2
       (/ 1.0000000_DP   /), &   ! eta
       (/ 0.53298909_DP  /), &   ! delta
       (/ 2.3890327_DP   /), &   ! alpha
       (/ 1              /), &   ! beta
       (/ 0.20173476_DP  /), &   ! n
       (/ 730418.72_DP   /), &   ! c
       (/ 1000000.0_DP   /), &   ! d
       (/ 1.0000000_DP   /), &   ! h
       (/ 26.000000_DP   /), &   ! d
       (/ -0.36500000_DP /), &   ! h
#ifdef SCREENING
       (/ 2.50_DP        /), &   ! r1
       (/ 2.50_DP*1.2_DP /), &   ! r2
       (/ 3.00_DP        /), &   ! or1
       (/ 3.00_DP*2.0_DP /), &   ! or2
       (/ 3.00_DP        /), &   ! bor1
       (/ 3.00_DP*2.0_DP /), &   ! bor2
       (/ 1.0_DP         /), &   ! Cmin
       (/ 3.0_DP         /)  &   ! Cmax
#else
       (/ 2.70_DP        /), &   ! r1
       (/ 3.30_DP        /)  &   ! r2
#endif
       )

  type(BOP_DB_TYPE), parameter, private  :: kumagai_db(1) = (/ &
       Kumagai_CompMaterSci_39_457_Si &
       /)
