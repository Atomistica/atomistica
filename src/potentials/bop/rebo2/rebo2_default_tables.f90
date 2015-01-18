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

module rebo2_default_tables
  use supplib

  implicit none

contains

  !**********************************************************************
  ! Construct the Fcc-table
  !********************************************************************** 
  subroutine rebo2_default_Fcc_table(F, dFdi, dFdj, dFdk)
    implicit none

    real(DP), intent(out)  :: F(0:4, 0:4, 0:9)
    real(DP), intent(out)  :: dFdi(0:4, 0:4, 0:9)
    real(DP), intent(out)  :: dFdj(0:4, 0:4, 0:9)
    real(DP), intent(out)  :: dFdk(0:4, 0:4, 0:9)

    ! ---

    real(DP)  :: x
    integer   :: i, j, k

    ! ---

    F(:, :, :)     = 0.0_DP
    dFdi(:, :, :)  = 0.0_DP
    dFdj(:, :, :)  = 0.0_DP
    dFdk(:, :, :)  = 0.0_DP


    !
    ! Values from Table 4
    !

    F(1, 1,   0)  =  0.105000_DP
    F(1, 1,   1)  = -0.0041775_DP
    F(1, 1, 2:8)  = -0.0160856_DP
    F(2, 2,   0)  =  0.09444957_DP
    F(2, 2,   1)  =  0.02200000_DP
    F(2, 2,   2)  =  0.03970587_DP
    F(2, 2,   3)  =  0.03308822_DP
    F(2, 2,   4)  =  0.02647058_DP
    F(2, 2,   5)  =  0.01985293_DP
    F(2, 2,   6)  =  0.01323529_DP
    F(2, 2,   7)  =  0.00661764_DP
    F(2, 2,   8)  =  0.0_DP
    F(0, 1,   0)  =  0.04338699_DP
    F(0, 1,   1)  =  0.0099172158_DP

!--
    F(0, 1, 1:8)  =  0.0099172158_DP
!--

    F(0, 2,   0)  =  0.0493976637_DP
    F(0, 2,   1)  = -0.011942669_DP

!--
    F(0, 2, 2:8)  = F(0, 1, 1)
!--

    F(0, 3, 0:8)  = -0.119798935_DP

!--
    F(0, 3, 0:1)  = -0.119798935_DP
    F(0, 3, 2:8)  =  F(0, 1, 1)
!--

    F(1, 2,   0)  =  0.0096495698_DP
    F(1, 2,   1)  =  0.030_DP
    F(1, 2,   2)  = -0.0200_DP
    F(1, 2,   3)  = -0.0233778774_DP
    F(1, 2,   4)  = -0.0267557548_DP
    F(1, 2, 5:8)  = -0.030133632_DP
!-- Refit for proper graphene elastic constants
!    F(1, 2, 5:8)  = -0.030133632_DP - 2*0.09968441349_DP
!--
    F(1, 3, 1:8)  = -0.124836752_DP
    F(2, 3, 0:8)  = -0.044709383_DP

!-- Refit for proper graphene elastic constants
!    F(2, 2,   8)  = g_graphene_F
!    F(2, 2,   8)  = 0.0_DP

    do i = 3, 7
       F(2, 2, i) = F(2, 2, 2) + (i-2)*( F(2, 2, 8) - F(2, 2, 2) )/6
    enddo

    do i = 3, 4
       F(1, 2, i) = F(1, 2, 2) + (i-2)*( F(1, 2, 5) - F(1, 2, 2) )/3
!       write (*, *) i, F(1, 2, i)
    enddo
!--

    dFdi(2, 1,   0)  = -0.052500_DP
    dFdi(2, 1, 4:8)  = -0.054376_DP
    dFdi(2, 3,   0)  =  0.0_DP
    dFdi(2, 3, 1:5)  =  0.062418_DP
    dFdk(2, 2, 3:7)  = -0.006618_DP
    dFdi(2, 3, 6:8)  =  0.062418_DP
    dFdk(1, 1,   1)  = -0.060543_DP
    dFdk(1, 2,   3)  = -0.020044_DP
    dFdk(1, 2,   4)  = -0.020044_DP


    !
    ! Symmetrize values
    !

    do k = 0, 9
       do i = 0, 3
          do j = i+1, 3
             x = F(i, j, k) + F(j, i, k)
             F(i, j, k) = x
             F(j, i, k) = x

             x = dFdi(i, j, k) + dFdj(j, i, k)
             dFdi(i, j, k) = x
             dFdj(j, i, k) = x

             x = dFdi(j, i, k) + dFdj(i, j, k)
             dFdi(j, i, k) = x
             dFdj(i, j, k) = x

             x = dFdk(i, j, k) + dFdk(j, i, k)
             dFdk(i, j, k) = x
             dFdk(j, i, k) = x
          enddo
       enddo
    enddo

#ifdef ZERO_TABLES
    F     = 0.0_DP
    dFdi  = 0.0_DP
    dFdj  = 0.0_DP
    dFdk  = 0.0_DP
#endif

  endsubroutine rebo2_default_Fcc_table


  !**********************************************************************
  ! Construct the Fch-table
  !********************************************************************** 
  subroutine rebo2_default_Fch_table(F)
    implicit none

    real(DP), intent(out)  :: F(0:4, 0:4, 0:9)

    ! ---

    real(DP)  :: x
    integer   :: i, j, k

    ! ---

    F(:, :, :)     = 0.0_DP


    !
    ! Values from Table 9
    !

    F(0, 2, 4:8)  = -0.0090477875161288110_DP
    F(1, 3, 0:8)  = -0.213_DP
    F(1, 2, 0:8)  = -0.25_DP
    F(1, 1, 0:8)  = -0.5_DP


    !
    ! Symmetrize values
    !

    do k = 0, 9
       do i = 0, 2
          do j = i+1, 3

             x = F(i, j, k) + F(j, i, k)
             F(i, j, k) = x
             F(j, i, k) = x

          enddo
       enddo
    enddo

#ifdef ZERO_TABLES
    F = 0.0_DP
#endif

  endsubroutine rebo2_default_Fch_table


  !**********************************************************************
  ! Construct the Fhh-table
  !********************************************************************** 
  subroutine rebo2_default_Fhh_table(F)
    implicit none

    real(DP), intent(out)  :: F(0:4, 0:4, 0:9)

    ! ---

    F(:, :, :)     = 0.0_DP


    !
    ! Values from Table 6
    !

    F(1, 1, 0)  = 0.249831916_DP

#ifdef ZERO_TABLES
    F = 0.0_DP
#endif

  endsubroutine rebo2_default_Fhh_table


  !**********************************************************************
  ! Construct the Pcc-table
  !********************************************************************** 
  subroutine rebo2_default_Pcc_table(P)
    implicit none

    real(DP), intent(out)  :: P(0:5, 0:5)

    ! ---

    P(:, :)     = 0.0_DP


    !
    ! Values from Table 8
    !

    P(1, 1)  = 0.003026697473481_DP
    P(2, 0)  = 0.007860700254745_DP
    P(3, 0)  = 0.016125364564267_DP
    P(1, 2)  = 0.003179530830731_DP
    P(2, 1)  = 0.006326248241119_DP

!-- Refit for proper graphen elastic constants
!    P(0, 2)  = g_graphene_P
!--

#ifdef ZERO_TABLES
    P = 0.0_DP
#endif

  endsubroutine rebo2_default_Pcc_table


  !**********************************************************************
  ! Construct the Pch-table
  !********************************************************************** 
  subroutine rebo2_default_Pch_table(P)
    implicit none

    real(DP), intent(out)  :: P(0:5, 0:5)

    ! ---

    P(:, :)     = 0.0_DP


    !
    ! Values from Table 8
    !

    P(1, 0)  =  0.2093367328250380_DP
    P(2, 0)  = -0.064449615432525_DP
    P(3, 0)  = -0.303927546346162_DP
    P(0, 1)  =  0.01_DP
    P(0, 2)  = -0.1220421462782555_DP
    P(1, 1)  = -0.1251234006287090_DP
    P(2, 1)  = -0.298905245783_DP
    P(0, 3)  = -0.307584705066_DP
    P(1, 2)  = -0.3005291724067579_DP

#ifdef ZERO_TABLES
    P = 0.0_DP
#endif

  endsubroutine rebo2_default_Pch_table


  !**********************************************************************
  ! Construct the Tcc-table
  !********************************************************************** 
  subroutine rebo2_default_Tcc_table(T)
    implicit none

    real(DP), intent(out)  :: T(0:4, 0:4, 0:9)

    ! ---

    T(:, :, :)     = 0.0_DP


    !
    ! Values from Table 5
    !

    T(2, 2, 0)  = -0.070280085_DP
!    T(2, 2, 8)  = -0.00809675_DP

    T(2, 2, 1:8)  = -0.00809675_DP

#ifdef ZERO_TABLES
    T = 0.0_DP
#endif
    
  endsubroutine rebo2_default_Tcc_table

endmodule rebo2_default_tables
