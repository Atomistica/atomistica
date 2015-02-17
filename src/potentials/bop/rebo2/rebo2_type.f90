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
! Declaration of the datatype which contains the materials parameters
! and local neighbor lists.
!**********************************************************************

  !
  ! Atom types
  !

  integer, parameter   :: rebo2_C_  = 1
  integer, parameter   :: rebo2_H_  = 3

  integer, parameter   :: possible_elements(2) = (/ rebo2_C_, rebo2_H_ /)

  integer, parameter   :: C_C = 1
  integer, parameter   :: C_H = 3
  integer, parameter   :: H_H = 6

  integer, parameter   :: H_ = 1
  integer, parameter   :: C_ = 6

  type g_coeff_t
     real(DP)  :: c(6, 3)
  endtype g_coeff_t

  public :: BOP_TYPE
  type BOP_TYPE

     character(MAX_EL_STR)  :: elements = "C,H"
     integer                :: els

     integer, allocatable   :: internal_el(:)

     !
     ! === THIS SECTION CONTAINS PARAMETERS ===
     !

#ifdef SCREENING

     !
     ! Screening parameters
     !

     real(DP)             :: Cmin  = 1.00_DP
     real(DP)             :: Cmax  = 2.00_DP

     real(DP)             :: screening_threshold  = log(1d-6)
     real(DP)             :: dot_threshold        = 1e-10

#endif

     !
     ! ============ C-C interaction ============
     !

     !
     ! Attractive function (C-C)
     ! [ Table 2 from Brenner2002 ]
     !

     real(DP)             :: cc_B1     = 12388.79197798_DP
     real(DP)             :: cc_B2     = 17.56740646509_DP
     real(DP)             :: cc_B3     = 30.71493208065_DP
     real(DP)             :: cc_beta1  = 4.7204523127_DP
     real(DP)             :: cc_beta2  = 1.4332132499_DP
     real(DP)             :: cc_beta3  = 1.3826912506_DP

     !
     ! Repulsive function (C-C)
     ! [ Table 2 from Brenner2002 ]
     !

     real(DP)             :: cc_Q     = 0.3134602960833_DP
     real(DP)             :: cc_A     = 10953.544162170_DP
     real(DP)             :: cc_alpha = 4.7465390606595_DP
     
     !
     ! g(cos(theta))  (C-C)
     ! [ Table 3 from Brenner2002 ]
     !
     
     !                                            ----------|----------|----------|----------|----------|----------|
     real(DP)             :: cc_g_theta(6)   = (/    -1.0_DP, -1.0_DP/2, -1.0_DP/3,    0.0_DP,  1.0_DP/2,    1.0_DP  /)
     real(DP)             :: cc_g_g1(6)      = (/      -0.01,   0.05280,   0.09733,   0.37545,    2.0014,       8.0  /)
     real(DP)             :: cc_g_dg1(6)     = (/    0.10400,   0.17000,   0.40000,       0.0,       0.0,       0.0  /)
     real(DP)             :: cc_g_d2g1(6)    = (/    0.00000,   0.37000,   1.98000,       0.0,       0.0,       0.0  /)
     real(DP)             :: cc_g_g2(6)      = (/        0.0,       0.0,   0.09733,  0.271856,  0.416335,       1.0  /)

     !
     ! ============ C-H interaction ============
     !

     !
     ! Attractive function (C-H)
     ! [ Table 7 from Brenner2002 ]
     !

     real(DP)             :: ch_B1     = 32.3551866587_DP
     real(DP)             :: ch_beta1  = 1.43445805925_DP

     !
     ! Repulsive function (C-H)
     ! [ Table 7 from Brenner2002 ]
     !

     real(DP)             :: ch_Q     = 0.340775728_DP
     real(DP)             :: ch_A     = 149.94098723_DP
     real(DP)             :: ch_alpha = 4.10254983_DP

     !
     ! ============ H-H interaction ============
     !

     !
     ! Attractive function (H-H)
     ! [ Table 6 from Brenner2002 ]
     !

     real(DP)             :: hh_B1     = 29.632593_DP
     real(DP)             :: hh_beta1  = 1.71589217_DP

     !
     ! Repulsive function (H-H)
     ! [ Table 6 from Brenner2002 ]
     !

     real(DP)             :: hh_Q     = 0.370471487045_DP
     real(DP)             :: hh_A     = 32.817355747_DP
     real(DP)             :: hh_alpha = 3.536298648_DP

     !
     ! g(cos(theta))  (H-H)
     ! [ Table 6 from Brenner2002 ]
     !

     !                                            ----------------|----------------|----------------|----------------|----------------|----------------|
     !  real(DP)             :: hh_g_theta(6)   = (/          -1.0_DP, -0.866025403_DP,       -1.0_DP/2,          0.0_DP,        1.0_DP/2,          1.0_DP  /)
     !  real(DP)             :: hh_g_g(6)       = (/        11.235870,       12.164186,       16.811574,       19.065124,       19.704059,       19.991787  /)
     
     ! This is directly from the MD code of Brenner's group...
     ! ...there is no way to get this from the paper.

     real(DP)  :: SPGH(6,3) = &
          reshape( &
          (/ 270.467795364007301_DP ,1549.701314596994564_DP  &
          ,3781.927258631323866_DP,4582.337619544424228_DP,2721.538161662818368_DP, &
          630.658598136730774_DP,16.956325544514659_DP,-21.059084522755980_DP, &
          -102.394184748124742_DP,-210.527926707779059_DP,-229.759473570467513_DP, &
          -94.968528666251945_DP,19.065031149937783_DP,2.017732531534021_DP, &
          -2.566444502991983_DP,3.291353893907436_DP,-2.653536801884563_DP, &
          0.837650930130006_DP /), &
          (/ 6, 3 /) )
     integer  :: IGH(25) = &
          (/ 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, &
             2, 2, 2, 2, &
             1, 1, 1 /)

!     integer   :: IGH(25) = &
!          (/18*3,4*2,3*1/)

!     DATA SPGH / 270.467795364007301,1549.701314596994564  &
!          ,3781.927258631323866,4582.337619544424228,2721.538161662818368, &
!          630.658598136730774,16.956325544514659,-21.059084522755980, &
!          -102.394184748124742,-210.527926707779059,-229.759473570467513, &
!          -94.968528666251945,19.065031149937783,2.017732531534021, &
!          -2.566444502991983,3.291353893907436,-2.653536801884563, &
!          0.837650930130006/

!     DATA IGH/18*3,4*2,3*1/

     !  real(DP)             :: hh_g_coeff(6)

     !
     ! Other parameters
     !
     
     real(DP)             :: hhh_lambda = 4.0_DP

     real(DP)             :: cc_re = 1.4_DP
     real(DP)             :: ch_re = 1.09_DP
     real(DP)             :: hh_re = 0.7415886997_DP

#ifdef SCREENING
     real(DP)             :: cc_in_r1 = 1.95_DP
     real(DP)             :: cc_in_r2 = 2.25_DP

     real(DP)             :: cc_ar_r1 = 2.179347_DP
     real(DP)             :: cc_ar_r2 = 2.819732_DP
     real(DP)             :: cc_bo_r1 = 1.866344_DP
     real(DP)             :: cc_bo_r2 = 2.758372_DP
     real(DP)             :: cc_nc_r1 = 1.217335_DP
     real(DP)             :: cc_nc_r2 = 4.000000_DP
#else
     real(DP)             :: cc_in_r1 = 1.70_DP
     real(DP)             :: cc_in_r2 = 2.00_DP
#endif

     !  --- Using the original dihedral terms
     !  real(DP)             :: cc_ar_r1 = 2.157379_DP
     !  real(DP)             :: cc_ar_r2 = 2.817738_DP
     !  real(DP)             :: cc_bo_r1 = 1.883182_DP
     !  real(DP)             :: cc_bo_r2 = 2.713410_DP
     !  real(DP)             :: cc_nc_r1 = 1.128100_DP
     !  real(DP)             :: cc_nc_r2 = 4.000000_DP

     !  --- Using the alt. dihedral terms

     real(DP)             :: ch_r1 = 1.30_DP
     real(DP)             :: ch_r2 = 1.80_DP

     real(DP)             :: hh_r1 = 1.10_DP
     real(DP)             :: hh_r2 = 1.70_DP

     logical(C_BOOL)      :: with_dihedral = .false.

     !
     ! === THIS SECTION CONTAINS *DERIVED* DATA
     !

#ifdef SCREENING
     real(DP)  :: dC
     real(DP)  :: C_dr_cut
#endif

     type(g_coeff_t)  :: cc_g1_coeff
     type(g_coeff_t)  :: cc_g2_coeff

     real(DP)  :: conalp,conear(6,6)
     real(DP)  :: conpe(3),conan(3),conpf(3)
     real(DP)  :: cut_in_h(10), cut_in_h2(10), cut_in_m(10), cut_in_l(10)
#ifdef SCREENING
     real(DP)  :: cut_ar_h(10), cut_ar_h2(10), cut_ar_m(10), cut_ar_l(10)
     real(DP)  :: cut_bo_h(10), cut_bo_h2(10), cut_bo_m(10), cut_bo_l(10)
     real(DP)  :: cut_nc_h(10), cut_nc_h2(10), cut_nc_m(10), cut_nc_l(10)
#endif

     real(DP)  :: max_cut_sq(10)

     !
     ! Lookup tables
     !

     type(table3d_t)  :: Fcc
     type(table3d_t)  :: Fch
     type(table3d_t)  :: Fhh

     type(table2d_t)  :: Pcc
     type(table2d_t)  :: Pch

     type(table3d_t)  :: Tcc

     !
     ! Splines
     !

     integer                :: spl_n = 1000
     real(DP)               :: spl_x0 = 0.1
     type(simple_spline_t)  :: spl_VA(10)
     type(simple_spline_t)  :: spl_VR(10)
     type(simple_spline_t)  :: spl_fCin(10)
#ifdef SCREENING
     type(simple_spline_t)  :: spl_fCar(10)
     type(simple_spline_t)  :: spl_fCbo(10)
     type(simple_spline_t)  :: spl_fCnc(10)
#endif

     !
     ! Counters
     !

     logical  :: neighbor_list_allocated  = .false.
     logical  :: tables_allocated         = .false.
     integer  :: it                       = 0

     !
     ! Quick and dirty hack
     !

!     real(DP)  :: g_graphene_P = -0.35_DP
!     real(DP)  :: g_graphene_F = -0.39873_DP

#ifdef NUM_NEIGHBORS
     ! Precomputed number of neighbors
     real(DP), allocatable  :: nn(:, :)
#endif

     !
     ! Internal neighbor lists
     !

     integer, allocatable   :: neb_seed(:)
     integer, allocatable   :: neb_last(:)

     integer, allocatable   :: neb(:)
     integer, allocatable   :: nbb(:)
#ifndef LAMMPS
     integer, allocatable   :: dcell(:)
#endif

     integer, allocatable   :: bndtyp(:)
     real(DP), allocatable  :: bndlen(:)
     real(DP), allocatable  :: bndnm(:, :)
     real(DP), allocatable  :: cutfcnar(:), cutdrvar(:)

#ifdef SCREENING
     real(DP), allocatable  :: cutfcnbo(:), cutdrvbo(:)
     real(DP), allocatable  :: cutfcnnc(:), cutdrvnc(:)
     ! "screened" neighbor list (all neighbors of a bond which sit in the
     ! screening cutoff)
     integer, allocatable   :: sneb_seed(:)
     integer, allocatable   :: sneb_last(:)
     integer, allocatable   :: sneb(:)
     integer(NEIGHPTR_T), allocatable   :: sbnd(:)

     ! for force calculation
     real(DP), allocatable  :: sfacbo(:)
     real(DP), allocatable  :: sfacnc(:)

     real(DP), allocatable  :: cutdrarik(:), cutdrarjk(:)
     real(DP), allocatable  :: cutdrboik(:), cutdrbojk(:)
     real(DP), allocatable  :: cutdrncik(:), cutdrncjk(:)
#endif

     !
     ! From the input file
     !
     
     logical(C_BOOL) :: zero_tables = .false.

     real(DP)  :: in_Fcc(0:4, 0:4, 0:9) = 0.0_DP
     real(DP)  :: in_dFdi(0:4, 0:4, 0:9) = 0.0_DP
     real(DP)  :: in_dFdj(0:4, 0:4, 0:9) = 0.0_DP
     real(DP)  :: in_dFdk(0:4, 0:4, 0:9) = 0.0_DP

     real(DP)  :: in_Fch(0:4, 0:4, 0:9) = 0.0_DP
     real(DP)  :: in_Fhh(0:4, 0:4, 0:9) = 0.0_DP

     real(DP)  :: in_Pcc(0:5, 0:5) = 0.0_DP
     real(DP)  :: in_Pch(0:5, 0:5) = 0.0_DP
     real(DP)  :: in_Tcc(0:4, 0:4, 0:9) = 0.0_DP

  endtype BOP_TYPE

!  type(BOP_TYPE), save  :: BOP_NAME##_default_parameters

  public :: init
  interface init
     module procedure INIT_FUNC
  endinterface

  public :: del
  interface del
     module procedure DEL_FUNC
  endinterface

  public :: bind_to
  interface bind_to
     module procedure BIND_TO_FUNC
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure COMPUTE_FUNC
  endinterface

  public :: register, REGISTER_FUNC
  interface register
     module procedure REGISTER_FUNC
  endinterface register

  public :: rebo2_default_Fcc_table, rebo2_default_Fch_table, rebo2_default_Fhh_table
  public :: rebo2_default_Pcc_table, rebo2_default_Pch_table, rebo2_default_Tcc_table
