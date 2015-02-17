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

  !>
  !! The BOP class
  !<
  public :: BOP_TYPE
  type BOP_TYPE

     !
     ! Current database
     !

     type(BOP_DB_TYPE)  :: db = Juslin_J_Appl_Phys_98_123520_WCH

     !
     ! String to a reference
     !
     character(JUSLIN_MAX_REF)  :: ref = "*"

     integer       :: Z2db(MAX_Z)

     !
     ! Precomputed constants
     !

     real(DP)  :: c_sq(JUSLIN_MAX_PAIRS)
     real(DP)  :: d_sq(JUSLIN_MAX_PAIRS)
     real(DP)  :: c_d(JUSLIN_MAX_PAIRS) 

     real(DP)  :: VR_f(JUSLIN_MAX_PAIRS)
     real(DP)  :: expR(JUSLIN_MAX_PAIRS)
     real(DP)  :: VA_f(JUSLIN_MAX_PAIRS)
     real(DP)  :: expA(JUSLIN_MAX_PAIRS)

     !
     ! Cutoff parameters
     !

     real(DP)  :: cut_in_h(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_in_h2(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_in_l(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_in_fca(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_in_fc(JUSLIN_MAX_PAIRS)

#ifdef SCREENING

! The other cutoff are identical!
#define cut_ar_h  cut_out_h

     real(DP)  :: cut_out_h(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_out_l(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_out_fca(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_out_fc(JUSLIN_MAX_PAIRS)

     real(DP)  :: cut_bo_h(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_bo_l(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_bo_fca(JUSLIN_MAX_PAIRS)
     real(DP)  :: cut_bo_fc(JUSLIN_MAX_PAIRS)

     real(DP)  :: max_cut_sq(JUSLIN_MAX_PAIRS)

     real(DP)  :: Cmin(JUSLIN_MAX_PAIRS)
     real(DP)  :: Cmax(JUSLIN_MAX_PAIRS)
     real(DP)  :: dC(JUSLIN_MAX_PAIRS)
     real(DP)  :: C_dr_cut(JUSLIN_MAX_PAIRS)

     real(DP)  :: screening_threshold  = log(1d-6)
     real(DP)  :: dot_threshold        = 1e-10
#endif

     !
     ! Bond-order stuff
     !

     real(DP)  :: bo_exp(JUSLIN_MAX_PAIRS)
     real(DP)  :: bo_fac(JUSLIN_MAX_PAIRS)
     real(DP)  :: bo_exp1(JUSLIN_MAX_PAIRS)

     !
     ! Counters
     !

     logical  :: neighbor_list_allocated  = .false.
     integer  :: it                       = 0

     !
     ! Internal neighbor lists
     !

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
     ! "screened" neighbor list (all neighbors of a bond which sit in the
     ! screening cutoff)
     integer, allocatable   :: sneb_seed(:)
     integer, allocatable   :: sneb_last(:)
     integer, allocatable   :: sneb(:)
     integer, allocatable   :: sbnd(:)

     ! for force calculation
     real(DP), allocatable  :: sfacbo(:)

     real(DP), allocatable  :: cutdrarik(:), cutdrarjk(:)
     real(DP), allocatable  :: cutdrboik(:), cutdrbojk(:)
#endif

  endtype BOP_TYPE


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

  public :: register
  interface register
     module procedure REGISTER_FUNC
  endinterface register


