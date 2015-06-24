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

! @meta
!   public:directory
!   classtype:juslin_scr_t classname:JuslinScr interface:potentials
!   features:mask,per_at,per_bond
! @endmeta

!>
!! Screened version of Juslin's W-C-H potential
!!
!! Implementation of the specific functional forms of Juslin's W-C-H potential,
!! including screening functions.
!! See: Juslin, Erhart, Traskelin, Nord, Henriksson, Nordlund, Salonen, Albe,
!! J. Appl. Phys. 98, 123520 (2005)
!<

#include "macros.inc"

module juslin_scr
  use libAtoms_module

  use ptrdict

  use logging

  use timer

  use particles
  use neighbors

  implicit none

  private

#define SCREENING
#define EXP_BOP

#define JUSLIN_MAX_REF     JUSLIN_SCR_MAX_REF
#define JUSLIN_MAX_EL      JUSLIN_SCR_MAX_EL
#define JUSLIN_MAX_PAIRS   JUSLIN_SCR_MAX_PAIRS

#define BOP_NAME             juslin_scr
#define BOP_NAME_STR         "juslin_scr"
#define BOP_STR              "JuslinScr"
#define BOP_KERNEL           juslin_scr_kernel
#define BOP_TYPE             juslin_scr_t
#define BOP_DB               juslin_db_scr
#define BOP_DB_TYPE          juslin_db_scr_t

#define REGISTER_FUNC        juslin_scr_register
#define INIT_FUNC            juslin_scr_init
#define DEL_FUNC             juslin_scr_del
#define BIND_TO_FUNC         juslin_scr_bind_to
#define COMPUTE_FUNC         juslin_scr_energy_and_forces
#define FORCE_FUNC           juslin_scr_force

#include "juslin_params.f90"

#include "juslin_type.f90"

contains

#include "juslin_module.f90"

#include "../bop_kernel.f90"

#include "juslin_func.f90"

#include "juslin_registry.f90"

endmodule juslin_scr
