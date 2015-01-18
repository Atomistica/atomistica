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
!   classtype:brenner_scr_t classname:BrennerScr interface:potentials
!   features:mask,per_at,per_bond
! @endmeta

!>
!! Screened Abell-Tersoff-Brenner type potentials
!!
!! Screened Abell-Tersoff-Brenner type potentials with the Morse-style pair
!! terms used by Brenner. Note: This potential does not contain the correction
!! tables for treatment of pi-orbitals, etc.
!<

#include "macros.inc"

module brenner_scr
  use libAtoms_module

  use ptrdict

  use logging

  use timer

  use particles
  use neighbors

  implicit none

  private

#define SCREENING
#define CUTOFF_T             exp_cutoff_t

#define BRENNER_MAX_REF      BRENNER_SCR_MAX_REF
#define BRENNER_MAX_EL       BRENNER_SCR_MAX_EL
#define BRENNER_MAX_PAIRS    BRENNER_SCR_MAX_PAIRS

#define BOP_NAME             brenner_scr
#define BOP_NAME_STR         "brenner_scr"
#define BOP_STR              "BrennerScr"
#define BOP_KERNEL           brenner_scr_kernel
#define BOP_TYPE             brenner_scr_t
#define BOP_DB               brenner_db_scr
#define BOP_DB_TYPE          brenner_db_scr_t

#define REGISTER_FUNC        brenner_scr_register
#define INIT_FUNC            brenner_scr_init
#define DEL_FUNC             brenner_scr_del
#define BIND_TO_FUNC         brenner_scr_bind_to
#define COMPUTE_FUNC         brenner_scr_energy_and_forces

#include "brenner_params.f90"

#include "brenner_type.f90"

contains

#include "brenner_module.f90"

#include "../bop_kernel.f90"

#include "brenner_func.f90"

#include "brenner_registry.f90"

endmodule brenner_scr
