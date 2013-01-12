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

! @meta
!   public:directory
!   classtype:brenner_scr_t classname:BrennerScr interface:potentials
! @endmeta

!>
!! Screened Tersoff-Brenner type potentials
!!
!! Screened Tersoff-Brenner type potentials with the Morse-style pair terms
!! introduced by Brenner. Note: This potential does not contain the correction
!! tables for treatment of pi-orbitals, etc.
!<

#include "macros.inc"

module brenner_scr
  use libAtoms_module

  use logging

  use timer

  use particles
  use neighbors

  implicit none

  private

#define SCREENING
#define EXP_CUTOFF

#define BRENNER_MAX_REF     BRENNER_SCR_MAX_REF
#define BRENNER_MAX_EL      BRENNER_SCR_MAX_EL
#define BRENNER_MAX_PAIRS   BRENNER_SCR_MAX_PAIRS

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
#define GET_CUTOFF_FUNC      brenner_scr_get_cutoff
#define BIND_TO_FUNC         brenner_scr_bind_to
#define COMPUTE_FUNC         brenner_scr_energy_and_forces
#define FORCE_FUNC           brenner_scr_force

#include "brenner_params.f90"

#include "brenner_type.f90"

contains

#include "brenner_module.f90"

#include "../bop_kernel.f90"

#include "brenner_func.f90"

endmodule brenner_scr
