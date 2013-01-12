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
!   classtype:albe_bop_scr_t classname:AlbeBOPScr interface:potentials
! @endmeta

!>
!! Screened Albe-style Tersoff-Brenner potentials
!!
!! Screened Albe-style Tersoff-Brenner potentials
!<

#include "macros.inc"

module albe_bop_scr
  use libAtoms_module

  use logging

  use timer

  use particles
  use neighbors

  implicit none

  private

#define SCREENING
#define EXP_CUTOFF

#define ALBE_BOP_MAX_REF     ALBE_BOP_SCR_MAX_REF
#define ALBE_BOP_MAX_EL      ALBE_BOP_SCR_MAX_EL
#define ALBE_BOP_MAX_PAIRS   ALBE_BOP_SCR_MAX_PAIRS

#define BOP_NAME             albe_bop_scr
#define BOP_NAME_STR         "albe_bop_scr"
#define BOP_STR              "AlbeBOPScr"
#define BOP_KERNEL           albe_bop_scr_kernel
#define BOP_TYPE             albe_bop_scr_t
#define BOP_DB               albe_bop_db_scr
#define BOP_DB_TYPE          albe_bop_db_scr_t

#define REGISTER_FUNC        albe_bop_scr_register
#define INIT_FUNC            albe_bop_scr_init
#define DEL_FUNC             albe_bop_scr_del
#define GET_CUTOFF_FUNC      albe_bop_scr_get_cutoff
#define BIND_TO_FUNC         albe_bop_scr_bind_to
#define COMPUTE_FUNC         albe_bop_scr_energy_and_forces
#define FORCE_FUNC           albe_bop_scr_force

#include "albe_bop_params.f90"

#include "albe_bop_type.f90"

contains

#include "albe_bop_module.f90"

#include "../bop_kernel.f90"

#include "albe_bop_func.f90"

endmodule albe_bop_scr
