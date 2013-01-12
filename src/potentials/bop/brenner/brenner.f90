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
!   classtype:albe_bop_t classname:AlbeBOP interface:potentials
! @endmeta

!>
!! Albe-style Tersoff-Brenner potentials
!!
!! Albe-style Tersoff-Brenner potentials
!<

#include "macros.inc"

module albe_bop
  use libAtoms_module

  use logging

  use timer

  use particles
  use neighbors

  implicit none

  private

#define BOP_NAME             albe_bop
#define BOP_NAME_STR         "albe_bop"
#define BOP_STR              "AlbeBOP"
#define BOP_KERNEL           albe_bop_kernel
#define BOP_TYPE             albe_bop_t
#define BOP_DB               albe_bop_db
#define BOP_DB_TYPE          albe_bop_db_t

#define REGISTER_FUNC        albe_bop_register
#define INIT_FUNC            albe_bop_init
#define DEL_FUNC             albe_bop_del
#define GET_CUTOFF_FUNC      albe_bop_get_cutoff
#define BIND_TO_FUNC         albe_bop_bind_to
#define COMPUTE_FUNC         albe_bop_energy_and_forces
#define FORCE_FUNC           albe_bop_force

#include "albe_bop_params.f90"

#include "albe_bop_type.f90"

contains

#include "albe_bop_module.f90"

#include "../bop_kernel.f90"

#include "albe_bop_func.f90"

endmodule albe_bop
