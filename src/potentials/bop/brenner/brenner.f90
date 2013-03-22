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
!   classtype:brenner_t classname:Brenner interface:potentials
! @endmeta

!>
!! Abell-Tersoff-Brenner type potentials
!!
!! Abell-Tersoff-Brenner type potentials with the Morse-style pair terms
!! used by Brenner. Note: This potential does not contain the correction 
!! tables for treatment of pi-orbitals, etc.
!<

#include "macros.inc"

module brenner
  use libAtoms_module

  use ptrdict

  use logging

  use timer

  use particles
  use neighbors

  implicit none

  private

#define BOP_NAME             brenner
#define BOP_NAME_STR         "brenner"
#define BOP_STR              "Brenner"
#define BOP_KERNEL           brenner_kernel
#define BOP_TYPE             brenner_t
#define BOP_DB               brenner_db
#define BOP_DB_TYPE          brenner_db_t

#define REGISTER_FUNC        brenner_register
#define INIT_FUNC            brenner_init
#define DEL_FUNC             brenner_del
#define BIND_TO_FUNC         brenner_bind_to
#define COMPUTE_FUNC         brenner_energy_and_forces

#include "brenner_params.f90"

#include "brenner_type.f90"

contains

#include "brenner_module.f90"

#include "../bop_kernel.f90"

#include "brenner_func.f90"

#include "brenner_registry.f90"

endmodule brenner
