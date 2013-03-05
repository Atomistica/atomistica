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
!    public:directory
!    classtype:tersoff_t classname:Tersoff interface:potentials
! @endmeta

!>
!! Tersoff's potential
!!
!! Tersoff's potential
!! See: Tersoff, Phys. Rev. Lett. 56, 632 (1986)
!! Tersoff, Phys. Rev. Lett. 61, 2879 (1988)
!! Tersoff, Phys. Rev. B 37, 6991 (1988)
!! Tersoff, Phys. Rev. B 38, 9902 (1988)
!! Tersoff, Phys. Rev. B 39, 5566 (1989)
!<

#include "macros.inc"

module tersoff
  use supplib

  use particles
  use neighbors

  implicit none

  private

#define BOP_NAME             tersoff
#define BOP_NAME_STR         "tersoff"
#define BOP_STR              "Tersoff"
#define BOP_KERNEL           tersoff_kernel
#define BOP_TYPE             tersoff_t
#define BOP_DB_TYPE          tersoff_db_t

#define REGISTER_FUNC        tersoff_register
#define INIT_FUNC            tersoff_init
#define DEL_FUNC             tersoff_del
#define GET_CUTOFF_FUNC      tersoff_get_cutoff
#define BIND_TO_FUNC         tersoff_bind_to
#define COMPUTE_FUNC         tersoff_energy_and_forces

#include "tersoff_params.f90"

#include "tersoff_type.f90"

contains

#include "tersoff_module.f90"

#include "../bop_kernel.f90"

#include "tersoff_func.f90"

#include "tersoff_registry.f90"

endmodule tersoff
