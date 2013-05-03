!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================
! @meta
!    public:directory
!    classtype:kumagai_t classname:Kumagai interface:potentials
! @endmeta

!>
!! Kumagai-Izumi-Hara-Sakai potential
!!
!! Kumagai-Izumi-Hara-Sakai potential
!! See: Kumagai, Izumi, Hara, Sakai, Comp. Mater. Sci. 39, 457 (2007)
!<

#include "macros.inc"

module kumagai
  use supplib

  use particles
  use neighbors

  implicit none

  private

#define BOP_NAME             kumagai
#define BOP_NAME_STR         "kumagai"
#define BOP_STR              "Kumagai"
#define BOP_KERNEL           kumagai_kernel
#define BOP_TYPE             kumagai_t
#define BOP_DB_TYPE          kumagai_db_t

#define REGISTER_FUNC        kumagai_register
#define INIT_FUNC            kumagai_init
#define DEL_FUNC             kumagai_del
#define BIND_TO_FUNC         kumagai_bind_to
#define COMPUTE_FUNC         kumagai_energy_and_forces

#include "kumagai_params.f90"

#include "kumagai_type.f90"

contains

#include "kumagai_module.f90"

#include "../bop_kernel.f90"

#include "kumagai_func.f90"

#include "kumagai_registry.f90"

endmodule kumagai
