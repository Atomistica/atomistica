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
!    classtype:kumagai_scr_t classname:KumagaiScr interface:potentials
! @endmeta

!>
!! Screened Kumagai-Izumi-Hara-Sakai potential
!!
!! Screened Kumagai-Izumi-Hara-Sakai potential
!! See: Kumagai, Izumi, Hara, Sakai, Comp. Mater. Sci. 39, 457 (2007)
!! Pastewka, Klemenz, Gumbsch, Moseler, arXiv:1301.2142
!<

#include "macros.inc"

module kumagai_scr
  use supplib

  use particles
  use neighbors

  implicit none

  private

#define SCREENING
#define EXP_CUTOFF

#define KUMAGAI_MAX_REF      KUMAGAI_SCR_MAX_REF
#define KUMAGAI_MAX_EL       KUMAGAI_SCR_MAX_EL
#define KUMAGAI_MAX_PAIRS    KUMAGAI_SCR_MAX_PAIRS

#define BOP_NAME             kumagai_scr
#define BOP_NAME_STR         "kumagai_scr"
#define BOP_STR              "KumagaiScr"
#define BOP_KERNEL           kumagai_scr_kernel
#define BOP_TYPE             kumagai_scr_t
#define BOP_DB_TYPE          kumagai_scr_db_t

#define REGISTER_FUNC        kumagai_scr_register
#define INIT_FUNC            kumagai_scr_init
#define DEL_FUNC             kumagai_scr_del
#define BIND_TO_FUNC         kumagai_scr_bind_to
#define COMPUTE_FUNC         kumagai_scr_energy_and_forces

#include "kumagai_params.f90"

#include "kumagai_type.f90"

contains

#include "kumagai_module.f90"

#include "../bop_kernel.f90"

#include "kumagai_func.f90"

#include "kumagai_registry.f90"

endmodule kumagai_scr
