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
!    classtype:tersoff_scr_t classname:TersoffScr interface:potentials
! @endmeta

!>
!! Screened Tersoff potential
!!
!! Screened Tersoff potential
!! See: Tersoff, Phys. Rev. Lett. 56, 632 (1986)
!! Tersoff, Phys. Rev. Lett. 61, 2879 (1988)
!! Tersoff, Phys. Rev. B 37, 6991 (1988)
!! Tersoff, Phys. Rev. B 38, 9902 (1988)
!! Tersoff, Phys. Rev. B 39, 5566 (1989)
!! Pastewka, Klemenz, Gumbsch, Moseler, arXiv:1301.2142
!<

#include "macros.inc"

module tersoff_scr
  use libAtoms_module

  use ptrdict

  use logging
  use timer

  use particles
  use neighbors

  implicit none

  private

#define SCREENING
#define EXP_CUTOFF

#define TERSOFF_MAX_REF      TERSOFF_SCR_MAX_REF
#define TERSOFF_MAX_EL       TERSOFF_SCR_MAX_EL
#define TERSOFF_MAX_PAIRS    TERSOFF_SCR_MAX_PAIRS

#define BOP_NAME             tersoff_scr
#define BOP_NAME_STR         "tersoff_scr"
#define BOP_STR              "TersoffScr"
#define BOP_KERNEL           tersoff_scr_kernel
#define BOP_TYPE             tersoff_scr_t
#define BOP_DB_TYPE          tersoff_scr_db_t

#define REGISTER_FUNC        tersoff_scr_register
#define INIT_FUNC            tersoff_scr_init
#define DEL_FUNC             tersoff_scr_del
#define GET_CUTOFF_FUNC      tersoff_scr_get_cutoff
#define BIND_TO_FUNC         tersoff_scr_bind_to
#define COMPUTE_FUNC         tersoff_scr_energy_and_forces

#include "tersoff_params.f90"

#include "tersoff_type.f90"

contains

#include "tersoff_module.f90"

#include "../bop_kernel.f90"

#include "tersoff_func.f90"

endmodule tersoff_scr
