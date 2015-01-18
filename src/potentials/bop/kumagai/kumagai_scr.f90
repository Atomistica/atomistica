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
!   classtype:kumagai_scr_t classname:KumagaiScr interface:potentials
!   features:mask,per_at,per_bond
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
#define CUTOFF_T             exp_cutoff_t

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
