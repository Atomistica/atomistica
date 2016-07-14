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
!   classtype:kumagai_t classname:Kumagai interface:potentials
!   features:mask,per_at,per_bond
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

#define CUTOFF_T             trig_off_t

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
