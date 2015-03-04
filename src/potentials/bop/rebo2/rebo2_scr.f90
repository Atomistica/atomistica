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
!    shared:directory
!    dependencies:rebo2_default_tables.f90
!    classtype:rebo2_scr_t classname:Rebo2Scr interface:potentials
!    features:per_at,per_bond
! @endmeta

!>
!! The screened second generation reactive empirical bond-order potential
!! (REBO2+S)
!!
!! The screened second generation reactive empirical bond-order potential
!! (REBO2+S)
!! See: Brenner et al., J. Phys.: Condens. Matter 14, 783 (2002)
!! Pastewka, Pou, Perez, Gumbsch, Moseler, Phys. Rev. B 78, 161402(R) (2008)
!<

#include "macros.inc"
#include "filter.inc"

module rebo2_scr
  use, intrinsic :: iso_c_binding

  use supplib

  use particles
  use filter
  use neighbors

  use table2d
  use table3d

  use rebo2_default_tables

  implicit none

  private

#define SCREENING

#define ALT_DIHEDRAL

#define NUM_NEIGHBORS

#define BOP_NAME             rebo2_scr
#define BOP_NAME_STR         "rebo2_scr"
#define BOP_STR              "Rebo2Scr"
#define BOP_KERNEL           rebo2_scr_kernel
#define BOP_TYPE             rebo2_scr_t

#define REGISTER_FUNC        rebo2_scr_register
#define INIT_FUNC            rebo2_scr_init
#define INTERNAL_INIT_FUNC   rebo2_scr_internal_init
#define DEL_FUNC             rebo2_scr_del
#define COMPUTE_FUNC         rebo2_scr_energy_and_forces

#include "rebo2_type.f90"

contains

#include "rebo2_db.f90"

#include "rebo2_module.f90"

#include "bop_kernel_rebo2.f90"

#include "rebo2_func.f90"

#include "rebo2_registry.f90"

endmodule rebo2_scr
