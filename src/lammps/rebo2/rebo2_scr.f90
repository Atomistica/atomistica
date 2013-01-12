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
!    shared:directory
!    dependencies:rebo2_default_tables.f90
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
  use libAtoms_module

  use logging
  use simple_spline

  use table2d
  use table3d

  use rebo2_default_tables

#ifdef _MP
  use parallel_3d
#endif

  implicit none

  private

  integer, parameter :: MAX_EL_STR = 16

#define SCREENING

#define ALT_DIHEDRAL

#define NUM_NEIGHBORS

#define BOP_NAME             rebo2_scr
#define BOP_NAME_STR         "rebo2_scr"
#define BOP_STR              "Rebo2Scr"
#define BOP_KERNEL           rebo2_scr_kernel
#define BOP_TYPE             rebo2_scr_t

#define CREATE_FUNC          rebo2_scr_create
#define DESTROY_FUNC         rebo2_scr_destroy
#define GET_CUTOFF_FUNC      rebo2_scr_get_cutoff

#define REGISTER_FUNC        rebo2_scr_register
#define INIT_FUNC            rebo2_scr_init
#define INIT_DEFAULT_FUNC    rebo2_scr_init
#define DEL_FUNC             rebo2_scr_del
#define BIND_TO_FUNC         rebo2_scr_init
#define COMPUTE_FUNC         rebo2_scr_energy_and_forces

#include "rebo2_type.f90"

contains

#include "rebo2_db.f90"

#include "rebo2_module.f90"

#include "bop_kernel_rebo2.f90"

#include "rebo2_func.f90"

endmodule rebo2_scr
