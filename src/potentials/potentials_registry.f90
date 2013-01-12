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
#include "macros.inc"

!>
!! Registry
!!
!! Contains function to register each classes property
!<
module potentials_registry
  use libAtoms_module

  use particles
  use neighbors
  use filter

#include "potentials.inc"

  use bop_registry

  implicit none

endmodule potentials_registry

