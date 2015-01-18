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

module supplib
  use c_f
  use error_module
  use MPI_context_module
  use units_module
  use periodictable_module
  use io
!  use histogram1d_module
  use linearalgebra
  use logging
  use misc
!  use rng
!  use math
#ifndef LAMMPS
  use data
#endif
!  use special_functions
  use simple_spline
  use nonuniform_spline
  use timer
  use tls
  use ptrdict
!  use signal_handler
  use cutoff
  use histogram1d_module
endmodule supplib
