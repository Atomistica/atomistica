!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
!! and others. See the AUTHORS file in the top-level Atomistica directory.
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
!   shared
!   classtype:variable_charge_t classname:VariableCharge interface:potentials
! @endmeta

!> Charge transfer
!!
!! Module for the evaluation of charge transfer. Does only work
!! in combination with a Coulomb (and possibly charge overlap) module. Uses
!! either a conjugate gradient algorithm, anderson mixing or
!! Car-Parrinello to optimize the charges.
!!
!! The method is described in:
!!
!!   A. K. Rappe and W. A. Goddard III, J. Phys. Chem. 95, 3358 (1991)
!!   F. H. Streitz and J. W. Mintmire, Phys. Rev. B 50, 11996 (1994)
