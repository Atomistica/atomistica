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

!>
!! General tight-binding datastructure
!<

#include "macros.inc"

module dense_hamiltonian_type
  use, intrinsic :: iso_c_binding

  implicit none

  private

  !
  ! This type contains all the information about the tight-binding
  ! part of the simulation, i.e., Hamiltionians, ...
  !
  ! This datastructure is fully C interoperable. See accompanying
  ! header file dense_hamiltonian.h
  !

  ! IF YOU MODIFY THIS STRUCTURE, *ALWAYS* ALSO MODIFY THE CORRESPONDING
  ! STRUCTURE IN dense_hamiltonian.h
  public :: dense_hamiltonian_t
  type, bind(C) :: dense_hamiltonian_t

     !
     ! References to other objects
     !
     ! IN REALITY:
     ! type(particles_t)  :: p
     ! type(materials_t)  :: mat
     !

     type(C_PTR)        :: p = C_NULL_PTR
     type(C_PTR)        :: mat = C_NULL_PTR

     !
     ! Which elements to we treat? I.e. the filter
     !

     integer(C_INT)     :: f

     !
     ! Number of kpoints, orbitals
     !

     integer(C_INT)     :: nat          ! number of atoms
     integer(C_INT)     :: nk = 1       ! number of k-points/spins
     integer(C_INT)     :: norb         ! number of orbitals
     integer(C_INT)     :: norbloc      ! number of orbitals (local to this process)

     type(C_PTR)        :: el           ! list of elements (from particles_t)

     real(C_DOUBLE)     :: mu   ! Fermi level

     real(C_DOUBLE)     :: cutoff

     ! IN REALITY:
     ! WF_T(DP)  :: H(:, :, :)
     ! WF_T(DP)  :: S(:, :, :)
     ! WF_T(DP)  :: rho(:, :, :)
     ! WF_T(DP)  :: e(:, :, :)

     type(C_PTR)        :: H = C_NULL_PTR     ! Hamiltonian (for each k-point)
     type(C_PTR)        :: S = C_NULL_PTR     ! the overlap matrix
     type(C_PTR)        :: rho = C_NULL_PTR   ! The density matrix (rho_ll)
     type(C_PTR)        :: e = C_NULL_PTR     ! H_rl * rho_ll

     !
     ! Band-structure and repulsive energies
     !

     real(C_DOUBLE)     :: ebs
     real(C_DOUBLE)     :: erep

     !
     ! Additional particle information
     !
     ! IN REALITY
     ! real(DP)              :: n(:)
     ! type(notb_element_t)  :: at(:)
     !

     type(C_PTR)        :: n = C_NULL_PTR    ! Number of electrons, charges
     type(C_PTR)        :: at = C_NULL_PTR

  endtype dense_hamiltonian_t

  interface
     subroutine dense_hamiltonian_allocate(this, nat, norb) bind(C)
       use, intrinsic :: iso_c_binding
       type(C_PTR),    value  :: this
       integer(C_INT), value  :: nat, norb
     endsubroutine dense_hamiltonian_allocate

     subroutine dense_hamiltonian_deallocate(this) bind(C)
       use, intrinsic :: iso_c_binding
       type(C_PTR),    value  :: this
     endsubroutine dense_hamiltonian_deallocate
  endinterface

  public :: dense_hamiltonian_allocate, dense_hamiltonian_deallocate

endmodule dense_hamiltonian_type
