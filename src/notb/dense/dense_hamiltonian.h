/* ======================================================================
   Atomistica - Interatomic potential library and molecular dynamics code
   https://github.com/Atomistica/atomistica

   Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
   See the AUTHORS file in the top-level Atomistica directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */

#ifndef __DENSE_HAMILTONIAN_H
#define __DENSE_HAMILTONIAN_H

#include "materials.h"

/*
 * This type contains all the information about the tight-binding
 * part of the simulation, i.e., Hamiltionians, ...
 *
 * This datastructure is fully F90 interoperable. See accompanying
 * file dense_hamiltonian_type.f90
 */

/*
 * IF YOU MODIFY THIS STRUCTURE, *ALWAYS* ALSO MODIFY THE CORRESPONDING
 * STRUCTURE IN dense_hamiltonian_type.f90
 */

struct dense_hamiltonian_t {
  /*
   * References to other objects
   *
   * IN FORTRAN:
   * type(particles_t)  :: p
   * type(materials_t)  :: mat
   */
  
  void *p, *mat;

  /*
   * Which elements to we treat? I.e. the filter.
   */

  int f;

  /*
   * Number of kpoints, orbitals
   */

  int nat; /* number of atoms */
  int nk; /* number of k-points */
  int norb; /* number of orbitals */
  int norbloc;  /* number of orbitals local to this process */

  int *el; /* elements, points to particles_t elements */

  double mu; /* Fermi level */
  double cutoff; /* ??? */

  /*
   * IN FORTRAN:
   * double(DP)  :: H(:, :, :)     ! => NULL()
   * double(DP)  :: S(:, :, :)     ! => NULL()
   * double(DP)  :: rho(:, :, :)   ! => NULL()
   * double(DP)  :: e(:, :, :)     ! => NULL()
   */

  double *H; /* the Hamiltonian (for each k-point) */
  double *S; /* the overlap matrix */
  double *rho; /* The density matrix (rho_ll) */
  double *e; /* H_rl * rho_ll */
       
  /*
   * Band-structure and repulsive energies
   */

  double ebs, erep;

  /*
   * Additional particle information
   *
   * IN FORTRAN:
   * real(DP)              :: n(:)
   * type(notb_element_t)  :: at(:)
   */

  double *n; /* Number of electrons, charges */
  struct notb_element_t *at;
};


/*
 * Allocation, deallocation methods
 * We need to allocate, deallocate in C because Fortran pointer information is
 * lost when passing through c_loc.
 */

extern "C" {

void dense_hamiltonian_allocate(struct dense_hamiltonian_t *self, int nat, 
				int norb);
void dense_hamiltonian_deallocate(struct dense_hamiltonian_t *self);

}

#endif
