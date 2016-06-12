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

#include <stdlib.h>

#include "dense_hamiltonian.h"

#include "materials.h"


extern "C" void
dense_hamiltonian_allocate(struct dense_hamiltonian_t *self, int nat, int norb)
{
  int nk = self->nk;

  if (self->nat < nat || self->norb < norb) {
    dense_hamiltonian_deallocate(self);
  }

  self->nat = nat;
  self->norb = norb;

  if (!self->H)
    self->H = (double *) malloc(norb*norb*nk*sizeof(double));
  if (!self->S)
    self->S = (double *) malloc(norb*norb*nk*sizeof(double));
  if (!self->rho)
    self->rho = (double *) malloc(norb*norb*nk*sizeof(double));
  if (!self->e)
    self->e = (double *) malloc(norb*norb*nk*sizeof(double));

  if (!self->n)
    self->n = (double *) malloc(nat*sizeof(double));
  if (!self->at)
    self->at = (struct notb_element_t *)
      malloc(nat*sizeof(struct notb_element_t));
}


extern "C" void
dense_hamiltonian_deallocate(struct dense_hamiltonian_t *self)
{
  if (self->H)  free(self->H);
  if (self->S)  free(self->S);
  if (self->rho)  free(self->rho);
  if (self->e)  free(self->e);

  if (self->n)  free(self->n);
  if (self->at)  free(self->at);

  self->H = NULL;
  self->S = NULL;
  self->rho = NULL;
  self->e = NULL;
  self->n = NULL;
  self->at = NULL;
}

