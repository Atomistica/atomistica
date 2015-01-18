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
#ifndef __NEIGHBORS_H
#define __NEIGHBORS_H

#include "ptrdict.h"

#include "particles.h"


typedef struct {
  PyObject_HEAD

  /* Pointer to Fortran object */
  void *f90obj;

} neighbors_t;


extern PyTypeObject neighbors_type;


int get_neighbors_size(neighbors_t *, particles_t *);
int get_number_of_all_neighbors(neighbors_t *, particles_t *);


void f_neighbors_new(void **);
void f_neighbors_free(void **);

void f_neighbors_init(void *, int);
void f_neighbors_del(void *);
void f_neighbors_set(void *, int, double, double, double);
void f_neighbors_request_interaction_range(void *, double);
void f_neighbors_update(void *, void *, int *);
void f_neighbors_find_neighbor(void *, int, int, int *, int *);

int f_get_neighbors_size(void *);
int f_get_coordination(void *, int, double);
int f_get_coordination_numbers(void *, double, int *);

int f_get_number_of_neighbors(void *, int);
int f_get_number_of_all_neighbors(void *);
void f_get_neighbors(void *, int, int *, double *);
void f_get_seed(void *, int *);
void f_get_all_neighbors(void *, int *, int *, double *);
void f_get_all_neighbors_vec(void *, int *, int *, double *, double *);

void f_pack_per_bond_scalar(void *, double *, double *);
void f_pack_per_bond_3x3(void *, double *, double *);

void f_neighbors_set_tag(void *, void *);
void f_neighbors_get_tag(void *, void **);

#endif
