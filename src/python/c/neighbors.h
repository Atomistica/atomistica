/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
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
void f_get_all_neighbors(void *, int *, int *, double *);
void f_get_all_neighbors_vec(void *, int *, int *, double *, double *);

void f_pack_per_bond_scalar(void *, double *, double *);
void f_pack_per_bond_3x3(void *, double *, double *);

#endif
