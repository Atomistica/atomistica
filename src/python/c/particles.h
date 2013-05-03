/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
#ifndef __PARTICLES_H
#define __PARTICLES_H

#include <Python.h>

typedef struct {
  PyObject_HEAD

  BOOL initialized;

  /* Pointer to Fortran object */
  void *f90obj;

  /* Pointer to the associated data object */
  void *f90data;

} particles_t;


extern PyTypeObject particles_type;


/* particles_t */
void f_particles_new(void **);
void f_particles_free(void *);

void f_particles_init(void *);
void f_particles_allocate(void *, int, int *);
void f_particles_del(void *);
void f_particles_update_elements(void *);

void f_particles_set_cell(void *, double *, BOOL *, int *);
void f_particles_set_lees_edwards(void *, double *, double *, int *);

void f_particles_inbox(void *);

void f_particles_i_changed_positions(void *);

void f_particles_get_data(void *, void **);

PyObject *particles_update_elements(particles_t *self, PyObject *args);

#endif
