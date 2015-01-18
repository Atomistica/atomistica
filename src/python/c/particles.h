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

void f_particles_set_tag(void *, void *);
void f_particles_get_tag(void *, void **);

int f_particles_get_nel(void *);

PyObject *particles_update_elements(particles_t *self, PyObject *args);

#endif
