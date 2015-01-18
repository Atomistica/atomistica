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

#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL ATOMISTICA_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

#include "error.h"

#include "atomisticamodule.h"


/*
 * Callback methods: These are called from Fortran to obtain the electrostatic
 * potential, field etc.
 */

void py_coulomb_set_hubbard_u(PyObject *self, void *p, double *U, int *error)
{
  particles_t *py_p;
  PyObject *py_U, *r;
  npy_intp dims[1];

  INIT_ERROR(error);

#ifdef DEBUG
  printf("[py_coulomb_set_Hubbard_U] %s %p %p\n",
	 PyString_AsString(PyObject_Str(self)), U, error);
#endif

  f_particles_get_tag(p, (void**) &py_p);
  assert(py_p->f90obj == p);

  dims[0] = f_particles_get_nel(py_p->f90obj);
  py_U = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, U);

  r = PyObject_CallMethod(self, "set_Hubbard_U", "(OO)", py_p, py_U);

  Py_DECREF(py_U);
  PASS_PYTHON_ERROR(error, r);
  Py_DECREF(r);
}


void py_coulomb_potential_and_field(PyObject *self, void *p, void *nl,
				    double *q, double *phi, double *epot,
				    double *E, double *wpot, int *error)
{
  particles_t *py_p;
  neighbors_t *py_nl;
  PyObject *py_q, *py_phi, *py_epot, *py_E, *py_wpot, *r;

  int nat;
  npy_intp dims[2];

#ifdef DEBUG
  printf("[py_coulomb_potential_and_field]\n");
#endif

  f_particles_get_tag(p, (void**) &py_p);
  assert(py_p->f90obj == p);
  nat = data_get_len(py_p->f90data);

  f_neighbors_get_tag(nl, (void**) &py_nl);
  assert(py_nl->f90obj == nl);

  dims[0] = nat;
  dims[1] = 3;

  py_q = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, q);
  py_phi = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, phi);
  py_E = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, E);

  dims[0] = 1;
  py_epot = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, epot);

  dims[0] = 3;
  dims[1] = 3;
  py_wpot = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, wpot);

  r = PyObject_CallMethod(self, "potential_and_field", "(OOOOOOO)", py_p,
                          py_nl, py_q, py_phi, py_epot, py_E, py_wpot);

  Py_DECREF(py_q);
  Py_DECREF(py_phi);
  Py_DECREF(py_E);
  Py_DECREF(py_epot);
  Py_DECREF(py_wpot);
  PASS_PYTHON_ERROR(error, r);
  Py_DECREF(r);
}


void py_coulomb_potential(PyObject *self, void *p, void *nl, double *q,
			  double *phi, int *error)
{
  particles_t *py_p;
  neighbors_t *py_nl;
  PyObject *py_q, *py_phi, *r;

  int nat;
  npy_intp dims[2];
  
#ifdef DEBUG
  printf("[py_coulomb_potential] self = %s\n",
	 PyString_AsString(PyObject_Str(self)));
#endif

  f_particles_get_tag(p, (void**) &py_p);
  assert(py_p->f90obj == p);
  nat = data_get_len(py_p->f90data);

  f_neighbors_get_tag(nl, (void**) &py_nl);
  assert(py_nl->f90obj == nl);

  dims[0] = nat;
  dims[1] = 3;

  py_q = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, q);
  py_phi = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, phi);

  r = PyObject_CallMethod(self, "potential", "(OOOO)", py_p, py_nl, py_q,
			  py_phi);

  Py_DECREF(py_q);
  Py_DECREF(py_phi);
  PASS_PYTHON_ERROR(error, r);
  Py_DECREF(r);
}
