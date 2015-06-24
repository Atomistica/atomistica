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

#include "atomisticamodule.h"


/* Python object types:
   Neighbors - Neighbor list
*/


/* Neighbors methods and class
 */

static PyObject *
neighbors_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    neighbors_t *self;

    self = (neighbors_t *)type->tp_alloc(type, 0);
    if (self != NULL) {
      f_neighbors_new(&self->f90obj);
      f_neighbors_set_tag(self->f90obj, self);
    }

    return (PyObject *)self;
}


static void
neighbors_dealloc(neighbors_t *self)
{
#ifdef DEBUG
  printf("[neighbors_dealloc] self = %p\n", self);
#endif
  f_neighbors_free(self->f90obj);
  Py_TYPE(self)->tp_free((PyObject*) self);
}


static int
neighbors_init(neighbors_t *self, PyObject *args, PyObject *kwds)
{
  int avgn;

  if (!PyArg_ParseTuple(args, "i", &avgn))
    return -1; 

  f_neighbors_init(self->f90obj, avgn);

  return 0;
}


static int
_neighbors_update(neighbors_t *self, particles_t *p)
{
  int ierror;

  ierror = ERROR_NONE;
  f_neighbors_update(self->f90obj, p->f90obj, &ierror);
  if (error_to_py(ierror))
    return -1;

  return 0;
}


int
get_neighbors_size(neighbors_t *self, particles_t *p)
{
  int size = f_get_neighbors_size(self->f90obj);
  if (size == 0) {
    _neighbors_update(self, p);
    size = f_get_neighbors_size(self->f90obj);
  }
  return size;
}


int
get_number_of_all_neighbors(neighbors_t *self, particles_t *p)
{
  if (f_get_neighbors_size(self->f90obj) == 0) {
    _neighbors_update(self, p);
  }
  return f_get_number_of_all_neighbors(self->f90obj);
}


static PyObject *
neighbors_request_interaction_range(neighbors_t *self, PyObject *args)
{
  double cutoff;

  if (!PyArg_ParseTuple(args, "d", &cutoff))
    return NULL;

  f_neighbors_request_interaction_range(self->f90obj, cutoff);

  Py_RETURN_NONE;
}


static PyObject *
neighbors_update(neighbors_t *self, PyObject *args)
{
  particles_t *p;

  if (!PyArg_ParseTuple(args, "O!", &particles_type, &p))
    return NULL;

  if (_neighbors_update(self, p))
    return NULL;

  Py_RETURN_NONE;
}


static PyObject *
neighbors_get_neighbors(neighbors_t *self, PyObject *args, PyObject *kwargs)
{
  particles_t *p;
  int i = -1;
  PyObject *vec = NULL, *seed = NULL;

  static char *kwlist[] = { "p", "i", "seed", "vec", NULL };
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!|iO!O!", kwlist,
                                   &particles_type, &p, &i,
                                   &PyBool_Type, &seed,
                                   &PyBool_Type, &vec))
    return NULL;

  if (i >= 0 && vec) {
    PyErr_SetString(PyExc_RuntimeError, "Vector distance not yet supported for per atom neighbors.");
    return NULL;
  }

  if (_neighbors_update(self, p))
    return NULL;

  PyObject *ret;

  if (i >= 0) {
    npy_intp dims[1];
    dims[0] = f_get_number_of_neighbors(self->f90obj, i);

    PyObject *j = PyArray_ZEROS(1, dims, NPY_INT, 1);
    PyObject *r = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
    f_get_neighbors(self->f90obj, i, PyArray_DATA(j), PyArray_DATA(r));

    ret = Py_BuildValue("OO", j, r);

    Py_DECREF(j);
    Py_DECREF(r);
  }
  else {
    npy_intp dims[2];

    PyObject *seed_arr = NULL;
    if (seed && seed == Py_True) {
      dims[0] = data_get_len(p->f90data);
      seed_arr = PyArray_ZEROS(1, dims, NPY_INT, 1);
      f_get_seed(self->f90obj, PyArray_DATA(seed_arr));
    }

    dims[0] = f_get_number_of_all_neighbors(self->f90obj);

    PyObject *i = PyArray_ZEROS(1, dims, NPY_INT, 1);
    PyObject *j = PyArray_ZEROS(1, dims, NPY_INT, 1);
    PyObject *r = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
    if (vec && vec == Py_True) {
      dims[1] = 3;
      PyObject *rvec = PyArray_ZEROS(2, dims, NPY_DOUBLE, 0);

      f_get_all_neighbors_vec(self->f90obj, PyArray_DATA(i), PyArray_DATA(j),
                              PyArray_DATA(rvec), PyArray_DATA(r));

      if (seed_arr) {
        ret = Py_BuildValue("OOOOO", i, j, seed_arr, rvec, r);
      }
      else {
        ret = Py_BuildValue("OOOO", i, j, rvec, r);
      }

      Py_DECREF(rvec);
    }
    else {
      f_get_all_neighbors(self->f90obj, PyArray_DATA(i), PyArray_DATA(j),
                          PyArray_DATA(r));

      if (seed_arr) {
        ret = Py_BuildValue("OOOO", i, j, seed, r);
      }
      else {
        ret = Py_BuildValue("OOO", i, j, r);
      }
    }

    Py_DECREF(i);
    Py_DECREF(j);
    Py_DECREF(r);
    if (seed_arr)
      Py_DECREF(seed_arr);
  }

  return ret;
}


static PyObject *
neighbors_find_neighbor(neighbors_t *self, PyObject *args)
{
  int i, j, n1, n2;
  PyObject *r;
  particles_t *p;

  if (!PyArg_ParseTuple(args, "O!ii", &particles_type, &p, &i, &j))
    return NULL;

  if (_neighbors_update(self, p))
    return NULL;

  /* Indices in Fortran start at 1. */
  i++;
  j++;

  f_neighbors_find_neighbor(self->f90obj, i, j, &n1, &n2);

  r = PyTuple_New(2);
  /* Indices in C and Python start at 0. */
#if PY_MAJOR_VERSION >= 3
  PyTuple_SET_ITEM(r, 0, PyLong_FromLong(n1-1));
  PyTuple_SET_ITEM(r, 1, PyLong_FromLong(n2-1));
#else
  PyTuple_SET_ITEM(r, 0, PyInt_FromLong(n1-1));
  PyTuple_SET_ITEM(r, 1, PyInt_FromLong(n2-1));
#endif
  return r;
}


static PyObject *
neighbors_get_coordination_numbers(neighbors_t *self, PyObject *args)
{
  particles_t *p;
  double cutoff;

  if (!PyArg_ParseTuple(args, "O!d", &particles_type, &p, &cutoff))
    return NULL;

  if (_neighbors_update(self, p))
    return NULL;

  npy_intp dims[1];
  dims[0] = data_get_len(p->f90data);

  PyObject *c = PyArray_ZEROS(1, dims, NPY_INT, 1);
  f_get_coordination_numbers(self->f90obj, cutoff, PyArray_DATA(c));

  return c;
}


static PyObject *
neighbors_coordination(neighbors_t *self, PyObject *args)
{
  int i, c;
  double cutoff;
  particles_t *p;

  if (!PyArg_ParseTuple(args, "O!id", &particles_type, &p, &i, &cutoff))
    return NULL;

  if (_neighbors_update(self, p))
    return NULL;

  i += 1;

  c = f_get_coordination(self->f90obj, i, cutoff);

#if PY_MAJOR_VERSION >= 3
  return PyLong_FromLong(c);
#else
  return PyInt_FromLong(c);
#endif
}


/* Methods declaration */

static PyMethodDef neighbors_methods[] = {
  { "request_interaction_range",
    (PyCFunction) neighbors_request_interaction_range, METH_VARARGS,
    "Set the cutoff radius." },
  { "update", (PyCFunction) neighbors_update, METH_VARARGS,
    "Update the neighbor list." },
  { "get_neighbors", (PyCFunction) neighbors_get_neighbors,
    METH_VARARGS | METH_KEYWORDS,
    "Return complete neighbors list in a three tuple giving, atom i, atom j, "
    "and the distance r." },
  { "find_neighbor", (PyCFunction) neighbors_find_neighbor, METH_VARARGS,
    "Find the neighbor index." },
  { "get_coordination_numbers", (PyCFunction) neighbors_get_coordination_numbers,
    METH_VARARGS, "Coordination count for all atoms." },
  { "coordination", (PyCFunction) neighbors_coordination, METH_VARARGS,
    "Coordination count for a certain atom." },
  { NULL, NULL, 0, NULL }  /* Sentinel */
};


PyTypeObject neighbors_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_atomistica.Neighbors",                    /*tp_name*/
    sizeof(neighbors_t),                        /*tp_basicsize*/
    0,                                          /*tp_itemsize*/
    (destructor)neighbors_dealloc,              /*tp_dealloc*/
    0,                                          /*tp_print*/
    0,                                          /*tp_getattr*/
    0,                                          /*tp_setattr*/
    0,                                          /*tp_compare*/
    0,                                          /*tp_repr*/
    0,                                          /*tp_as_number*/
    0,                                          /*tp_as_sequence*/
    0,                                          /*tp_as_mapping*/
    0,                                          /*tp_hash */
    0,                                          /*tp_call*/
    0,                                          /*tp_str*/
    0,                                          /*tp_getattro*/
    0,                                          /*tp_setattro*/
    0,                                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /*tp_flags*/
    "Neighbors objects",                        /* tp_doc */
    0,		                                /* tp_traverse */
    0,		                                /* tp_clear */
    0,		                                /* tp_richcompare */
    0,		                                /* tp_weaklistoffset */
    0,		                                /* tp_iter */
    0,		                                /* tp_iternext */
    neighbors_methods,                          /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)neighbors_init,                   /* tp_init */
    0,                                          /* tp_alloc */
    neighbors_new,                              /* tp_new */
};

