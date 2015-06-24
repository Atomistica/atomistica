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

/* Backport type definitions from Python 2.5's object.h */ 
#if PY_VERSION_HEX < 0x02050000 
typedef int Py_ssize_t; 
typedef Py_ssize_t (*lenfunc)(PyObject *); 
typedef PyObject *(*ssizeargfunc)(PyObject *, Py_ssize_t); 
typedef PyObject *(*ssizessizeargfunc)(PyObject *, Py_ssize_t, Py_ssize_t); 
typedef int(*ssizeobjargproc)(PyObject *, Py_ssize_t, PyObject *); 
typedef int(*ssizessizeobjargproc)(PyObject *, Py_ssize_t, Py_ssize_t, PyObject *); 
#endif /* PY_VERSION_HEX */ 

/* Python object types:
   Particles - List of particles, including velocities and forces
*/


/* Helper methods
 */

static void
particles_get_ptr(particles_t *self)
{
  f_particles_get_data(self->f90obj, &self->f90data);
#ifdef DEBUG
  printf("[particles_get_ptr] self->f90obj = %p, self->f90data = %p\n",
   self->f90obj, self->f90data);
#endif
}


PyObject *
particles_update_elements(particles_t *self, PyObject *args)
{
  f_particles_update_elements(self->f90obj);

  particles_get_ptr(self);

  Py_RETURN_NONE;
}


/* Particles methods and class
 */

static PyObject *
particles_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  particles_t *self;
  
  self = (particles_t *)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->initialized = 0;
    f_particles_new(&self->f90obj);
    f_particles_set_tag(self->f90obj, self);
  }

  return (PyObject *)self;
}


static void
particles_dealloc(particles_t *self)
{
#ifdef DEBUG
  printf("[particles_dealloc] self = %p\n", self);
#endif
  if (self->initialized)
    f_particles_del(self->f90obj);
  f_particles_free(self->f90obj);
  Py_TYPE(self)->tp_free((PyObject*) self);
}


static int
particles_init(particles_t *self, PyObject *args)
{
#ifdef DEBUG
  printf("[particles_init] %p %p\n", self, args);
#endif

  f_particles_init(self->f90obj);

  return 0;
}


/* Methods */

static PyObject*
particles_allocate(particles_t *self, PyObject *args)
{
  int nat;
  int ierror = ERROR_NONE;

  if (!PyArg_ParseTuple(args, "i", &nat))
    return NULL;

  f_particles_allocate(self->f90obj, nat, &ierror);
  if (error_to_py(ierror))
    return NULL;

  self->initialized = TRUE;

  particles_get_ptr(self);

  Py_RETURN_NONE;
}


static PyObject*
particles_set_lees_edwards(particles_t *self, PyObject *value)
{
  PyObject *dx_arr, *dv_arr = NULL;
  double *dx, dv[3];

  int ierror = ERROR_NONE;

  if (!PyArg_ParseTuple(value, "O|O", &dx_arr, &dv_arr))
    return NULL;
    
  if (dv_arr == Py_None)
    dv_arr = NULL;
    
  dx_arr = PyArray_FROMANY(dx_arr, NPY_DOUBLE, 1, 1, 0);
    
  if (dv_arr) {
    dv_arr = PyArray_FROMANY(dv_arr, NPY_DOUBLE, 1, 1, 0);
    if (!dv_arr) {
      Py_DECREF(dx_arr);
      return NULL;
    }
  }

  if (PyArray_DIM(dx_arr, 0) != 3) {
    PyErr_SetString(PyExc_TypeError, "dx needs to be 3 vector.");
    Py_DECREF(dx_arr);
    if (dv_arr) {
      Py_DECREF(dv_arr);
    }
    return NULL;
  }
  dx = (double *)  PyArray_DATA(dx_arr);

  dv[0] = 0.0;
  dv[1] = 0.0;
  dv[2] = 0.0;
  if (dv_arr) {
    if (PyArray_DIM(dv_arr, 0) != 3) {
      PyErr_SetString(PyExc_TypeError, "dv needs to be 3 vector.");
      Py_DECREF(dx_arr);
      Py_DECREF(dv_arr);
      return NULL;
    }
    dv[0] = ((double *)  PyArray_DATA(dv_arr))[0];
    dv[1] = ((double *)  PyArray_DATA(dv_arr))[1];
    dv[2] = ((double *)  PyArray_DATA(dv_arr))[2];
  }
  
  f_particles_set_lees_edwards(self->f90obj, dx, dv, &ierror);
  if (error_to_py(ierror))
      return NULL;

  Py_DECREF(dx_arr);
  if (dv_arr) {
    Py_DECREF(dv_arr);
  }

  Py_RETURN_NONE;
}


static PyObject *
data_array_by_name(particles_t *self, char *key)
{
  int data_type;
  BOOL ex;
  int ierror = ERROR_NONE;

  void *array;

  PyObject *r;

  npy_intp dims[3];
#ifndef SEP_XYZ
  npy_intp strides[3];
#endif

  char errstr[100];

#ifdef DEBUG
  printf("[data_array_by_name] self = %p, key = %p\n", self, key);
  printf("[data_array_by_name] self->f90obj = %p, self->f90data = %p\n",
   self->f90obj, self->f90data);
  printf("[data_array_by_name] key = %s\n", key);
#endif

  ex = f_data_exists(self->f90data, key, &data_type);

#ifdef DEBUG
  printf("[data_array_by_name] ex = %i\n", ex);
#endif

  r  = NULL;

  if (ex) {

    switch (data_type) {

    case TYPE_REAL:
      real_ptr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0] = data_get_len(self->f90data);
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_REAL, dim = %i\n", dims[0]);
#endif
      r  = PyArray_New(&PyArray_Type, 1, dims, NPY_DOUBLE, NULL, array, 0,
           NPY_FARRAY, NULL);
      break;

    case TYPE_INTEGER:
      integer_ptr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0] = data_get_len(self->f90data);
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_INTEGER, dim = %i\n", dims[0]);
#endif
      r  = PyArray_New(&PyArray_Type, 1, dims, NPY_INT, NULL, array, 0,
           NPY_FARRAY, NULL);
      break;

    case TYPE_REAL3:
      realx_ptr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0]     = data_get_len(self->f90data);
      dims[1]     = 3;
      strides[0]  = 3*NPY_SIZEOF_DOUBLE;
      strides[1]  = NPY_SIZEOF_DOUBLE;
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_REAL3, dim = %i %i, strides = %i %i\n",
       dims[0], dims[1], strides[0], strides[1]);
#endif
      r  = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides, array, 0,
           NPY_BEHAVED, NULL);
      break;

    case TYPE_REAL3x3:
#ifdef DEBUG
      printf("[data_array_by_name] Type is REAL3x3\n");
#endif
      realxxx_ptr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0] = data_get_len(self->f90data);
      dims[1] = 3;
      dims[2] = 3;
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_REAL3x3, dim = %i %i %i\n", dims[0],
       dims[1], dims[2]);
#endif
      r  = PyArray_New(&PyArray_Type, 3, dims, NPY_DOUBLE, NULL, array, 0,
           NPY_FARRAY, NULL);
      break;

    case TYPE_REAL_ATTR:
      real_attr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

#ifdef DEBUG
      printf("[data_array_by_name] TYPE_REAL_ATTR\n");
#endif
      r  = PyFloat_FromDouble(*((double*) array));
      break;

    case TYPE_REAL3_ATTR:
      real3_attr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0] = 3;
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_REAL3_ATTR, dim = %i\n", dims[0]);
#endif
      r  = PyArray_New(&PyArray_Type, 1, dims, NPY_DOUBLE, NULL, array, 0,
           NPY_FARRAY, NULL);
      break;

    case TYPE_REAL3x3_ATTR:
      real3x3_attr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0] = 3;
      dims[1] = 3;
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_REAL3_ATTR, dim = %i %i\n", dims[0],
       dims[1]);
#endif
      r  = PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, NULL, array, 0,
           NPY_FARRAY, NULL);
      break;

    case TYPE_INTEGER_ATTR:
      integer_attr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

#ifdef DEBUG
      printf("[data_array_by_name] TYPE_INTEGER_ATTR\n");
#endif
#if PY_MAJOR_VERSION >= 3
      r  = PyLong_FromLong(*((int*) array));
#else
      r  = PyInt_FromLong(*((int*) array));
#endif
      break;

    case TYPE_INTEGER3_ATTR:
      integer3_attr_by_name(self->f90data, key, &array, &ierror);
      if (error_to_py(ierror))
        return NULL;

      dims[0] = 3;
#ifdef DEBUG
      printf("[data_array_by_name] TYPE_INTEGER3_ATTR, dim = %i\n", dims[0]);
#endif
      r  = PyArray_New(&PyArray_Type, 1, dims, NPY_INT, NULL, array, 0,
                       NPY_FARRAY, NULL);
      break;

    default:
      sprintf(errstr, "InternalError: Unknown type returned for field or "
              "attribute '%s'.", key);
      PyErr_SetString(PyExc_KeyError, errstr);
      r  = NULL;

    }

  }

  return r;
}


static Py_ssize_t
particles_len(particles_t *self)
{
  return data_get_len(self->f90data);
}


static PyObject *
particles_getitem(particles_t *self, PyObject *key)
{
  PyObject *r;
  char errstr[100];

#if PY_MAJOR_VERSION >= 3
  if (!PyUnicode_Check(key)) {
#else
  if (!PyString_Check(key)) {
#endif
    PyErr_SetString(PyExc_ValueError, "Key must be a string.");
    return NULL;
  }

#ifdef DEBUG
  printf("[particles_getitem] key = %s\n", PyString_AS_STRING(key));
#endif

#if PY_MAJOR_VERSION >= 3
  PyObject *bkey = PyUnicode_AsASCIIString(key);
  r = (PyObject*) data_array_by_name(self, PyBytes_AS_STRING(bkey));
#else
  r = (PyObject*) data_array_by_name(self, PyString_AS_STRING(key));
#endif

  if (!r) {
#if PY_MAJOR_VERSION >= 3
    sprintf(errstr, "No field or attribute '%s' defined for this object.",
            PyBytes_AS_STRING(bkey));
#else
    sprintf(errstr, "No field or attribute '%s' defined for this object.",
            PyString_AS_STRING(key));
#endif
    PyErr_SetString(PyExc_KeyError, errstr);
  };

#if PY_MAJOR_VERSION >= 3
  Py_DECREF(bkey);
#endif
  return r;
}


static PyObject *
particles_getattro(particles_t *self, PyObject *key)
{
  PyObject *r;

#ifdef DEBUG
  printf("[particles_getattro] %p %p\n", self, key);
  printf("[particles_getattro] key = %s\n", PyString_AS_STRING(key));
#endif
  
  r = NULL;
  if (self->initialized) {
#if PY_MAJOR_VERSION >= 3
    PyObject *bkey = PyUnicode_AsASCIIString(key);
    r = (PyObject*) data_array_by_name(self, PyBytes_AS_STRING(bkey));
    Py_DECREF(bkey);
#else
    r = (PyObject*) data_array_by_name(self, PyString_AS_STRING(key));
#endif
  }

#ifdef DEBUG
  printf("[particles_getattro] r = %p\n", r);
#endif

  if (!r) {
    r = PyObject_GenericGetAttr((PyObject*) self, key);
  }

  return r;
}


static PyObject *
particles_inbox(particles_t *self, PyObject *args)
{
  f_particles_inbox(self->f90obj);

  Py_RETURN_NONE;
}


static PyObject *
particles_I_changed_positions(particles_t *self, PyObject *args)
{
  f_particles_i_changed_positions(self->f90obj);

  Py_RETURN_NONE;
}


/* Get-seters */

static PyObject*
particles_set_cell(particles_t *self, PyObject *args)
{
  PyObject *Abox_obj, *pbc_obj = NULL;
  PyArrayObject *Abox_arr;
  PyArrayObject *pbc_arr = NULL;
  double *Abox;
  npy_bool *pbc;
  BOOL pbc_for[3];
  int ierror = ERROR_NONE;

  if (!PyArg_ParseTuple(args, "O|O", &Abox_obj, &pbc_obj))
    return NULL;

  Abox_arr = (PyArrayObject *) PyArray_FROMANY(Abox_obj, NPY_DOUBLE, 2, 2,
                                               NPY_C_CONTIGUOUS);
  if (!Abox_arr)
    return NULL;
  Abox = DOUBLEP(Abox_arr);

  pbc_for[0] = 1;
  pbc_for[1] = 1;
  pbc_for[2] = 1;
  if (pbc_obj) {
    pbc_arr = (PyArrayObject *) PyArray_FROMANY(pbc_obj, NPY_BOOL, 1, 1,
                                                NPY_C_CONTIGUOUS);
    if (!pbc_arr)
      return NULL;
    pbc = (npy_bool *) BOOLP(pbc_arr);

    pbc_for[0] = pbc[0];
    pbc_for[1] = pbc[1];
    pbc_for[2] = pbc[2];
  }

#ifdef DEBUG
  printf("[particles_set_cell] pbc_for %i, %i, %i\n", pbc_for[0], pbc_for[1],
   pbc_for[2]);
#endif
  f_particles_set_cell(self->f90obj, Abox, pbc_for, &ierror);
  if (error_to_py(ierror))
    return NULL;

  Py_RETURN_NONE;
}


/* Methods declaration */

static PyMethodDef particles_methods[] = {
  { "allocate",
    (PyCFunction) particles_allocate, METH_VARARGS,
    "Allocate the particles object to a certain length." },
  { "inbox",
    (PyCFunction) particles_inbox, METH_NOARGS,
    "Wrap all particles into the box." },
  { "set_cell",
    (PyCFunction) particles_set_cell, METH_VARARGS,
    "Set the simulation cell." },
  { "set_lees_edwards",
    (PyCFunction) particles_set_lees_edwards, METH_VARARGS,
    "Enable Lees-Edwards boundary conditions." },
  { "update_elements",
    (PyCFunction) particles_update_elements, METH_NOARGS,
    "Update internal list of elements." },
  { "I_changed_positions",
    (PyCFunction) particles_I_changed_positions, METH_NOARGS,
    "Notifiy the Particles object that the positions have been changed." },
  { NULL, NULL, 0, NULL }  /* Sentinel */
};


static PyMappingMethods particles_as_mapping = {
  (lenfunc)particles_len,          /*mp_length*/
  (binaryfunc)particles_getitem,   /*mp_subscript*/
  NULL                             /*mp_ass_subscript*/
};



PyTypeObject particles_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_atomistica.Particles",                    /*tp_name*/
    sizeof(particles_t),                        /*tp_basicsize*/
    0,                                          /*tp_itemsize*/
    (destructor)particles_dealloc,              /*tp_dealloc*/
    0,                                          /*tp_print*/
    0,                                          /*tp_getattr*/
    0,                                          /*tp_setattr*/
    0,                                          /*tp_compare*/
    0,                                          /*tp_repr*/
    0,                                          /*tp_as_number*/
    0,                                          /*tp_as_sequence*/
    &particles_as_mapping,                      /*tp_as_mapping*/
    0,                                          /*tp_hash */
    0,                                          /*tp_call*/
    0,                                          /*tp_str*/
    (getattrofunc)particles_getattro,           /*tp_getattro*/
    0,                                          /*tp_setattro*/
    0,                                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /*tp_flags*/
    "Particles objects",                        /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    particles_methods,                    /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                        /* tp_dict */
    0,                                      /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    (initproc)particles_init,             /* tp_init */
    0,                                          /* tp_alloc */
    particles_new,                        /* tp_new */
};

