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

#include "potential.h"

#include "py_f.h"
#include "particles.h"
#include "neighbors.h"


/* Python object types:
   Potential - Single potential
*/

/* Potential methods and class
 */


static PyObject *
potential_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  potential_t *self;
  section_t *zero;
  int i;
  char *name;
  char errstr[200];

#ifdef DEBUG
  printf("[potential_new] %p %p %p\n", type, args, kwds);
#endif

  self = (potential_t *)type->tp_alloc(type, 0);
  if (self != NULL) {

    /* FIXME: the offset *12* assumes the namespace is always _atomistica.* */
    name  = strdup(Py_TYPE(self)->tp_name + 12);

#ifdef DEBUG
    printf("[potential_new] Potential name: %s (%s)\n", name,
           Py_TYPE(self)->tp_name);
#endif

    self->f90class  = NULL;
    for (i = 0; i < N_POTENTIAL_CLASSES; i++) {
      if (!strcmp(name, potential_classes[i].name))
        self->f90class  = &potential_classes[i];
    }
      
    if (!self->f90class) {
      sprintf(errstr, "Internal error: Potential not found: %s\n", name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return NULL;
    }

    zero = NULL;
    self->f90class->new_instance(&self->f90obj, zero, &self->f90members);
#ifdef DEBUG
    printf("[potential_new] pointer = %p\n", self->f90obj);
#endif
  }

  return (PyObject *) self;
}


static void
potential_dealloc(potential_t *self)
{
#ifdef DEBUG
  printf("[potential_dealloc] %p\n", self);
#endif

  self->f90class->free_instance(self->f90obj);

  Py_TYPE(self)->tp_free((PyObject*) self);
}


static int
potential_init(potential_t *self, PyObject *args, PyObject *kwargs)
{
  int ierror = ERROR_NONE;

#ifdef DEBUG
  printf("[potential_init] %p %p %p\n", self, args, kwargs);
#endif

  if (kwargs) {
    if (pydict_to_ptrdict(kwargs, self->f90members))
      return -1;
  }

  self->f90class->init(self->f90obj, &ierror);

  if (error_to_py(ierror))
    return -1;

#ifdef DEBUG
  printf("{potential_init}\n");
#endif

  return 0;
}


/* Attribute set/getters */

static PyObject *
potential_getattro(potential_t *self, PyObject *pyname)
{
  char *name;
  property_t *p;

#if PY_MAJOR_VERSION >= 3
  if (!PyUnicode_Check(pyname)) {
#else
  if (!PyString_Check(pyname)) {
#endif
    PyErr_SetString(PyExc_ValueError, "Key must be a string.");
    return NULL;
  }

#if PY_MAJOR_VERSION >= 3
  PyObject *pybname = PyUnicode_AsASCIIString(pyname);
  name = PyBytes_AS_STRING(pybname);
#else
  name = PyString_AS_STRING(pyname);
#endif

#ifdef DEBUG
  printf("[potential_getattro] %p %p %s\n", self, pyname, name);
#endif

  /* Search potential parameter data */
  p = self->f90members->first_property;

  while (p != NULL && strcmp(p->name, name)) {
#ifdef DEBUG
    printf("[potential_getattro] p->name = %s\n", p->name);
#endif
    p = p->next;
  }

  if (p) {
#if PY_MAJOR_VERSION >= 3
    Py_DECREF(pybname);
#endif
    return property_to_pyobject(p);
}

  /* Not in parameter data? Search in state data */
  if (self->f90class->get_dict) {
#ifdef DEBUG
    printf("[potential_getattro] searching state dictionary\n");
#endif

    int ierror = ERROR_NONE;
    section_t *s = ptrdict_register_section(NULL, self->f90class->name, "");
    PyObject *r = NULL;
    
    self->f90class->get_dict(self->f90obj, s, &ierror);

    if (error_to_py(ierror)) {
#if PY_MAJOR_VERSION >= 3
      Py_DECREF(pybname);
#endif
      return NULL;
    }

    p = s->first_property;
    while (p != NULL && strcmp(p->name, name)) {
#ifdef DEBUG
      printf("[potential_getattro] p->name = %s\n", p->name);
#endif
      p = p->next;
    }

    if (p) r = property_to_pyobject(p);
    ptrdict_cleanup(s);
    if (r) {
#if PY_MAJOR_VERSION >= 3
      Py_DECREF(pybname);
#endif
      return r;
    }
  }

#if PY_MAJOR_VERSION >= 3
  Py_DECREF(pybname);
#endif
  /* Fall back to default getattro */
  return PyObject_GenericGetAttr((PyObject *) self, pyname);
}


/* Methods */

static PyObject *
potential_str(potential_t *self, PyObject *args)
{
#if PY_MAJOR_VERSION >= 3
  return PyUnicode_FromString(self->f90class->name);
#else
  return PyString_FromString(self->f90class->name);
#endif
}


static PyObject *
potential_register_data(potential_t *self, PyObject *args)
{
  particles_t *a;
  int ierror = ERROR_NONE;

  if (!PyArg_ParseTuple(args, "O!", &particles_type, &a))
    return NULL; 

  self->f90class->register_data(self->f90obj, a->f90obj, &ierror);

  if (error_to_py(ierror))
    return NULL;

  Py_RETURN_NONE;
}


static PyObject *
potential_bind_to(potential_t *self, PyObject *args)
{
  particles_t *a;
  neighbors_t *n;
  int ierror = ERROR_NONE;

#ifdef DEBUG
  printf("[potential_bind_to] self = %p\n", self);
#endif

  if (!PyArg_ParseTuple(args, "O!O!", &particles_type, &a, &neighbors_type, &n))
    return NULL; 

  self->f90class->bind_to(self->f90obj, a->f90obj, n->f90obj, &ierror);

  if (error_to_py(ierror))
    return NULL;

  Py_RETURN_NONE;
}


static PyObject *
potential_set_Coulomb(potential_t *self, PyObject *args)
{
  int ierror = ERROR_NONE;
  PyObject *coul;

#ifdef DEBUG
  printf("[potential_set_Coulomb] self = %p\n", self);
#endif

  if (!self->f90class->set_Coulomb) {
    char errstr[1024];
    sprintf(errstr, "Potential %s does not require a Coulomb solver",
            self->f90class->name);
    PyErr_SetString(PyExc_RuntimeError, errstr);
    return NULL;
  }

  if (!PyArg_ParseTuple(args, "O", &coul))
    return NULL; 

  self->f90class->set_Coulomb(self->f90obj, coul, &ierror);

  if (error_to_py(ierror))
    return NULL;

  Py_RETURN_NONE;
}


static PyObject *
potential_get_per_bond_property(potential_t *self, PyObject *args)
{
  int ierror = ERROR_NONE;
  PyObject *coul;

#ifdef DEBUG
  printf("[potential_get_per_bond_property] self = %p\n", self);
#endif

  if (!self->f90class->get_per_bond_property) {
    char errstr[1024];
    sprintf(errstr, "Potential %s does not implement get_per_bond_property.",
            self->f90class->name);
    PyErr_SetString(PyExc_RuntimeError, errstr);
    return NULL;
  }

  particles_t *a;
  neighbors_t *n;
  char *propname;
  if (!PyArg_ParseTuple(args, "O!O!s", &particles_type, &a, &neighbors_type, &n, &propname))
    return NULL;

  /* Get property from potential object */
  npy_intp dims[1];
  dims[0] = get_neighbors_size(n, a);
  PyObject *tmp = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
  self->f90class->get_per_bond_property(self->f90obj, a->f90obj, n->f90obj,
                                        propname, PyArray_DATA(tmp), &ierror);
  if (error_to_py(ierror)) {
    Py_DECREF(tmp);
    return NULL;
  }

  /* Pack neighbor list into Python format */
  dims[0] = get_number_of_all_neighbors(n, a);
  PyObject *prop = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
  f_pack_per_bond_scalar(n->f90obj, PyArray_DATA(tmp), PyArray_DATA(prop));

  /* Release temporary buffer */
  Py_DECREF(tmp);

  return prop;
}


static PyObject *
potential_energy_and_forces(potential_t *self, PyObject *args, PyObject *kwargs)
{
  static char *kwlist[] = {
    "particles",
    "neighbors",
    "forces",
    "mask",
    "charges",
    "epot_per_at",
    "epot_per_bond",
    "f_per_bond",
    "wpot_per_at",
    "wpot_per_bond",
    NULL
  };

  npy_intp dims[3];
  npy_intp strides[3];

  particles_t *a;
  neighbors_t *n;
  PyObject *return_epot_per_at = NULL;
  PyObject *return_epot_per_bond = NULL;
  PyObject *return_f_per_bond = NULL;
  PyObject *return_wpot_per_at = NULL;
  PyObject *return_wpot_per_bond = NULL;

  int ierror = ERROR_NONE;

  double epot;
  PyObject *f = NULL;
  PyObject *wpot;
  PyObject *q = NULL;
  PyObject *mask = NULL;
  PyObject *epot_per_at = NULL;
  PyObject *epot_per_bond = NULL;
  PyObject *f_per_bond = NULL;
  PyObject *wpot_per_at = NULL;
  PyObject *wpot_per_bond = NULL;

  double *epot_per_at_ptr = NULL;
  double *epot_per_bond_ptr = NULL;
  double *f_per_bond_ptr = NULL;
  double *wpot_per_at_ptr = NULL;
  double *wpot_per_bond_ptr = NULL;

  PyObject *r;

  int i, nat;

  /* --- */

#ifdef DEBUG
  printf("[potential_energy_and_forces] self = %p\n", self);
#endif

  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|O!O!O!O!O!O!O!O!",
                                  kwlist, &particles_type, &a, &neighbors_type,
                                  &n, &PyArray_Type, &f, &PyArray_Type, &mask,
                                  &PyArray_Type, &q,
                                  &PyBool_Type, &return_epot_per_at, 
                                  &PyBool_Type, &return_epot_per_bond, 
                                  &PyBool_Type, &return_f_per_bond,
                                  &PyBool_Type, &return_wpot_per_at, 
                                  &PyBool_Type, &return_wpot_per_bond))
    return NULL;

  epot = 0.0;
  nat = data_get_len(a->f90data);

  if (q) {
    if (PyArray_DIM(q, 0) != nat) {
      char errstr[1024];
      sprintf(errstr, "Length of charge array (= %"NPY_INTP_FMT") does not "
              "equal number of atoms (= %i).", PyArray_DIM(q, 0), nat);
      PyErr_SetString(PyExc_RuntimeError, errstr);
      return NULL;
    }
  }

  if (mask) {
    if (mask == Py_None) {
      Py_DECREF(Py_None);
    }
    else {
      mask = PyArray_FROMANY(mask, NPY_INT, 1, 1, 0);
      if (!mask) {
        PyErr_SetString(PyExc_RuntimeError, "Could not convert mask array to "
                        "integer.");
        return NULL;
      }
      if (PyArray_DIM(mask, 0) != nat) {
        char errstr[1024];
        sprintf(errstr, "Length of mask array (= %"NPY_INTP_FMT") does not "
                "equal number of atoms (= %i).", PyArray_DIM(mask, 0), nat);
        PyErr_SetString(PyExc_RuntimeError, errstr);
        Py_DECREF(mask);
        return NULL;
      }
    }
  }

  if (f) {
    if (PyArray_DIM(f, 0) != nat) {
      char errstr[1024];
      sprintf(errstr, "Length of force array (= %"NPY_INTP_FMT") does not equal "
              " number of atoms (= %i).", PyArray_DIM(f, 0), nat);
      PyErr_SetString(PyExc_RuntimeError, errstr);
      if (mask) Py_DECREF(mask);
      return NULL;
    }
    Py_INCREF(f);
  }
  else {
    dims[0] = nat;
    dims[1] = 3;
    strides[0] = dims[1]*NPY_SIZEOF_DOUBLE;
    strides[1] = NPY_SIZEOF_DOUBLE;
    f = (PyObject*) PyArray_New(&PyArray_Type, 2, dims, NPY_DOUBLE, strides,
                                NULL, 0, NPY_FARRAY, NULL);
    memset(PyArray_DATA(f), 0, dims[0]*dims[1]*NPY_SIZEOF_DOUBLE);
  }

  dims[0] = 3;
  dims[1] = 3;
  wpot = PyArray_ZEROS(2, dims, NPY_DOUBLE, 1);

  if (return_epot_per_at) {

    if (return_epot_per_at == Py_True) {
      dims[0] = nat;
      epot_per_at = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
      epot_per_at_ptr = PyArray_DATA(epot_per_at);
    } else {
      epot_per_at  = Py_None;
      Py_INCREF(Py_None);
    }

  }

  if (return_epot_per_bond) {

    if (return_epot_per_bond == Py_True) {
      dims[0] = get_neighbors_size(n, a);
      epot_per_bond = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);
      epot_per_bond_ptr = PyArray_DATA(epot_per_bond);
    } else {
      epot_per_bond = Py_None;
      Py_INCREF(Py_None);
    }

  }

  if (return_f_per_bond) {
    
    if (return_f_per_bond == Py_True) {
      dims[0] = get_neighbors_size(n, a);
      dims[1] = 3;
      strides[0] = dims[1]*NPY_SIZEOF_DOUBLE;
      strides[1] = NPY_SIZEOF_DOUBLE;
      f_per_bond = (PyObject*) PyArray_New(&PyArray_Type, 2, dims,
                                           NPY_DOUBLE, strides, NULL, 0,
                                           NPY_FARRAY, NULL);
      f_per_bond_ptr = PyArray_DATA(f_per_bond);
      memset(f_per_bond_ptr, 0, dims[0]*dims[1]*NPY_SIZEOF_DOUBLE);
    } else {
      f_per_bond  = Py_None;
      Py_INCREF(Py_None);
    }

  }

  if (return_wpot_per_at) {

    if (return_wpot_per_at == Py_True) {
      dims[0]            = nat;
      dims[1]            = 3;
      dims[2]            = 3;
      strides[0]         = dims[1]*dims[2]*NPY_SIZEOF_DOUBLE;
      strides[1]         = dims[2]*NPY_SIZEOF_DOUBLE;
      strides[2]         = NPY_SIZEOF_DOUBLE;
      wpot_per_at      = (PyObject*) PyArray_New(&PyArray_Type, 3, dims,
                         NPY_DOUBLE, strides, NULL, 0,
                         NPY_FARRAY, NULL);
      wpot_per_at_ptr  = PyArray_DATA(wpot_per_at);
      memset(wpot_per_at_ptr, 0, dims[0]*dims[1]*dims[2]*NPY_SIZEOF_DOUBLE);
    } else {
      wpot_per_at  = Py_None;
      Py_INCREF(Py_None);
    }

  }

  if (return_wpot_per_bond) {

    if (return_wpot_per_bond == Py_True) {
      dims[0]            = get_neighbors_size(n, a);
      dims[1]            = 3;
      dims[2]            = 3;
      strides[0]         = dims[1]*dims[2]*NPY_SIZEOF_DOUBLE;
      strides[1]         = dims[2]*NPY_SIZEOF_DOUBLE;
      strides[2]         = NPY_SIZEOF_DOUBLE;
      wpot_per_bond      = (PyObject*) PyArray_New(&PyArray_Type, 3, dims,
                           NPY_DOUBLE, strides, NULL, 
                           0, NPY_FARRAY, NULL);
      wpot_per_bond_ptr  = PyArray_DATA(wpot_per_bond);
      memset(wpot_per_bond_ptr, 0, dims[0]*dims[1]*dims[2]*NPY_SIZEOF_DOUBLE);
    } else {
      wpot_per_bond  = Py_None;
      Py_INCREF(Py_None);
    }

  }

#ifdef DEBUG
  printf("[potential_energy_and_forces] self->f90class->name = %s\n",
     self->f90class->name);
  printf("[potential_energy_and_forces] self->f90obj = %p\n",
     self->f90obj);
  printf("[potential_energy_and_forces] a->f90obj = %p\n",
     a->f90obj);
  printf("[potential_energy_and_forces] n->f90obj = %p\n",
     n->f90obj);
  printf("[potential_energy_and_forces] self->f90class->energy_and_forces = %p\n",
     self->f90class->energy_and_forces);
#endif

  double *mask_data = NULL;
  if (mask)  mask_data = PyArray_DATA(mask);
  double *q_data = NULL;
  if (q)  q_data = PyArray_DATA(q);
  self->f90class->energy_and_forces(self->f90obj, a->f90obj, n->f90obj,
                                    q_data, &epot, PyArray_DATA(f),
                                    PyArray_DATA(wpot), mask_data, epot_per_at_ptr,
                                    epot_per_bond_ptr, f_per_bond_ptr,
                                    wpot_per_at_ptr, wpot_per_bond_ptr,
                                    &ierror);

  /*
   * Now we need to reorder the per-bond properties such that some Python
   * script can actually make sense out of the data.
   */

  if (epot_per_bond_ptr) {
    dims[0] = get_number_of_all_neighbors(n, a);
    PyObject *tmp = PyArray_ZEROS(1, dims, NPY_DOUBLE, 1);

    f_pack_per_bond_scalar(n->f90obj, epot_per_bond_ptr, PyArray_DATA(tmp));

    Py_DECREF(epot_per_bond);
    epot_per_bond = tmp;
  }

  if (wpot_per_bond_ptr) {
    dims[0] = get_number_of_all_neighbors(n, a);
    dims[1] = 3;
    dims[2] = 3;
    PyObject *tmp = PyArray_ZEROS(3, dims, NPY_DOUBLE, 0);

    f_pack_per_bond_3x3(n->f90obj, wpot_per_bond_ptr, PyArray_DATA(tmp));

    Py_DECREF(wpot_per_bond);
    wpot_per_bond = tmp;
  }

#ifdef DEBUG
  printf("[potential_energy_and_forces] epot = %f\n", epot);
#endif

  if (error_to_py(ierror))
    return NULL;

  /* --- Compose return tuple --- */

  i = 3;
  if (epot_per_at)    i++;
  if (epot_per_bond)  i++;
  if (f_per_bond)     i++;
  if (wpot_per_at)    i++;
  if (wpot_per_bond)  i++;

  r  = PyTuple_New(i);
  if (!r)
    return NULL;

  PyTuple_SET_ITEM(r, 0, PyFloat_FromDouble(epot));
  PyTuple_SET_ITEM(r, 1, f);
  PyTuple_SET_ITEM(r, 2, wpot);

  i = 2;
  if (epot_per_at) {
    i++;
    PyTuple_SET_ITEM(r, i, epot_per_at);
  }
  if (epot_per_bond) {
    i++;
    PyTuple_SET_ITEM(r, i, epot_per_bond);
  }
  if (f_per_bond) {
    i++;
    PyTuple_SET_ITEM(r, i, f_per_bond);
  }
  if (wpot_per_at) {
    i++;
    PyTuple_SET_ITEM(r, i, wpot_per_at);
  }
  if (wpot_per_bond) {
    i++;
    PyTuple_SET_ITEM(r, i, wpot_per_bond);
  }

#ifdef DEBUG
  printf("{potential_energy_and_forces}\n");
#endif

  Py_XDECREF(mask);

  return r;
}


/* Methods declaration */

static PyMethodDef potential_methods[] = {
  { "register_data", (PyCFunction) potential_register_data, METH_VARARGS,
    "Register internal data fields with a particles object." },
  { "bind_to", (PyCFunction) potential_bind_to, METH_VARARGS,
    "Bind this potential to a certain Particles and Neighbors object. This is "
    "to be called if either one changes." },
  { "set_Coulomb", (PyCFunction) potential_set_Coulomb, METH_VARARGS,
    "Set the object that handles Coulomb callbacks." },
  { "get_per_bond_property", (PyCFunction) potential_get_per_bond_property,
    METH_VARARGS, "Return a named property that is defined per bond." },
  { "energy_and_forces", (PyCFunctionWithKeywords) potential_energy_and_forces,
    METH_VARARGS | METH_KEYWORDS,
    "Compute the forces and return the potential energy." },
  { NULL, NULL, 0, NULL }  /* Sentinel */
};


PyTypeObject potential_type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_atomistica.Potential",                    /*tp_name*/
    sizeof(potential_t),                        /*tp_basicsize*/
    0,                                          /*tp_itemsize*/
    (destructor)potential_dealloc,              /*tp_dealloc*/
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
    (reprfunc)potential_str,                    /*tp_str*/
    (getattrofunc)potential_getattro,           /*tp_getattro*/
    0,                                          /*tp_setattro*/
    0,                                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /*tp_flags*/
    "Potential objects",                        /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    0,                                          /* tp_iter */
    0,                                          /* tp_iternext */
    potential_methods,                          /* tp_methods */
    0,                                          /* tp_members */
    0,                                          /* tp_getset */
    0,                                          /* tp_base */
    0,                                          /* tp_dict */
    0,                                          /* tp_descr_get */
    0,                                          /* tp_descr_set */
    0,                                          /* tp_dictoffset */
    (initproc)potential_init,                   /* tp_init */
    0,                                          /* tp_alloc */
    potential_new,                              /* tp_new */
};
