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

#include "py_f.h"

#include "error.h"

#define MAX_STR 100

#define min(x, y)  ( (x) < (y) ? x : y )

int
error_to_py(int ierror)
{
  char errstr[ERRSTRLEN];

  if (ierror != ERROR_NONE) {
    get_full_error_string(errstr);
    PyErr_SetString(PyExc_RuntimeError, errstr);
    return 1;
  } else {
    return 0;
  }   
}


void
py_to_error(char *file, int line, int *ierror)
{
  PyObject *ptype, *pvalue, *ptraceback;
  PyErr_Fetch(&ptype, &pvalue, &ptraceback); 

#if PY_MAJOR_VERSION >= 3
  PyObject *bpvalue = PyUnicode_AsASCIIString(pvalue);
  c_push_error_with_info(PyBytes_AS_STRING(bpvalue), file, line,
                         ERROR_UNSPECIFIED);
  Py_DECREF(bpvalue);
#else
  c_push_error_with_info(PyString_AS_STRING(pvalue), file, line,
                         ERROR_UNSPECIFIED);
#endif

  PyErr_Clear();

  if (ierror != NULL) {
    *ierror = ERROR_UNSPECIFIED;
  }
  else {
    c_error_abort(ERROR_UNSPECIFIED);
  }
}


void
pystring_to_fstring(PyObject *pystr, char *forstr, int len)
{
  char *str;
  int j;

#if PY_MAJOR_VERSION >= 3
  PyObject *pybstr = PyUnicode_AsASCIIString(pystr);
  str = PyBytes_AS_STRING(pybstr);
#else
  str = PyString_AS_STRING(pystr);
#endif
  strncpy(forstr, str, len);
  j = strlen(str);
  if (j < len)  memset(forstr+j, ' ', len-j);
#if PY_MAJOR_VERSION >= 3
  Py_DECREF(pybstr);
#endif
}


void
cstring_to_fstring(char *cstr, int clen, char *forstr, int forlen)
{
  int j;

  strncpy(forstr, cstr, min(clen, forlen));
  j = min(strlen(cstr), clen);
  if (j < forlen)  memset(forstr+j, ' ', forlen-j);
}


PyObject*
fstring_to_pystring(char *forstr, int len)
{
  char str[MAX_STR];
  int j;

  strncpy(str, forstr, min(MAX_STR, len));
  j = min(len, strlen(str)-1);
  while (j > 0 && str[j] == ' ')  j--;
  str[j+1] = 0;

#if PY_MAJOR_VERSION >= 3
  return PyUnicode_FromString(str);
#else
  return PyString_FromString(str);
#endif
}


int
pyobject_to_property(PyObject *value, property_t *p)
{
  BOOL b;
  double d;
  int i, j, k;
  char *str, *str2;

  char errstr[1024];

  PyObject *t, *bvalue;
  PyArrayObject *arr;

#ifdef DEBUG
  printf("[pyobject_to_property] Property: %s, kind %i\n", p->name, p->kind);
#endif

  if (!p->ptr) {
    sprintf(errstr, "Pointer for property %s is NULL.", p->name);
    PyErr_SetString(PyExc_RuntimeError, errstr);
    return -1;
  }

  switch (p->kind) {
  case PK_INT:
#if PY_MAJOR_VERSION >= 3
    i = PyLong_AsLong(value);
#else
    i = PyInt_AsLong(value);
#endif
    if (i == -1 && PyErr_Occurred())
      return -1;
    *((int*) p->ptr) = i;
    break;
  case PK_DOUBLE:
    d = PyFloat_AsDouble(value);
    if (PyErr_Occurred())
      return -1;
    *((double*) p->ptr) = d;
    break;
  case PK_BOOL:
    if (!PyBool_Check(value))
      return -1;
    b = value == Py_True;
    *((BOOL*) p->ptr) = b;
    break;
  case PK_STRING:
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(value)) {
#else
    if (!PyString_Check(value)) {
#endif
      sprintf(errstr,
              "Property '%s' of section '%s' should be a string.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    }
#if PY_MAJOR_VERSION >= 3
    bvalue = PyUnicode_AsASCIIString(value);
    str = PyBytes_AS_STRING(bvalue);
#else
    str = PyString_AS_STRING(value);
#endif
    strncpy((char*) p->ptr, str, p->tag-1);
#if PY_MAJOR_VERSION >= 3
    Py_DECREF(bvalue);
#endif
    break;
  case PK_FORTRAN_STRING:
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(value)) {
#else
    if (!PyString_Check(value)) {
#endif
      sprintf(errstr,
              "Property '%s' of section '%s' should be a string.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    }
    pystring_to_fstring(value, (char*) p->ptr, p->tag);
    break;
  case PK_POINT:
    if (!PyTuple_Check(value)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a tuple.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    }
    if (PyTuple_Size(value) != 3) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 3-tuple of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    }
    for (i = 0; i < 3; i++) {
      t = PyTuple_GET_ITEM(value, i);
      if (!PyFloat_Check(t)) {
        sprintf(errstr,
                "Property '%s' of section '%s' should be a 3-tuple of "
                "floats.\n",
                p->name, p->parent->name);
        PyErr_SetString(PyExc_ValueError, errstr);
        return -1;
      }
      ((double *) p->ptr)[i] = PyFloat_AS_DOUBLE(t);
    }
    break;
  case PK_INTPOINT:
    if (!PyTuple_Check(value)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 3-tuple of "
              "integers.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    if (PyTuple_Size(value) != 3) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 3-tuple of "
              "integers.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    }
    for (i = 0; i < 3; i++) {
      t = PyTuple_GET_ITEM(value, i);
#if PY_MAJOR_VERSION >= 3
      if (!PyLong_Check(t)) {
#else
      if (!PyInt_Check(t)) {
#endif
        sprintf(errstr,
                "Property '%s' of section '%s' should be a 3-tuple of "
                "integers.\n",
                p->name, p->parent->name);
        PyErr_SetString(PyExc_ValueError, errstr);
        return -1;
      }
#if PY_MAJOR_VERSION >= 3
      ((int *) p->ptr)[i] = PyLong_AS_LONG(t);
#else
      ((int *) p->ptr)[i] = PyInt_AS_LONG(t);
#endif
    }
    
    break;
  case PK_LIST:
    if (PyFloat_Check(value)) {
      *p->tag5 = 1;
      *((double*) p->ptr) = PyFloat_AS_DOUBLE(value);
    } else {
      PyArray_Converter(value, (PyObject**) &arr);
      if (!PyArray_ISFLOAT(arr)) {
        sprintf(errstr,
                "Property '%s' of section '%s' should be a list of "
                "floats.\n",
                p->name, p->parent->name);
        PyErr_SetString(PyExc_TypeError, errstr);
        return -1;
      }
      if (arr->nd == 1) {
        *p->tag5 = PyArray_DIM(arr, 0);
      } else {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_TypeError, "Array needs to be scalar or "
                        "one-dimensional.");
        return -1;
      }
      /* Type conversion madness */
      switch (arr->descr->type_num) {
      case NPY_FLOAT:
        for (i = 0; i < *p->tag5; i++) {
          ((double *) p->ptr)[i] = ((npy_float *) PyArray_DATA(arr))[i];
        }
        break;
      case NPY_DOUBLE:
        for (i = 0; i < *p->tag5; i++) {
          ((double *) p->ptr)[i] = ((npy_double *) PyArray_DATA(arr))[i];
        }
        break;
      default:
        PyErr_SetString(PyExc_TypeError, "Don't know how to convert from "
                        "numpy float type.");
        return -1;
      }
      Py_DECREF(arr);
    }
    break;
  case PK_INT_LIST:
#if PY_MAJOR_VERSION >= 3
    if (PyLong_Check(value)) {
#else
    if (PyInt_Check(value)) {
#endif
      *p->tag5 = 1;
#if PY_MAJOR_VERSION >= 3
      *((int*) p->ptr) = PyLong_AS_LONG(value);
#else
      *((int*) p->ptr) = PyInt_AS_LONG(value);
#endif
    } else {
      PyArray_Converter(value, (PyObject**) &arr);
      if (!PyArray_ISINTEGER(arr)) {
        sprintf(errstr,
                "Property '%s' of section '%s' should be a list of "
                "integers.\n",
                p->name, p->parent->name);
        PyErr_SetString(PyExc_TypeError, errstr);
        return -1;
      }
      if (arr->nd == 1) {
        *p->tag5 = PyArray_DIM(arr, 0);
      } else {
        Py_DECREF(arr);
        PyErr_SetString(PyExc_TypeError, "Array needs to be scalar or "
                        "one-dimensional.");
        return -1;
      }
      /* Type conversion madness */
      switch (arr->descr->type_num) {
      case NPY_INT:
        for (i = 0; i < *p->tag5; i++) {
          ((int *) p->ptr)[i] = ((npy_int*) PyArray_DATA(arr))[i];
        }
        break;
      case NPY_LONG:
        for (i = 0; i < *p->tag5; i++) {
          ((int *) p->ptr)[i] = ((npy_long*) PyArray_DATA(arr))[i];
        }
        break;
      default:
        PyErr_SetString(PyExc_TypeError, "Don't know how to convert from "
                        "numpy integer type.");
        return -1;
      }
      Py_DECREF(arr);
    }
    break;
  case PK_FORTRAN_STRING_LIST:
#if PY_MAJOR_VERSION >= 3
    if (PyUnicode_Check(value)) {
#else
    if (PyString_Check(value)) {
#endif
      *p->tag5 = 1;
      pystring_to_fstring(value, (char*) p->ptr, p->tag);
    } else {
      arr = PyArray_FROMANY(value, NPY_STRING, 1, 1, NPY_DEFAULT);
      if (!arr)
        return -1;
      *p->tag5 = PyArray_DIM(arr, 0);
      k = PyArray_STRIDE(arr, 0);
      str2 = (char *) p->ptr;
      for (i = 0; i < *p->tag5; i++) {
        str = PyArray_GETPTR1(arr, i);
        cstring_to_fstring(str, k, str2, p->tag);
        str2 += p->tag;
      }
      Py_DECREF(arr);
    }
    break;
  case PK_ENUM:
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(value)) {
#else
    if (!PyString_Check(value)) {
#endif
      sprintf(errstr,
              "Property '%s' of section '%s' should be a string.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    }
#if PY_MAJOR_VERSION >= 3
    bvalue = PyUnicode_AsASCIIString(value);
    str = PyBytes_AS_STRING(bvalue);
#else
    str = PyString_AS_STRING(value);
#endif
    j = -1;
    for (i = 0; i < p->tag; i++) {
      if (!strcmp(str, p->tag4 + i*p->tag2))
        j = i;
    }
    if (j < 0) {
      sprintf(errstr, "[ptrdict_set_property] Error: Could not find key "
              "'%s' in property '%s' of section '%s'.\n",
              str, p->name, p->parent->name);
#if PY_MAJOR_VERSION >= 3
      Py_DECREF(bvalue);
#endif
      PyErr_SetString(PyExc_ValueError, errstr);
      return -1;
    } else {
      *((int*) p->ptr) = j;
    };
#if PY_MAJOR_VERSION >= 3
    Py_DECREF(bvalue);
#endif
    break;
  case PK_ARRAY1D:
    if (!PyArray_Check(value)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 1d array "
              "of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    arr = (PyArrayObject *) value;
    if (!PyArray_ISFLOAT(arr)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 1d array "
              "of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    if (arr->nd != 1) {
      PyErr_SetString(PyExc_TypeError, "Array needs to be 1-dimensional.");
      return -1;
    }
    if (PyArray_DIM(arr, 0) != p->tag) {
      sprintf(errstr, "Wrong dimensions: Array needs to be of length %i.", 
              p->tag);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    /* Type conversion madness */
    switch (arr->descr->type_num) {
    case NPY_FLOAT:
      for (i = 0; i < p->tag; i++) {
        ((double*) p->ptr)[i] = ((npy_float *) PyArray_DATA(arr))[i];
      }
      break;
    case NPY_DOUBLE:
      for (i = 0; i < p->tag; i++) {
        ((double*) p->ptr)[i] = ((npy_double *) PyArray_DATA(arr))[i];
      }
      break;
    default:
      PyErr_SetString(PyExc_TypeError, "Don't know how to convert from "
                      "numpy float type.");
      return -1;
    }
    break;
  case PK_ARRAY2D:
    if (!PyArray_Check(value)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 2d array "
              "of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    arr = (PyArrayObject *) value;
    if (!PyArray_ISFLOAT(arr)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 2d array "
              "of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    if (arr->nd != 2) {
      PyErr_SetString(PyExc_TypeError, "Array needs to be 2-dimensional.");
      return -1;
    }
    if (PyArray_DIM(arr, 0) != p->tag || PyArray_DIM(arr, 1) != p->tag2) {
      sprintf(errstr, "Wrong dimensions: Array needs to be %ix%i.", p->tag, p->tag2);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    /* Type conversion madness */
    switch (arr->descr->type_num) {
    case NPY_FLOAT:
      for (i = 0; i < p->tag; i++) {
        for (j = 0; j < p->tag2; j++) {
          ((double*) p->ptr)[i + j*p->tag] = ((npy_float *) PyArray_DATA(arr))[j + i*p->tag2];
        }
      }
      break;
    case NPY_DOUBLE:
      for (i = 0; i < p->tag; i++) {
        for (j = 0; j < p->tag2; j++) {
          ((double*) p->ptr)[i + j*p->tag] = ((npy_double *) PyArray_DATA(arr))[j + i*p->tag2];
        }
      }
      break;
    default:
      PyErr_SetString(PyExc_TypeError, "Don't know how to convert from "
                      "numpy float type.");
      return -1;
    }
    break;
  case PK_ARRAY3D:
    if (!PyArray_Check(value)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 3d array "
              "of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    arr = (PyArrayObject *) value;
    if (!PyArray_ISFLOAT(arr)) {
      sprintf(errstr,
              "Property '%s' of section '%s' should be a 3d array "
              "of floats.\n",
              p->name, p->parent->name);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    if (arr->nd != 3) {
      PyErr_SetString(PyExc_TypeError, "Array needs to be 3-dimensional.");
      return -1;
    }
    if (PyArray_DIM(arr, 0) != p->tag || PyArray_DIM(arr, 1) != p->tag2 ||
        PyArray_DIM(arr, 2) != p->tag3) {
      sprintf(errstr, "Wrong dimensions: Array needs to be %ix%ix%i.",
              p->tag, p->tag2, p->tag3);
      PyErr_SetString(PyExc_TypeError, errstr);
      return -1;
    }
    /* Type conversion madness */
    switch (arr->descr->type_num) {
    case NPY_FLOAT:
      for (i = 0; i < p->tag; i++) {
        for (j = 0; j < p->tag2; j++) {
          for (k = 0; k < p->tag3; k++) {
            ((double*) p->ptr)[i + (j + k*p->tag2)*p->tag] = 
              ((npy_float *) PyArray_DATA(arr))[k + (j + i*p->tag2)*p->tag3];
          }
        }
      }
      break;
    case NPY_DOUBLE:
      for (i = 0; i < p->tag; i++) {
        for (j = 0; j < p->tag2; j++) {
          for (k = 0; k < p->tag3; k++) {
            ((double*) p->ptr)[i + (j + k*p->tag2)*p->tag] = 
              ((npy_double *) PyArray_DATA(arr))[k + (j + i*p->tag2)*p->tag3];
          }
        }
      }
      break;
    default:
      PyErr_SetString(PyExc_TypeError, "Don't know how to convert from "
                      "numpy float type.");
      return -1;
    }
    //      memcpy(p->ptr, data, p->tag*p->tag2*p->tag3*sizeof(double));
    break;
  default:
    sprintf(errstr, "Internal error: Unknown type with id %i encountered in section.", p->kind);
    PyErr_SetString(PyExc_TypeError, errstr);
    return -1;
    break;
  }

  return 0;
}


int
pydict_to_ptrdict(PyObject *dict, section_t *s)
{
  char *key;
  char errstr[1024];

  property_t *p;
  section_t *child;

  Py_ssize_t pos;
  PyObject *pkey, *value;

#ifdef DEBUG
  printf("[pydict_to_ptrdict] %p %p\n", dict, s);
  printf("[pydict_to_ptrdict] size(dict) = %i\n", PyDict_Size(dict));
#endif

  pos = 0;
  while (PyDict_Next(dict, &pos, &pkey, &value)) {
#if PY_MAJOR_VERSION >= 3
    if (!PyUnicode_Check(pkey)) {
#else
    if (!PyString_Check(pkey)) {
#endif
      PyErr_SetString(PyExc_TypeError, "Dictionary key needs to be string.");
      return -1;
    }

#if PY_MAJOR_VERSION >= 3
    PyObject *bpkey = PyUnicode_AsASCIIString(pkey);
    key = PyBytes_AS_STRING(bpkey);
#else
    key = PyString_AS_STRING(pkey);
#endif

#ifdef DEBUG
    printf("[pydict_to_ptrdict] key = %s\n", key);
#endif

    // Look for property with name *key*
    p = s->first_property;
    while (p != NULL && strcmp(p->name, key)) {
      p = p->next;
    }

    if (p) {
      // Property found

      if (pyobject_to_property(value, p)) {
#if PY_MAJOR_VERSION >= 3
        Py_DECREF(bpkey);
#endif
        return -1;
      }
    } else {
#ifdef DEBUG
      printf("[pydict_to_ptrdict] Property '%s' not found.\n", key);
#endif

      // No property found, check if there is a section with that name
      child = s->first_child;

      while (child != NULL && strcmp(child->name, key)) {
        child = child->next;
      }

      if (child) {

#ifdef DEBUG
        printf("[pydict_to_ptrdict] Section '%s' found.\n", key);
#endif

        // Value should be a dictionary
        if (!PyDict_Check(value)) {
#if PY_MAJOR_VERSION >= 3
          Py_DECREF(bpkey);
#endif
          return -1;
        }

#ifdef DEBUG
        printf("[pydict_to_ptrdict] Child: %s\n", child->name);
#endif

        child->provided = TRUE;
        if (child->provided_notification)
          *child->provided_notification = TRUE;
      
        pydict_to_ptrdict(value, child);

      } else {

        /* Ignore this property if it starts with '__' */
        if (key[0] != '_' || key[1] != '_') {
          sprintf(errstr, "Could not find property '%s' of section '%s'.",
                  key, s->name);
#if PY_MAJOR_VERSION >= 3
          Py_DECREF(bpkey);
#endif
          PyErr_SetString(PyExc_TypeError, errstr);
          return -1;
        }
      }
    }

#if PY_MAJOR_VERSION >= 3
    Py_DECREF(bpkey);
#endif
  }

  return 0;
}


PyObject *
property_to_pyobject(property_t *p)
{
  if (!p->ptr) Py_RETURN_NONE;

  int i, j, k;
  double *data;

  npy_intp dims[3];

  PyObject **odata;
  PyArrayObject *arr;

  PyObject *r = NULL;

  switch (p->kind) {
  case PK_INT:
#if PY_MAJOR_VERSION >= 3
    r = PyLong_FromLong(*((int*) p->ptr));
#else
    r = PyInt_FromLong(*((int*) p->ptr));
#endif
    break;
  case PK_DOUBLE:
    r = PyFloat_FromDouble(*((double*) p->ptr));
    break;
  case PK_BOOL:
    r = PyBool_FromLong(*((BOOL*) p->ptr));
    break;
  case PK_LIST:
    if (*p->tag5 == 1) {
      r = PyFloat_FromDouble(*((double*) p->ptr));
    } else {
      dims[0] = *p->tag5;
      arr = (PyArrayObject*) PyArray_SimpleNew(1, (npy_intp*) dims,
                                               NPY_DOUBLE);
      data = (double *) PyArray_DATA(arr);
      for (i = 0; i < *p->tag5; i++) {
        data[i] = ((double*) p->ptr)[i];
      }
      r = (PyObject*) arr;
    }
    break;
  case PK_FORTRAN_STRING_LIST:
    if (*p->tag5 == 1) {
      r = fstring_to_pystring((char*) p->ptr, p->tag);
    } else {
      dims[0] = *p->tag5;
      arr = (PyArrayObject*) PyArray_SimpleNew(1, (npy_intp*) dims,
                                               NPY_OBJECT);
      odata = (PyObject **) PyArray_DATA(arr);
      for (i = 0; i < *p->tag5; i++) {
        odata[i] = fstring_to_pystring(((char*) p->ptr + i*p->tag), p->tag);
      }
      r = (PyObject*) arr;
    }
    break;
  case PK_ARRAY1D:
    dims[0] = p->tag;
    arr = (PyArrayObject*) PyArray_SimpleNew(1, (npy_intp*) dims, NPY_DOUBLE);
    data = (double *) PyArray_DATA(arr);
    for (i = 0; i < p->tag; i++) {
      data[i] = ((double*) p->ptr)[i];
    }
    //        memcpy(data, p->ptr, p->tag*p->tag2*sizeof(double));
    r = (PyObject*) arr;
    break;
  case PK_ARRAY2D:
    dims[0] = p->tag;
    dims[1] = p->tag2;
    arr = (PyArrayObject*) PyArray_SimpleNew(2, (npy_intp*) dims, NPY_DOUBLE);
    data = (double *) PyArray_DATA(arr);
    for (i = 0; i < p->tag; i++) {
      for (j = 0; j < p->tag2; j++) {
        data[j + i*p->tag2] = ((double*) p->ptr)[i + j*p->tag];
      }
    }
    //        memcpy(data, p->ptr, p->tag*p->tag2*sizeof(double));
    r = (PyObject*) arr;
    break;
  case PK_ARRAY3D:
    dims[0] = p->tag;
    dims[1] = p->tag2;
    dims[2] = p->tag3;
    arr = (PyArrayObject*) PyArray_SimpleNew(3, (npy_intp*) dims, NPY_DOUBLE);
    data = (double *) PyArray_DATA(arr);
    for (i = 0; i < p->tag; i++) {
      for (j = 0; j < p->tag2; j++) {
        for (k = 0; k < p->tag3; k++) {
          data[k + (j + i*p->tag2)*p->tag3] = 
            ((double*) p->ptr)[i + (j + k*p->tag2)*p->tag];
        }
      }
    }
    //        memcpy(data, p->ptr, p->tag*p->tag2*p->tag3*sizeof(double));
    r = (PyObject*) arr;
    break;
  default:
    PyErr_SetString(PyExc_TypeError, "Internal error: Unknown type encountered in section.");
    break;
  }

  return r;
}