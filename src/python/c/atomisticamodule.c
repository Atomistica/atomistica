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
#include <numpy/arrayobject.h>

#include <stddef.h>

#include "atomisticamodule.h"

#include "coulomb_factory_c.h"
#include "potentials_factory_c.h"

#include <fenv.h>

/*
#include "forcomp.c"
#include "pyf90obj.c"

#include "particles.c"
#include "neighbors.c"
#include "potential.c"
*/

#define MAX_COULOMB_NAME 100
#define MAX_POTENTIAL_NAME 100

/* Global methods
 */

int has_started = 0;

PyObject *
py_atomistica_startup(PyObject *self, PyObject *args)
{
  if (!has_started) {
    atomistica_startup(-1);

    has_started = 1;
  }

  Py_RETURN_NONE;
}


PyObject *
py_atomistica_shutdown(PyObject *self, PyObject *args)
{
  atomistica_shutdown();

  Py_RETURN_NONE;
}


PyObject *
py_set_logfile(PyObject *self, PyObject *args)
{
  if (!has_started) {
    PyErr_SetString(PyExc_RuntimeError, "Please run _atomistica.startup() "
                    "before changing the log file.");
    return NULL;
  }

  char *fn;
  if (!PyArg_ParseTuple(args, "s", &fn))
    return NULL;

  f_logging_start(fn);

  Py_RETURN_NONE;
}


static PyMethodDef module_methods[] = {
  { "startup", py_atomistica_startup, METH_NOARGS,
    "File which to write log information to." },
  { "shutdown", py_atomistica_shutdown, METH_NOARGS,
    "Write timings and close log file." },
  { "set_logfile", py_set_logfile, METH_VARARGS,
    "Set name of log file." },
  { "pair_distribution", py_pair_distribution, METH_VARARGS,
    "Compute pair distribution function." },
  { "angle_distribution", py_angle_distribution, METH_VARARGS,
    "Compute angular distribution functions." },
  { "bond_angles", py_bond_angles, METH_VARARGS,
    "Compute moments of the bond angle distribution (per-atom)." },
  { NULL, NULL, 0, NULL }  /* Sentinel */
};


/*
 * Module initialization
 */

#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

/*
 * Module declaration
 */

#if PY_MAJOR_VERSION >= 3
    #define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
        static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
        ob = PyModule_Create(&moduledef);
    #define MOD_RETURN_ERROR return NULL
#else
    #define MOD_INIT(name) PyMODINIT_FUNC init##name(void)
    #define MOD_DEF(ob, name, methods, doc) \
        ob = Py_InitModule3(name, methods, doc);
    #define MOD_RETURN_ERROR return
#endif

static PyTypeObject
coulomb_types[N_COULOMB_CLASSES];
static PyTypeObject
potential_types[N_POTENTIAL_CLASSES];

MOD_INIT(_atomistica)
{
    PyObject* m;
    int i;

#if 0
    /* Uncomment to enable floating-point exception */
    //    int feenableexcept();
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

    import_array();

    if (PyType_Ready(&particles_type) < 0)
      MOD_RETURN_ERROR;

    if (PyType_Ready(&neighbors_type) < 0)
      MOD_RETURN_ERROR;

    MOD_DEF(m, "_atomistica", module_methods,
            "Interface to the Atomistica interatomic potential library.")
    if (m == NULL)
      MOD_RETURN_ERROR;

    Py_INCREF(&particles_type);
    PyModule_AddObject(m, "Particles", (PyObject *) &particles_type);

    Py_INCREF(&neighbors_type);
    PyModule_AddObject(m, "Neighbors", (PyObject *) &neighbors_type);

    for (i = 0; i < N_COULOMB_CLASSES; i++) {
      coulomb_types[i] = coulomb_type;
      coulomb_types[i].tp_name = malloc(MAX_COULOMB_NAME);
      strncpy((char*) coulomb_types[i].tp_name, "_atomistica.",
              MAX_COULOMB_NAME);
      strncat((char*) coulomb_types[i].tp_name, coulomb_classes[i].name,
              MAX_COULOMB_NAME);

      if (PyType_Ready(&coulomb_types[i]) < 0)
        MOD_RETURN_ERROR;

      Py_INCREF(&coulomb_types[i]);
      PyModule_AddObject(m, coulomb_classes[i].name,
                         (PyObject *) &coulomb_types[i]);
    }

    for (i = 0; i < N_POTENTIAL_CLASSES; i++) {
      potential_types[i] = potential_type;
      potential_types[i].tp_name = malloc(MAX_POTENTIAL_NAME);
      strncpy((char*) potential_types[i].tp_name, "_atomistica.",
	      MAX_POTENTIAL_NAME);
      strncat((char*) potential_types[i].tp_name, potential_classes[i].name,
	      MAX_POTENTIAL_NAME);

      if (PyType_Ready(&potential_types[i]) < 0)
        MOD_RETURN_ERROR;

      Py_INCREF(&potential_types[i]);
      PyModule_AddObject(m, potential_classes[i].name,
			 (PyObject *) &potential_types[i]);
    }

#if PY_MAJOR_VERSION >= 3
    return m;
#endif
}
