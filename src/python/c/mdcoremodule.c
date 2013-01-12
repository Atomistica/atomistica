/* ======================================================================
   MDCORE - Interatomic potential library
   https://github.com/pastewka/mdcore
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL MDCORE_ARRAY_API
#include <numpy/arrayobject.h>

#include <stddef.h>

#include "mdcoremodule.h"

#include "potentials_factory_c.h"

#include <fenv.h>

/*
#include "forcomp.c"
#include "pyf90obj.c"

#include "particles.c"
#include "neighbors.c"
#include "potential.c"
*/

#define MAX_POTENTIAL_NAME 100

/* Global methods
 */

int has_started = 0;

PyObject *
py_mdcore_startup(PyObject *self, PyObject *args)
{
  if (!has_started) {
    mdcore_startup(-1);

    has_started = 1;
  }

  Py_RETURN_NONE;
}


PyObject *
py_mdcore_shutdown(PyObject *self, PyObject *args)
{
  mdcore_shutdown();

  Py_RETURN_NONE;
}


#if 0
PyObject *
read_atoms(PyObject *self, PyObject *args)
{
  char *at_fn;
  particles_t *a;
  int i;
  FOR_INTEGER ierror = ERROR_NONE;

  if (!PyArg_ParseTuple(args, "s", &at_fn))
    return NULL;

  a = (particles_t*) particles_type.tp_new(&particles_type, NULL, NULL);
  i = strlen(at_fn);
  CALL5(native_io, read_atoms, a->obj, at_fn, NULL, &ierror, i);
  if (error_to_py(ierror))
    return NULL;

  a->initialized = 1;

  particles_update_elements(a, NULL);

  return (PyObject*) a;
}
#endif


static PyMethodDef module_methods[] = {
  { "startup", py_mdcore_startup, METH_NOARGS,
    "File which to write log information to." },
  { "shutdown", py_mdcore_shutdown, METH_NOARGS,
    "Write timings and close log file." },
  { "pair_distribution", py_pair_distribution, METH_VARARGS,
    "Compute pair distribution function." },
  { "angle_distribution", py_angle_distribution, METH_VARARGS,
    "Compute angular distribution functions." },
  /*
  { "read_atoms", read_atoms, METH_VARARGS,
      "Read atom data from an .dat-file." },
  */
  { NULL, NULL, 0, NULL }  /* Sentinel */
};


/* Module initialization
 */

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif

static PyTypeObject
potential_types[N_CLASSES];

PyMODINIT_FUNC
init_mdcore(void)
{
    PyObject* m;
    int i;

#if 1
    /* Uncomment to enable floating-point exception */
    //    int feenableexcept();
    feenableexcept(FE_DIVBYZERO | FE_INVALID);
#endif

    import_array();

    if (PyType_Ready(&particles_type) < 0)
      return;

    if (PyType_Ready(&neighbors_type) < 0)
      return;

    m = Py_InitModule3("_mdcore", module_methods,
                       "Interface to the MDCORE molecular dynamics Fortran "
		       "kernel.");

    if (m == NULL)
      return;

    Py_INCREF(&particles_type);
    PyModule_AddObject(m, "Particles", (PyObject *) &particles_type);

    Py_INCREF(&neighbors_type);
    PyModule_AddObject(m, "Neighbors", (PyObject *) &neighbors_type);

    for (i = 0; i < N_CLASSES; i++) {
      potential_types[i] = potential_type;
      potential_types[i].tp_name = malloc(MAX_POTENTIAL_NAME);
      strncpy((char*) potential_types[i].tp_name, "mdcore.",
	      MAX_POTENTIAL_NAME);
      strncat((char*) potential_types[i].tp_name, potential_classes[i].name,
	      MAX_POTENTIAL_NAME);

      if (PyType_Ready(&potential_types[i]) < 0)
        return;

      Py_INCREF(&potential_types[i]);
      PyModule_AddObject(m, potential_classes[i].name,
			 (PyObject *) &potential_types[i]);
    }
}
