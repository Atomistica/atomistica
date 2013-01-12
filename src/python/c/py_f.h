/* ======================================================================
   MDCORE - Interatomic potential library
   https://github.com/pastewka/mdcore
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
#ifndef _PY_F_H
#define _PY_F_H

#include <Python.h>

#include "ptrdict.h"


#define ERROR_NONE          0

#define ERRSTRLEN       10000


#define DOUBLEP(a) ((double*)((a)->data))
#define BOOLP(a)   ((int*)((a)->data))


/* String conversion */

void cstring_to_fstring(char *cstr, int clen, char *forstr, int forlen);
void pystring_to_fstring(PyObject *pystr, char *forstr, int len);
PyObject* fstring_to_pystring(char *forstr, int len);


/* Obtain error information from libAtoms error module */

void get_full_error_string(char *);


/* Pass an error from an F90-object to the Python runtime */

int error_to_py(int ierror);


/* Initialize an F90-object from a Python dictionary */

int pydict_to_ptrdict(PyObject *dict, section_t *s);

#endif
