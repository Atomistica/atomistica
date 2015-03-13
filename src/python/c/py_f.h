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
void py_to_error(char *file, int line, int *ierror);


/* Initialize an F90-object from a Python dictionary */

int pyobject_to_property(PyObject *value, property_t *p);
int pydict_to_ptrdict(PyObject *dict, section_t *s);
PyObject *property_to_pyobject(property_t *p);

#endif
