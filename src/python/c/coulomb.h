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

#ifndef __COULOMB_H
#define __COULOMB_H

#include <Python.h>

#include "ptrdict.h"
#include "coulomb_factory_c.h"


typedef struct {
  PyObject_HEAD

  /* Pointer to F90-object */
  void *f90obj;

  /* Pointer to the F90-member descriptor */
  section_t *f90members;

  /* Pointer to the F90-class descriptor */
  coulomb_class_t *f90class;

} coulomb_t;


extern PyTypeObject coulomb_type;


#endif
