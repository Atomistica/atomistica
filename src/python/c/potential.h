/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
#ifndef __POTENTIAL_H
#define __POTENTIAL_H

#include <Python.h>

#include "ptrdict.h"
#include "potentials_factory_c.h"


typedef struct {
  PyObject_HEAD

  /* Pointer to F90-object */
  void *f90obj;

  /* Pointer to the F90-member descriptor */
  section_t *f90members;

  /* Pointer to the F90-class descriptor */
  potential_class_t *f90class;

} potential_t;


extern PyTypeObject potential_type;

#endif
