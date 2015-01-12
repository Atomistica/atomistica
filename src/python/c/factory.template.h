/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level Atomistica directory.

   Copyright (2005-2013) Fraunhofer IWM

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
/*
 * %(disclaimer)s
 */

#ifndef __%(name)s_DISPATCH_H_
#define __%(name)s_DISPATCH_H_

#include "ptrdict.h"

#define N_POTENTIAL_CLASSES %(n_classes)i

/*
 * Class definition
 */

typedef struct __%(name)s_class_t {

  char name[MAX_NAME+1];
  void (*new_instance)(void **, section_t *, section_t **);
  void (*free_instance)(void *);

  void (*register_data)(void *, void *, int *);
  void (*init)(void *, int *);
  void (*bind_to)(void *, void *, void *, int *);
  void (*set_Coulomb)(void *, void *, int *);
  void (*energy_and_forces)(void *, void *, void *, double *, double *,
                            double *, int *, double *, double *, double *,
                            double *, double *, double *, int *);

} %(name)s_class_t;

extern %(name)s_class_t %(name)s_classes[N_POTENTIAL_CLASSES];

#endif


