/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
/*
 * %(disclaimer)s
 */

#ifndef __%(name)s_DISPATCH_H_
#define __%(name)s_DISPATCH_H_

#include "ptrdict.h"

#define N_CLASSES %(n_classes)i

/*
 * Class definition
 */

typedef struct __%(name)s_class_t {

  char name[MAX_NAME+1];
  void (*new_instance)(void **, section_t *, section_t **);
  void (*free_instance)(void *);

  void (*register_data)(void *, void *, int *);
  void (*init)(void *);
  void (*bind_to)(void *, void *, void *, int *);
  void (*energy_and_forces)(void *, void *, void *, double *, double *,
			    double *, double *, double *, double *, double *,
			    double *, int *);

} %(name)s_class_t;

extern %(name)s_class_t %(name)s_classes[N_CLASSES];

#endif


