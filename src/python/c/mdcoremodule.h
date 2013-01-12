/* ======================================================================
   MDCORE - Interatomic potential library
   https://github.com/pastewka/mdcore
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
#ifndef __MDCOREMODULE_H_
#define __MDCOREMODULE_H_

#include <Python.h>

#include <stddef.h>

#include "ptrdict.h"
#include "py_f.h"


#include "particles.h"
#include "neighbors.h"
#include "potential.h"
#include "analysis.h"


#define TYPE_REAL_ATTR      1
#define TYPE_REAL3_ATTR     2
#define TYPE_REAL3x3_ATTR   3
#define TYPE_REAL           4
#define TYPE_INTEGER_ATTR   5
#define TYPE_INTEGER        6
#define TYPE_LOGICAL        7
#define TYPE_REAL3          8
#define TYPE_REAL6          9
#define TYPE_REAL3x3        10


/* Prototypes, implementation found in python_helper.f90 */

/* data_t */
BOOL data_exists(void *, char *, int *);
int data_get_len(void *);
void real_ptr_by_name(void *, char *, void **, int *);
void integer_ptr_by_name(void *, char *, void **, int *);
void realx_ptr_by_name(void *, char *, void **, int *);
void realxxx_ptr_by_name(void *, char *, void **, int *);
void real_attr_by_name(void *, char *, void **, int *);
void real3_attr_by_name(void *, char *, void **, int *);
void real3x3_attr_by_name(void *, char *, void **, int *);

/* general */
void units_init(int);
void mdcore_startup(int);
void mdcore_shutdown(void);
void timer_print_to_log(void);

#endif
