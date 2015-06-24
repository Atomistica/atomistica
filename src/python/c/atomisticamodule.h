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
#ifndef __ATOMISTICAMODULE_H_
#define __ATOMISTICAMODULE_H_

#include <Python.h>

#include <stddef.h>

#include "ptrdict.h"
#include "py_f.h"


#include "particles.h"
#include "neighbors.h"
#include "coulomb.h"
#include "potential.h"
#include "analysis.h"


#define TYPE_REAL_ATTR      1
#define TYPE_REAL3_ATTR     2
#define TYPE_REAL3x3_ATTR   3
#define TYPE_REAL           4
#define TYPE_INTEGER_ATTR   5
#define TYPE_INTEGER3_ATTR  6
#define TYPE_INTEGER        7
#define TYPE_LOGICAL        8
#define TYPE_REAL3          9
#define TYPE_REAL6          10
#define TYPE_REAL3x3        11


/* Prototypes, implementation found in python_helper.f90 */

/* data_t */
BOOL f_data_exists(void *, char *, int *);
int data_get_len(void *);
void real_ptr_by_name(void *, char *, void **, int *);
void integer_ptr_by_name(void *, char *, void **, int *);
void realx_ptr_by_name(void *, char *, void **, int *);
void realxxx_ptr_by_name(void *, char *, void **, int *);
void integer_attr_by_name(void *, char *, void **, int *);
void integer3_attr_by_name(void *, char *, void **, int *);
void real_attr_by_name(void *, char *, void **, int *);
void real3_attr_by_name(void *, char *, void **, int *);
void real3x3_attr_by_name(void *, char *, void **, int *);

/* general */
void units_init(int);
void atomistica_startup(int);
void atomistica_shutdown(void);
void timer_print_to_log(void);

#endif
