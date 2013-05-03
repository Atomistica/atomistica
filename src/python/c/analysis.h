/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
   ====================================================================== */
#ifndef __ANALYSIS_H
#define __ANALYSIS_H

#include <Python.h>

PyObject *py_pair_distribution(PyObject *, PyObject *);
PyObject *py_angle_distribution(PyObject *, PyObject *);

#endif
