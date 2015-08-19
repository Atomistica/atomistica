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

#ifndef __MATERIALS_H
#define __MATERIALS_H

#include <stdbool.h>

/*
 * IF YOU MODIFY THIS STRUCTURE, *ALWAYS* ALSO MODIFY THE CORRESPONDING
 * STRUCTURE IN materials.f90
 */

struct notb_element_t {

  _Bool exists;

  char name[2];     /* name of element */
  char cname[10];   /* common name of element */
  int elem;         /* number of element (official) */
  int no;           /* number of orbitals */
  int l[9];         /* angular momenta of orbitals */
  int lmax;         /* maximum angular momentum */
  double e[9];      /* orbital energies [ e(1:no) ] */
  double el_max;    /* max number of valence electrons on an atom */
  double U;         /* Hubbard U */
  double q0;        /* charge (nr of electrons in neutral) */
       
  /*
   * internal bookkeeping
   */

  int o1;           /* index of the first orbital */
  int enr;          /* element number in the internal book-keeping */

  /*
   * spin-related variables
   */

  _Bool spin;       /* spin parameter set */
  double W[9];      /* W parameter values */

};

#endif
