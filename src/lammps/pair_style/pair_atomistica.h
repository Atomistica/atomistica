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

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(mdcore,PairAtomistica)
PairStyle(atomistica,PairAtomistica)

#else

#ifndef LMP_PAIR_ATOMISTICA_H
#define LMP_PAIR_ATOMISTICA_H

#include "pair.h"

#include "ptrdict.h"
#include "potentials_factory_c.h"

namespace LAMMPS_NS {

class PairAtomistica : public Pair {
 public:
  PairAtomistica(class LAMMPS *);
  ~PairAtomistica();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 //private:
  char *name_;                     // name of the potential
  char *fn_;                       // file name with potential parameters
  int maxlocal_;                   // size of numneigh, firstneigh arrays

  intptr_t *Atomistica_seed_;
  intptr_t *Atomistica_last_;
  int *Atomistica_neighb_;
  int *Atomistica_neighb_endptr_;
  int Atomistica_nneighb_;

  int *mask_;

  // neighbor list range
  double rcmax_;
  double **rcmaxsq_;
  // width of the communication border
  double rcghost_;
  // interaction range
  double rangemax_;

  // pointer to the -member descriptor
  section_t *members_;

  // pointer to the -class descriptor
  potential_class_t *class_;

  // particles, neighbors and potential objects
  void *particles_,*neighbors_,*potential_;

  void Atomistica_neigh();
  void FAtomistica(int, int);

  void allocate();
};

int error2lmp(Error *, const char *, int, int);

}

extern "C" {

  void get_atomistica_pair_style_git_ident(char *);

  void particles_new(void **self);            // allocate particles object
  void particles_free(void *self);            // free particles object
  
  void particles_init(void *self);            // initialize particles object
  void particles_del(void *self);             // finalize particles object

  void particles_set_element(void *self, char *el_str, int nel, int el_no, 
                             int *Z, int *error);
  void particles_set_pointers(void *self, int nat, int natloc, int maxnatloc,
                              void *tag, void *el, void *r);

  void particles_get_border(void *self, double *border);
  void particles_get_interaction_range(void *self, int el1, int el2,
                                       double *border);


  void neighbors_new(void **self);            // allocate neighbors object
  void neighbors_free(void *self);            // free neighbors object
  
  void neighbors_init(void *self);            // initialize neighbors object
  void neighbors_del(void *self);             // finalize neighbors object

  void neighbors_set_pointers(void *self, int natloc, void *seed, void *last,
                              int neighbors_size, void *neighbors);

  void neighbors_get_cutoff(void *self, int el1, int el2, double *cutoff);
  void neighbors_dump_cutoffs(void *self, void *p);

  void atomistica_startup(int);
  void atomistica_shutdown(void);
  void timer_print_to_log(void);
  void get_full_error_string(char *);

}

#endif
#endif
