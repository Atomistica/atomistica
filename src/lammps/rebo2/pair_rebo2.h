/* ======================================================================
   MDCORE - Interatomic potential library
   https://github.com/pastewka/mdcore
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
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

PairStyle(rebo2,PairREBO2)

#else

#ifndef LMP_PAIR_REBO2_H
#define LMP_PAIR_REBO2_H

#include "pair.h"
#include <vector>

namespace LAMMPS_NS {

class PairREBO2 : public Pair {
 public:
  PairREBO2(class LAMMPS *);
  ~PairREBO2();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();

 private:
  int maxlocal;                    // size of numneigh, firstneigh arrays
  int *map;                        // 0 (C), 1 (H), or -1 (NULL) for each type

  double cut3rebo;                 // maximum distance for 3rd REBO neigh

  int screened;

  int *REBO_numneigh;              // # of pair neighbors for each atom
  int **REBO_firstneigh;           // ptr to 1st neighbor of each atom
  double *closestdistsq;           // closest owned atom dist to each ghost
  int *REBO_seed;
  int *REBO_last;
  int *REBO_neighb;
  int REBO_nneighb;
  int *REBO_el;

  double rcmax[2][2],rcmaxsq[2][2];

  void *f90obj;

  void REBO_neigh();
  void FREBO(int, int);

  void allocate();
};

}

/*
 * Fortran-90 function prototypes
 */

extern "C" {

void rebo2_create(void **self);            // allocate object
void rebo2_destroy(void *self);            // free object
void rebo2_init(void *self);               // constructor
void rebo2_del(void *self);                // destructor
void rebo2_get_cutoff(void *self, double rcmaxsq[2][2]); 
void rebo2_energy_and_forces(void *self,
			     int natloc, int nat, int *tag, int *Z, double *x,
			     int *seed, int *last, int nneighb, int *neighb,
			     double *epot, double *f, double wpot[3][3],
			     double *wpot_per_at, //double *wpot_per_bond,
			     int *ierror);   // compute forces


void rebo2_scr_create(void **self);            // allocate object
void rebo2_scr_destroy(void *self);            // free object
void rebo2_scr_init(void *self);               // constructor
void rebo2_scr_del(void *self);                // destructor
void rebo2_scr_get_cutoff(void *self, double rcmaxsq[2][2]); 
void rebo2_scr_energy_and_forces(void *self,
				 int natloc, int nat, int *tag, int *Z, double *x,
				 int *seed, int *last, int nneighb, int *neighb,
				 double *epot, double *f, double wpot[3][3],
				 double *wpot_per_at, //double *wpot_per_bond,
				 int *ierror);   // compute forces

}

#endif
#endif
