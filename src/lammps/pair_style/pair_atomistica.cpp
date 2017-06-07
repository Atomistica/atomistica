#define ATOMISTICA_PAIR_STYLE_GIT_IDENT "$Id$"

#ifndef DEFINE_GIT_IDENT

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

/* ----------------------------------------------------------------------
   Contributing author: Tim Kunze (FZDR), Lars Pastewka (Fh-IWM, JHU)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "pair_atomistica.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "update.h"

//#include "gfmd_grid.h"

#include "potentials_factory_c.h"

using namespace LAMMPS_NS;

#define ERROR_NONE 0
#define ERRSTRLEN 10000

/* ---------------------------------------------------------------------- */

extern "C" void get_atomistica_pair_style_git_ident(char *ident)
{
  strcpy(ident, ATOMISTICA_PAIR_STYLE_GIT_IDENT);
}

/* ---------------------------------------------------------------------- */

int LAMMPS_NS::error2lmp(Error *error, const char *fn, int line, int ierror)
{
  char errstr[ERRSTRLEN];

  if (ierror != ERROR_NONE) {
    get_full_error_string(errstr);
    error->all(fn,line,errstr);
    return 1;
  } else {
    return 0;
  }   
}

/* ---------------------------------------------------------------------- */

PairAtomistica::PairAtomistica(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;
  no_virial_fdotr_compute = 1;
  ghostneigh = 1;

  name_ = NULL;
  fn_ = NULL;

  rcmax_ = 0.0;
  rcghost_ = 0.0;
  rangemax_ = 0.0;

  maxlocal_ = 0;
  Atomistica_seed_ = NULL;
  Atomistica_last_ = NULL;
  Atomistica_nneighb_ = 0;
  Atomistica_neighb_ = NULL;

  mask_ = NULL;

  particles_new(&particles_);
  particles_init(particles_);

  neighbors_new(&neighbors_);
  neighbors_init(neighbors_);

  class_ = NULL;
  members_ = NULL;
  potential_ = NULL;

  atomistica_startup(-1);
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairAtomistica::~PairAtomistica()
{
  if (name_)  free(name_);
  if (fn_)  free(fn_);

  if (potential_) {
    class_->del(potential_);
    class_->free_instance(potential_);
  }

  memory->sfree(Atomistica_seed_);
  memory->sfree(Atomistica_last_);

  if (mask_)
    memory->sfree(mask_);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(rcmaxsq_);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
  }

  if (neighbors_) {
    neighbors_del(neighbors_);
    neighbors_free(neighbors_);
  }

  if (particles_) {
    particles_del(particles_);
    particles_free(particles_);
  }

  atomistica_shutdown();
}

/* ---------------------------------------------------------------------- */

void PairAtomistica::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  Atomistica_neigh();
  FAtomistica(eflag,vflag);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairAtomistica::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(rcmaxsq_,n+1,n+1,"pair:rcmaxsq");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairAtomistica::settings(int narg, char **arg)
{
  if (narg != 1 && narg != 2)
    error->all(FLERR,"pair_style atomistica expects potential name and "
           "configuration file as parameters");

  if (name_)  free(name_);
  if (fn_)  free(fn_);

  name_ = strdup(arg[0]);
  if (narg == 2)
    fn_ = strdup(arg[1]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairAtomistica::coeff(int narg, char **arg)
{
  int n = atom->ntypes;
  int map[n];

  if (!allocated)  allocate();

  if (narg != 2 + n) {
    char errstr[1024];
    sprintf(errstr,"Incorrect number of arguments for pair coefficients. "
        "There are %i atom types in this system.", n);
    error->all(FLERR,errstr);
  }

  // ensure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients; must be * *");

  // read args that map atom types to C and H
  // map[i] = which element (0,1) the Ith atom type is, -1 if NULL

  for (int i = 2; i < narg; i++) {
    map[i-1] = 0;
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-1] = -1;
      continue;
    } else {
      int Z, ierror;
      particles_set_element(particles_,arg[i],n,i-1,&Z,&ierror);
      error2lmp(error,FLERR,ierror);
      map[i-2] = Z;
    }
  }

  // clear setflag since coeff() called once with I,J = * *

  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i-1] >= 0 && map[j-1] >= 0) {
    setflag[i][j] = 1;
    count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients -> "
                 "count = 0");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairAtomistica::init_style()
{
  if (!allocated)
    error->all(FLERR,"Something wrong. pair atomistica not allocated.");

  if (strcmp(update->unit_style,"metal"))
    error->all(FLERR,"Pair style atomistica requires metal units");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style atomistica requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style atomistica requires newton pair on");

  // delete potential object, if already allocated

  if (potential_) {
    if (!class_)
      error->one(FLERR,"(Internal error) class_ is NULL, but potential_ is "
                       "not.");
    class_->del(potential_);
    class_->free_instance(potential_);
  }

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  // find potential class in Atomistica potential database

  class_ = NULL;
  for (int i = 0; i < N_POTENTIAL_CLASSES; i++) {
    if (!strcmp(name_,potential_classes[i].name))
      class_ = &potential_classes[i];
  }
  if (!class_) {
    char errstr[1024];
    sprintf(errstr,"Could not find potential '%s' in the Atomistica potential "
        "database",name_);
    error->all(FLERR,errstr);
  }

  // initialize  potential object

  section_t *zero = NULL;
  class_->new_instance(&potential_,zero,&members_);

  if (fn_) {
    ptrdict_read(members_, fn_);
  }

  // set pointers in particles object
  particles_set_pointers(particles_,atom->nlocal+atom->nghost,atom->nlocal,
                         atom->nmax,atom->tag,atom->type,&atom->x[0][0]);

  class_->init(potential_);

  int ierror;
  class_->bind_to(potential_,particles_,neighbors_,&ierror);
  error2lmp(error,FLERR,ierror);

  // dump all cutoffs to the Atomistica log file

  neighbors_dump_cutoffs(neighbors_,particles_);

  // determine width of ghost communication border

  particles_get_border(particles_,&rcghost_);
  comm->cutghostuser = MAX(comm->cutghostuser,rcghost_);

  rcmax_ = 0.0;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairAtomistica::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double rc, range;

  neighbors_get_cutoff(neighbors_,i,j,&rc);
  particles_get_interaction_range(particles_,i,j,&range);
  range = MAX(range, rc);

  rcmax_ = MAX(rc, rcmax_);
  rangemax_ = MAX(range, rangemax_);

  rcmaxsq_[i][j] = rcmaxsq_[j][i] = rc*rc;
  cutghost[i][j] = cutghost[j][i] = rc;

  return rc;
}

/* ----------------------------------------------------------------------
   create Atomistica neighbor list from main neighbor list
   Atomistica neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairAtomistica::Atomistica_neigh()
{
  int i,j,ii,jj,n,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (!list->ghost) {
    error->all(FLERR,"Atomistica needs neighbor list with ghost atoms.");
  }

  if (nall > maxlocal_) {
    maxlocal_ = atom->nmax;
    memory->sfree(Atomistica_seed_);
    memory->sfree(Atomistica_last_);
#ifdef GFMD_GRID_H
    if (atom->gfmd_flag)
      memory->sfree(mask_);
#endif
    Atomistica_seed_ = (intptr_t *)
      memory->smalloc(maxlocal_*sizeof(intptr_t),"Atomistica:Atomistica_seed");
    Atomistica_last_ = (intptr_t *)
      memory->smalloc(maxlocal_*sizeof(intptr_t),"Atomistica:Atomistica_last");
#ifdef GFMD_GRID_H
    if (atom->gfmd_flag)
      mask_ = (int *) memory->smalloc(maxlocal_*sizeof(int),"Atomistica::mask");
#endif
  }

  // set start values for neighbor array Atomistica_neighb
  for (i = 0; i < nall; i++) {
    Atomistica_seed_[i] = -1;
    Atomistica_last_[i] = -2;
    if (mask_)
      mask_[i] = 1;
  }

  inum = list->inum+list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // Map seed and last arrays to point to the appropriate position in the
  // native LAMMPS neighbor list. This avoids copying the full list.

  // Seed is reported relative to the lowest neighbor pointer value. What is
  // passed to Atomistica is a Fortran array that encloses all neighbors, from
  // the lowest to the highest pointer value. This is necessary because the
  // neigbor list in LAMMPS is not necessarily consecutive in memory.
  Atomistica_neighb_ = firstneigh[ilist[0]];
  Atomistica_neighb_endptr_ = NULL;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    Atomistica_neighb_ = std::min(Atomistica_neighb_, firstneigh[i]);
    Atomistica_neighb_endptr_ = std::max(Atomistica_neighb_endptr_,
                                         firstneigh[i]+numneigh[i]);
  }

  // Fill seed and last arrays.
  Atomistica_nneighb_ = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    Atomistica_seed_[i] = firstneigh[i]-Atomistica_neighb_+1;
    Atomistica_last_[i] = Atomistica_seed_[i]+numneigh[i]-1;
    Atomistica_nneighb_ += numneigh[i];
  }

#if 0
  // DEBUG: Check if neighbor list is symmetric
  for (i = 0; i < nall; i++) {
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    for (ii = Atomistica_seed_[i]-1; ii < Atomistica_last_[i]; ii++) {
      j = Atomistica_neighb_[ii]-1;

      // Check if i is neighbor of j
      n = 0;
      for (jj = Atomistica_seed_[j]-1; jj < Atomistica_last_[j]; jj++) {
    if (Atomistica_neighb_[jj]-1 == i) {
      n = 1;
    }
      }
      if (!n) {
    printf("i = %i, j = %i\n", i, j);
    printf("Neighbors of i\n");
    for (jj = Atomistica_seed_[i]-1; jj < Atomistica_last_[i]; jj++) {
      j = Atomistica_neighb_[jj]-1;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      printf("   %i  %f\n", j, sqrt(delx*delx+dely*dely+delz*delz));
    }
    printf("Neighbors of j\n");
    for (jj = Atomistica_seed_[j]-1; jj < Atomistica_last_[j]; jj++) {
      j = Atomistica_neighb_[jj]-1;
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      printf("   %i  %f\n", j, sqrt(delx*delx+dely*dely+delz*delz));
    }
    error->one(FLERR,"Neighbor list not symmetric");
      }
    }
  }
#endif
}

/* ----------------------------------------------------------------------
   Atomistica forces and energy
------------------------------------------------------------------------- */

void PairAtomistica::FAtomistica(int eflag, int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = NULL;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double epot,*epot_per_at,*wpot_per_at,wpot[3][3];

  memset(wpot, 0, 9*sizeof(double));

  epot_per_at = NULL;
  if (eflag_atom) {
    epot_per_at = &eatom[0];
  }

  wpot_per_at = NULL;
  if (vflag_atom) {
    wpot_per_at = &vatom[0][0];
  }

#ifdef GFMD_GRID_H
  if (atom->gfmd_flag) {
    for (int i = 0; i < nall; i++) mask_[i] = !FLAG_FROM_POW2_IDX(atom->gid[i]);
    mask = mask_;
  }
#endif

  // set pointers in particles object
  particles_set_pointers(particles_,nall,atom->nlocal,atom->nmax,tag,
                         type,&x[0][0]);

  // set pointers in neighbor list object
  neighbors_set_pointers(neighbors_,nall,Atomistica_seed_,Atomistica_last_,
                         Atomistica_neighb_endptr_-Atomistica_neighb_+1,Atomistica_neighb_);

  int ierror;
  epot = 0.0;
  class_->energy_and_forces(potential_,particles_,neighbors_,&epot,&f[0][0],
                            &wpot[0][0],mask,epot_per_at,wpot_per_at,&ierror);
  error2lmp(error,FLERR,ierror);

  if (evflag) {
    // update energies
    eng_vdwl += epot;

    // update virial
    virial[0] -= wpot[0][0];
    virial[1] -= wpot[1][1];
    virial[2] -= wpot[2][2];
    virial[3] -= 0.5*(wpot[1][0]+wpot[0][1]);
    virial[4] -= 0.5*(wpot[2][0]+wpot[0][2]);
    virial[5] -= 0.5*(wpot[2][1]+wpot[1][2]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairAtomistica::memory_usage()
{
  double bytes = 0.0;
  return bytes;
}

#endif
