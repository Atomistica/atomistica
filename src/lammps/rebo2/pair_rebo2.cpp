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
   Contributing author: Tim Kunze (FZDR), Lars Pastewka (Fh-IWM)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"
#include "pair_rebo2.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include<iostream>
//#include<vector>
#include "update.h"

using namespace LAMMPS_NS;

#define PGDELTA 1

/* ---------------------------------------------------------------------- */

PairREBO2::PairREBO2(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  one_coeff = 1;
  no_virial_fdotr_compute = 1;
  ghostneigh = 1;

  maxlocal = 0;
  screened = 0;
  REBO_seed = NULL;
  REBO_last = NULL;
  REBO_nneighb = 0;
  REBO_neighb = NULL;
  REBO_el = NULL;

  f90obj = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairREBO2::~PairREBO2()
{
  if (f90obj) {
    if (screened) {
      rebo2_scr_del(f90obj);
      rebo2_scr_destroy(f90obj);
    }
    else {
      rebo2_del(f90obj);
      rebo2_destroy(f90obj);
    }
  }

  memory->sfree(REBO_seed);
  memory->sfree(REBO_last);
  memory->sfree(REBO_neighb);
  memory->sfree(REBO_el);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cutghost);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairREBO2::compute(int eflag, int vflag)
{
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = vflag_atom = 0;

  REBO_neigh();
  FREBO(eflag,vflag);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairREBO2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cutghost,n+1,n+1,"pair:cutghost");

  // only sized by C,H = 2 types

  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairREBO2::settings(int narg, char **arg)
{
  if (narg != 0 && narg != 1) error->all(FLERR,"Illegal pair_style command");

  if (narg == 1)
    screened = force->inumeric(arg[0]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairREBO2::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg != 2 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients -> narg != 2 + "
	       "atom->ntypes");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients -> * *");

  // read args that map atom types to C and H
  // map[i] = which element (0,1) the Ith atom type is, -1 if NULL

  for (int i = 2; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-1] = -1;
      continue;
    } else if (strcmp(arg[i],"C") == 0) {
      map[i-1] = 0;
    } else if (strcmp(arg[i],"H") == 0) {
      map[i-1] = 1;
    } else error->all(FLERR,"Incorrect args for pair coefficients -> "
		      "undefined atom species");
  }

  // clear setflag since coeff() called once with I,J = * *

  int n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
	setflag[i][j] = 1;
	count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients -> "
			     "count = 0");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairREBO2::init_style()
{
  int i, j;
  double rc;

  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style REBO2 requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style REBO2 requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
  neighbor->requests[irequest]->ghost = 1;

  char name[16];

  // initialize f90 rebo2 object
  if (screened) {
    strcpy(name, "REBO2+S");
    rebo2_scr_create(&f90obj);
    rebo2_scr_init(f90obj);
    rebo2_scr_get_cutoff(f90obj, rcmaxsq);
  }
  else {
    strcpy(name, "REBO2");
    rebo2_create(&f90obj);
    rebo2_init(f90obj);
    rebo2_get_cutoff(f90obj, rcmaxsq);
  }

  rc = 0.0;
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      rcmax[i][j] = sqrt(rcmaxsq[i][j]);
      rc = MAX(rc, rcmax[i][j]);
    }
  }

  // determine ghost communication
  if (screened) {
    comm->cutghostuser = MAX(comm->cutghostuser, 4*rc);
  }
  else {
    comm->cutghostuser = MAX(comm->cutghostuser, 3*rc);
  }

  // dump information about cuoffs
  if (screen) {
    fprintf(screen, "%s, cut-offs: ", name);
    fprintf(screen, "C-C = %f, ", rcmax[0][0]);
    fprintf(screen, "H-H = %f, ", rcmax[1][1]);
    fprintf(screen, "C-H = %f, %f\n", rcmax[0][1], rcmax[1][0]);
  }
  if (logfile) {
    fprintf(logfile, "%s, cut-offs: ", name);
    fprintf(logfile, "C-C = %f, ", rcmax[0][0]);
    fprintf(logfile, "H-H = %f, ", rcmax[1][1]);
    fprintf(logfile, "C-H = %f, %f\n", rcmax[0][1], rcmax[1][0]);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairREBO2::init_one(int i, int j)
{
  //printf("--------> ENTERING INIT_ONE by Indizes %d,%d\n", i,j);
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double rc = rcmax[map[i]][map[j]];

  cutghost[i][j] = cutghost[j][i] = rc;

  return rc;
}

/* ----------------------------------------------------------------------
   create REBO neighbor list from main neighbor list
   REBO neighbor list stores neighbors of ghost atoms
------------------------------------------------------------------------- */

void PairREBO2::REBO_neigh()
{
  //printf("...entering REBO_neigh():\n");
  int i,j,ii,jj,n,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;

  if (!list->ghostflag) {
    error->all(FLERR,"REBO2 needs neighbor list with ghost atoms.");
  }

  if (nall > maxlocal) {
    maxlocal = atom->nmax;
    memory->sfree(REBO_seed);
    memory->sfree(REBO_last);
    memory->sfree(REBO_neighb);
    memory->sfree(REBO_el);
    REBO_seed = (int *)
      memory->smalloc(maxlocal*sizeof(int),"REBO2:REBO_seed");
    REBO_last = (int *)
      memory->smalloc(maxlocal*sizeof(int),"REBO2:REBO_last");
    REBO_nneighb = 50*maxlocal;
    REBO_neighb = (int *)
      memory->smalloc(REBO_nneighb*sizeof(int),"REBO2:REBO_neighb");      
    REBO_el = (int *)
      memory->smalloc(maxlocal*sizeof(int),"REBO2:REBO_el");
  }

  // set start values for neighbor array REBO_neighb
  for (i = 0; i < nall; i++) {
    REBO_seed[i] = -1;
    REBO_last[i] = -2;

    // write REBO2 array of atom types
    REBO_el[i] =  map[type[i]] + 1;
  }

  inum = list->inum+list->gnum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // copy LAMMPS style neighbor list into MDCORE style neighbor list
  // FIXME: this should be harmonized at some point
  n = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    REBO_seed[i] = n+1;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = map[type[i]];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      jtype = map[type[j]];
      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < rcmaxsq[itype][jtype]) {
	REBO_neighb[n] = j+1;
	n++;
      }
    }
    REBO_last[i] = n;

    if (n >= REBO_nneighb)
      error->one(FLERR,"Neighbor list overflow. Increase REBO_nneighb.");

  } // end for loop over all ii (inum)

#if 0
  // DEBUG: Check if neighbor list is symmetric
  for (i = 0; i < nall; i++) {
    for (ii = REBO_seed[i]-1; ii < REBO_last[i]; ii++) {
      j = REBO_neighb[ii]-1;

      // Check if i is neighbor of j
      n = 0;
      for (jj = REBO_seed[j]-1; jj < REBO_last[j]; jj++) {
	if (REBO_neighb[jj]-1 == i) {
	  n = 1;
	}
      }
      if (!n) {
	printf("i = %i, j = %i\n", i, j);
	printf("Neighbors of i\n");
	for (jj = REBO_seed[i]-1; jj < REBO_last[i]; jj++) {
	  printf("   %i\n", REBO_neighb[jj]-1);
	}
	printf("Neighbors of j\n");
	for (jj = REBO_seed[j]-1; jj < REBO_last[j]; jj++) {
	  printf("   %i\n", REBO_neighb[jj]-1);
	}
	error->one(FLERR,"Neighbor list not symmetric");
      }
    }
  }
#endif
}

/* ----------------------------------------------------------------------
   REBO forces and energy
------------------------------------------------------------------------- */

void PairREBO2::FREBO(int eflag, int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = atom->nlocal + atom->nghost;
  double *wpot_per_at, wpot[3][3];

  memset(wpot, 0, 9*sizeof(double));

  // Convert LAMMPS atom types to REBO2 array of atom types
  for (int i = 0; i < nall; i++) {
    REBO_el[i] = map[type[i]] + 1;
  }

  wpot_per_at = NULL;
  if (vflag_atom) {
    wpot_per_at = vatom[0];
  }

  if (screened) {
    rebo2_scr_energy_and_forces(f90obj,
				nlocal, nall, tag, REBO_el, x[0],
				REBO_seed, REBO_last, REBO_nneighb,
				REBO_neighb, &eng_vdwl, f[0], wpot, wpot_per_at,
				NULL);
  }
  else {
    rebo2_energy_and_forces(f90obj,
			    nlocal, nall, tag, REBO_el, x[0],
			    REBO_seed, REBO_last, REBO_nneighb,
			    REBO_neighb, &eng_vdwl, f[0], wpot, wpot_per_at,
			    NULL);
  }

  // update virial
  virial[0] += wpot[0][0];
  virial[1] += wpot[1][1];
  virial[2] += wpot[2][2];
  virial[3] += wpot[1][0];
  virial[4] += wpot[2][0];
  virial[5] += wpot[2][1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays 
------------------------------------------------------------------------- */

double PairREBO2::memory_usage()
{
  double bytes = 0.0;
  return bytes;
}

