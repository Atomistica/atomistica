#ifndef __COULOMB_INC
#define __COULOMB_INC

#include "have.inc"

#ifdef HAVE_CHARGE_OVERLAP
  use charge_overlap
#endif
#ifdef HAVE_CUTOFF_COULOMB
  use cutoff_coulomb
#endif
#ifdef HAVE_DIRECT_COULOMB
  use direct_coulomb
#endif

#endif
