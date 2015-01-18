#ifndef __STANDALONE_INC
#define __STANDALONE_INC

#include "have.inc"

#ifdef HAVE_CONSTANT_FORCE
  use constant_force
#endif
#ifdef HAVE_CONFINEMENT
  use confinement
#endif
#ifdef HAVE_HARMONIC_HOOK
  use harmonic_hook
#endif
#ifdef HAVE_FFM_TIP
  use ffm_tip
#endif

#ifdef HAVE_SLIDING_P
  use sliding_p
#endif

#ifdef MDCORE_MONOLITHIC
  use pair_potential
#endif

#ifdef HAVE_LENNARD_JONES
  use lennard_jones
#endif
#ifdef HAVE_HARMONIC
  use harmonic
#endif

#ifdef MDCORE_MONOLITHIC
  use potentials_standalone
#endif

#endif
