#ifndef __BOP_INC
#define __BOP_INC

#include "have.inc"

#ifdef HAVE_BRENNER
  use Brenner
#endif
#ifdef HAVE_BRENNER_SCR
  use Brenner_scr
#endif
#ifdef HAVE_KUMAGAI
  use kumagai
#endif
#ifdef HAVE_KUMAGAI_SCR
  use kumagai_scr
#endif
#ifdef HAVE_TERSOFF
  use tersoff
#endif
#ifdef HAVE_TERSOFF_SCR
  use tersoff_scr
#endif

#endif
