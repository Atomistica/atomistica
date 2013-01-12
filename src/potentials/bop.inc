#ifndef __BOP_INC
#define __BOP_INC

#include "have.inc"

#ifdef HAVE_ALBE_BOP
  use albe_bop
#endif
#ifdef HAVE_ALBE_BOP_SCR
  use albe_bop_scr
#endif
#ifdef HAVE_JUSLIN_BOP
  use juslin_bop
#endif
#ifdef HAVE_JUSLIN_BOP_SCR
  use juslin_bop_scr
#endif
#ifdef HAVE_BORN_MAYER
  use born_mayer
#endif
#ifdef HAVE_EXP_BOP
  use exp_bop
#endif
#ifdef HAVE_EXP_BOP_SCR
  use exp_bop_scr
#endif
#ifdef HAVE_EXP_BOP_ATT
  use exp_bop_att
#endif
#ifdef HAVE_KUMAGAI
  use kumagai
#endif
#ifdef HAVE_KUMAGAI_SCR
  use kumagai_scr
#endif
#ifdef HAVE_O2_BOP
  use o2_bop
#endif
#ifdef HAVE_O2_BOP_SCR
  use o2_bop_scr
#endif
#ifdef HAVE_REBO2
  use rebo2
#endif
#ifdef HAVE_REBO2_SCR
  use rebo2_scr
#endif
#ifdef HAVE_REBO2X
  use rebo2x
#endif
#ifdef HAVE_REBO2X_SCR
  use rebo2x_scr
#endif
#ifdef HAVE_REBO2_WITHOUT_TABLES
  use rebo2_without_tables
#endif
#ifdef HAVE_REBO2_WITHOUT_TABLES_SCR
  use rebo2_without_tables_scr
#endif
#ifdef HAVE_TERSOFF
  use tersoff
#endif
#ifdef HAVE_TERSOFF_SCR
  use tersoff_scr
#endif

#endif
