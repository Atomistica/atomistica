!! ======================================================================
!! MDCORE - Interatomic potential library
!! https://github.com/pastewka/mdcore
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================

! @meta
!   shared
! @endmeta

#include "macros.inc"

#define c_loc1(x) c_loc(x(lbound(x,1)))
#define c_loc11(x) c_loc(x(lbound(x,1),lbound(x,2)))
#define c_loc111(x) c_loc(x(lbound(x,1),lbound(x,2),lbound(x,3)))

!>
!! Registry
!!
!! Contains function to register each classes property
!<
module bop_registry
  use filter, only: MAX_EL_STR

#include "bop.inc"

  use ptrdict

  implicit none

  private

  ! Potentials
#ifdef HAVE_BRENNER
  public :: brenner_register
#endif
#ifdef HAVE_BRENNER_SCR
  public :: brenner_scr_register
#endif
#ifdef HAVE_KUMAGAI
  public :: kumagai_register
#endif
#ifdef HAVE_KUMAGAI_SCR
  public :: kumagai_scr_register
#endif
#ifdef HAVE_TERSOFF
  public :: tersoff_register
#endif
#ifdef HAVE_TERSOFF_SCR
  public :: tersoff_scr_register
#endif

contains

  !
  ! Bond-order potentials
  !

#ifdef HAVE_BRENNER
  subroutine brenner_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(brenner_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)          :: cfg
    type(c_ptr), intent(out)         :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Brenner"), &
         CSTR("Tersoff-Brenner type bond-order potential."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, BRENNER_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_string_property(m, c_loc(this%ref(1)), &
         BRENNER_MAX_REF, &
         CSTR("ref"), &
         CSTR("Reference string to choose a parameters set from the database."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%D0), BRENNER_MAX_PAIRS, c_loc(this%db%nD0), &
         CSTR("D0"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r0), BRENNER_MAX_PAIRS, c_loc(this%db%nr0), &
         CSTR("r0"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%S), BRENNER_MAX_PAIRS, c_loc(this%db%nS), &
         CSTR("S"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%beta), BRENNER_MAX_PAIRS, c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%gamma), BRENNER_MAX_PAIRS, c_loc(this%db%ngamma), &
         CSTR("gamma"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c), BRENNER_MAX_PAIRS, c_loc(this%db%nc), &
         CSTR("c"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%d), BRENNER_MAX_PAIRS, c_loc(this%db%nd), &
         CSTR("d"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), BRENNER_MAX_PAIRS, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mu), BRENNER_MAX_PAIRS, c_loc(this%db%nmu), &
         CSTR("mu"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%n), BRENNER_MAX_PAIRS, c_loc(this%db%nn), &
         CSTR("n"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%m), BRENNER_MAX_PAIRS, c_loc(this%db%nm), &
         CSTR("m"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), BRENNER_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), BRENNER_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))

    call ptrdict_register_integer_property(m, c_loc(this%nebmax), &
         CSTR("nebmax"), CSTR("Maximum number of neighbors (internal neighbor list)."))
    call ptrdict_register_integer_property(m, c_loc(this%nebavg), &
         CSTR("nebavg"), CSTR("Average number of neighbors (internal neighbor list)."))

  endsubroutine brenner_register
#endif

#ifdef HAVE_BRENNER_SCR
  subroutine brenner_scr_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(brenner_scr_t), target, intent(inout)  :: this
    type(c_ptr),                  intent(in)     :: cfg
    type(c_ptr),                  intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("BrennerScr"), &
         CSTR("Tersoff-Brenner type bond-order potential (screened)."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, BRENNER_SCR_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_string_property(m, c_loc(this%ref(1)), &
         BRENNER_MAX_REF, &
         CSTR("ref"), &
         CSTR("Reference string to choose a parameters set from the database."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%D0), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nD0), &
         CSTR("D0"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r0), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nr0), &
         CSTR("r0"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%S), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nS), &
         CSTR("S"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%beta), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%gamma),BRENNER_SCR_MAX_PAIRS,c_loc(this%db%ngamma), &
         CSTR("gamma"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nc), &
         CSTR("c"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%d), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nd), &
         CSTR("d"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mu), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nmu), &
         CSTR("mu"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%n), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nn), &
         CSTR("n"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%m), BRENNER_MAX_PAIRS, c_loc(this%db%nm), &
         CSTR("m"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or1), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nor1), &
         CSTR("or1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or2), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nor2), &
         CSTR("or2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor1), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nbor1), &
         CSTR("bor1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor2), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nbor2), &
         CSTR("bor2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmin), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nCmin), &
         CSTR("Cmin"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmax), BRENNER_SCR_MAX_PAIRS, c_loc(this%db%nCmax), &
         CSTR("Cmax"), CSTR("See functional form."))

    call ptrdict_register_integer_property(m, c_loc(this%nebmax), &
         CSTR("nebmax"), CSTR("Maximum number of neighbors (internal neighbor list)."))
    call ptrdict_register_integer_property(m, c_loc(this%nebavg), &
         CSTR("nebavg"), CSTR("Average number of neighbors (internal neighbor list)."))

  endsubroutine brenner_scr_register
#endif

#ifdef HAVE_KUMAGAI
  subroutine kumagai_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(kumagai_t), target, intent(inout)  :: this
    type(c_ptr),             intent(in)     :: cfg
    type(c_ptr),             intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Kumagai"), &
         CSTR("The Kumagai-Izumi-Hara-Sakai potential."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, KUMAGAI_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%A), KUMAGAI_MAX_PAIRS, c_loc(this%db%nA), &
         CSTR("A"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%B), KUMAGAI_MAX_PAIRS, c_loc(this%db%nB), &
         CSTR("B"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%lambda1), KUMAGAI_MAX_PAIRS, &
         c_loc(this%db%nlambda1), CSTR("lambda1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%lambda2), KUMAGAI_MAX_PAIRS, &
         c_loc(this%db%nlambda2), CSTR("lambda2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%eta), KUMAGAI_MAX_EL, c_loc(this%db%neta), &
         CSTR("eta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%delta), KUMAGAI_MAX_EL, c_loc(this%db%ndelta), &
         CSTR("delta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%alpha), KUMAGAI_MAX_PAIRS, c_loc(this%db%nalpha), &
         CSTR("alpha"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%beta), KUMAGAI_MAX_PAIRS, c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c1), KUMAGAI_MAX_EL, c_loc(this%db%nc1), &
         CSTR("c1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c2), KUMAGAI_MAX_EL, c_loc(this%db%nc2), &
         CSTR("c2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c3), KUMAGAI_MAX_EL, c_loc(this%db%nc3), &
         CSTR("c3"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c4), KUMAGAI_MAX_EL, c_loc(this%db%nc4), &
         CSTR("c4"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c5), KUMAGAI_MAX_EL, c_loc(this%db%nc5), &
         CSTR("c5"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), KUMAGAI_MAX_EL, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), KUMAGAI_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), KUMAGAI_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))

  endsubroutine kumagai_register
#endif

#ifdef HAVE_KUMAGAI_SCR
  subroutine kumagai_scr_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(kumagai_scr_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("KumagaiScr"), &
         CSTR("The Kumagai-Izumi-Hara-Sakai potential (screened)."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, KUMAGAI_SCR_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%A), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nA), &
         CSTR("A"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%B), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nB), &
         CSTR("B"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%lambda1), KUMAGAI_SCR_MAX_PAIRS, &
         c_loc(this%db%nlambda1), CSTR("lambda1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%lambda2), KUMAGAI_SCR_MAX_PAIRS, &
         c_loc(this%db%nlambda2), CSTR("lambda2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%eta), KUMAGAI_SCR_MAX_EL, c_loc(this%db%neta), &
         CSTR("eta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%delta), KUMAGAI_SCR_MAX_EL, c_loc(this%db%ndelta), &
         CSTR("delta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%alpha), KUMAGAI_MAX_PAIRS, c_loc(this%db%nalpha), &
         CSTR("alpha"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%beta), KUMAGAI_MAX_PAIRS, c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c1), KUMAGAI_SCR_MAX_EL, c_loc(this%db%nc1), &
         CSTR("c1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c2), KUMAGAI_SCR_MAX_EL, c_loc(this%db%nc2), &
         CSTR("c2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c3), KUMAGAI_SCR_MAX_EL, c_loc(this%db%nc3), &
         CSTR("c3"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c4), KUMAGAI_SCR_MAX_EL, c_loc(this%db%nc4), &
         CSTR("c4"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c5), KUMAGAI_SCR_MAX_EL, c_loc(this%db%nc5), &
         CSTR("c5"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), KUMAGAI_SCR_MAX_EL, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or1), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nor1), &
         CSTR("or1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or2), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nor2), &
         CSTR("or2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor1), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nbor1), &
         CSTR("bor1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor2), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nbor2), &
         CSTR("bor2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmin), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nCmin), &
         CSTR("Cmin"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmax), KUMAGAI_SCR_MAX_PAIRS, c_loc(this%db%nCmax), &
         CSTR("Cmax"), CSTR("See functional form."))

  endsubroutine kumagai_scr_register
#endif

#ifdef HAVE_TERSOFF
  subroutine tersoff_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(tersoff_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)            :: cfg
    type(c_ptr), intent(out)           :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Tersoff"), &
         CSTR("The Tersoff potential."))

!    call ptrdict_register_string_property(m, this%elements, MAX_EL_STR, CSTR("elements"), &
!         CSTR("Element for which to use this potential."))

!    call ptrdict_register_string_property(m, this%ref, GUPTA_MAX_REF, CSTR("ref"), &
!         CSTR("Reference string to choose a parameters set from the database."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, TERSOFF_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), TERSOFF_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), TERSOFF_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%A), TERSOFF_MAX_PAIRS, c_loc(this%db%nA), &
         CSTR("A"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%B), TERSOFF_MAX_PAIRS, c_loc(this%db%nB), &
         CSTR("B"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%xi), TERSOFF_MAX_PAIRS, c_loc(this%db%nxi), &
         CSTR("xi"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%lambda), TERSOFF_MAX_EL, c_loc(this%db%nlambda), &
         CSTR("lambda"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mu), TERSOFF_MAX_PAIRS, c_loc(this%db%nmu), &
         CSTR("mu"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mubo), TERSOFF_MAX_PAIRS, c_loc(this%db%nmubo), &
         CSTR("mubo"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%m), TERSOFF_MAX_PAIRS, c_loc(this%db%nm), &
         CSTR("m"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%beta), TERSOFF_MAX_EL, c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%n), TERSOFF_MAX_EL, c_loc(this%db%nn), &
         CSTR("n"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c), TERSOFF_MAX_EL, c_loc(this%db%nc), &
         CSTR("c"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%d), TERSOFF_MAX_EL, c_loc(this%db%nd), &
         CSTR("d"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), TERSOFF_MAX_EL, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))

  endsubroutine tersoff_register
#endif

#ifdef HAVE_TERSOFF_SCR
  subroutine tersoff_scr_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(tersoff_scr_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("TersoffScr"), &
         CSTR("The Tersoff potential (screened)."))

!    call ptrdict_register_string_property(m, this%elements, MAX_EL_STR, CSTR("elements"), &
!         CSTR("Element for which to use this potential."))

!    call ptrdict_register_string_property(m, this%ref, GUPTA_MAX_REF, CSTR("ref"), &
!         CSTR("Reference string to choose a parameters set from the database."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, TERSOFF_SCR_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%A), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nA), &
         CSTR("A"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%B), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nB), &
         CSTR("B"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%xi), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nxi), &
         CSTR("xi"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%lambda), TERSOFF_SCR_MAX_EL, c_loc(this%db%nlambda), &
         CSTR("lambda"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mu), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nmu), &
         CSTR("mu"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mubo), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nmubo), &
         CSTR("mubo"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%m), TERSOFF_MAX_PAIRS, c_loc(this%db%nm), &
         CSTR("m"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%beta), TERSOFF_SCR_MAX_EL, c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%n), TERSOFF_SCR_MAX_EL, c_loc(this%db%nn), &
         CSTR("n"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c), TERSOFF_SCR_MAX_EL, c_loc(this%db%nc), &
         CSTR("c"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%d), TERSOFF_SCR_MAX_EL, c_loc(this%db%nd), &
         CSTR("d"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), TERSOFF_SCR_MAX_EL, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or1), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nor1), &
         CSTR("or1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or2), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nor2), &
         CSTR("or2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor1), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nbor1), &
         CSTR("bor1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor2), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nbor2), &
         CSTR("bor2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmin), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nCmin), &
         CSTR("Cmin"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmax), TERSOFF_SCR_MAX_PAIRS, c_loc(this%db%nCmax), &
         CSTR("Cmax"), CSTR("See functional form."))

  endsubroutine tersoff_scr_register
#endif

endmodule bop_registry

