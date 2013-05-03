!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================

  subroutine REGISTER_FUNC(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(BOP_TYPE), target      :: this
    type(c_ptr),    intent(in)  :: cfg
    type(c_ptr),    intent(out) :: m

    ! ---

#ifdef SCREENING
    m = ptrdict_register_section(cfg, CSTR(BOP_STR), &
         CSTR("The Tersoff potential (screened)."))
#else
    m = ptrdict_register_section(cfg, CSTR(BOP_STR), &
         CSTR("The Tersoff potential."))
#endif

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, TERSOFF_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

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

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), TERSOFF_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), TERSOFF_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))
#ifdef SCREENING
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or1), TERSOFF_MAX_PAIRS, c_loc(this%db%nor1), &
         CSTR("or1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or2), TERSOFF_MAX_PAIRS, c_loc(this%db%nor2), &
         CSTR("or2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor1), TERSOFF_MAX_PAIRS, c_loc(this%db%nbor1), &
         CSTR("bor1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor2), TERSOFF_MAX_PAIRS, c_loc(this%db%nbor2), &
         CSTR("bor2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmin), TERSOFF_MAX_PAIRS, c_loc(this%db%nCmin), &
         CSTR("Cmin"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmax), TERSOFF_MAX_PAIRS, c_loc(this%db%nCmax), &
         CSTR("Cmax"), CSTR("See functional form."))
#endif

    call ptrdict_register_integer_property(m, c_loc(this%nebmax), &
         CSTR("nebmax"), CSTR("Internal per-loop neighbor list size."))
    call ptrdict_register_integer_property(m, c_loc(this%nebavg), &
         CSTR("nebavg"), CSTR("Internal global neighbor list size."))

  endsubroutine REGISTER_FUNC
