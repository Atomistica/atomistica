!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 2 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
         CSTR("The Tersoff potential (screened). Parameters are named according Tersoff, PRB 39, 5566 (1989)."))
#else
    m = ptrdict_register_section(cfg, CSTR(BOP_STR), &
         CSTR("The Tersoff potential. Parameters are named according Tersoff, PRB 39, 5566 (1989)."))
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
         c_loc1(this%db%lambda), TERSOFF_MAX_PAIRS, c_loc(this%db%nlambda), &
         CSTR("lambda"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%mu), TERSOFF_MAX_PAIRS, c_loc(this%db%nmu), &
         CSTR("mu"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%omega), TERSOFF_MAX_PAIRS, c_loc(this%db%nomega), &
         CSTR("omega"), CSTR("See functional form."))
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

  endsubroutine REGISTER_FUNC
