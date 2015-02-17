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
         CSTR("Juslin-Type bond-order potential (screened)."))
#else
    m = ptrdict_register_section(cfg, CSTR(BOP_STR), &
         CSTR("Juslin-Type bond-order potential."))
#endif

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, JUSLIN_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_string_property(m, c_loc(this%ref(1:1)), &
         JUSLIN_MAX_REF, CSTR("ref"), &
         CSTR("Reference string to choose a parameters set from the database."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%D0), JUSLIN_MAX_PAIRS, c_loc(this%db%nD0), &
         CSTR("D0"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r0), JUSLIN_MAX_PAIRS, c_loc(this%db%nr0), &
         CSTR("r0"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%S), JUSLIN_MAX_PAIRS, c_loc(this%db%nS), &
         CSTR("S"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%beta),JUSLIN_MAX_PAIRS,c_loc(this%db%nbeta), &
         CSTR("beta"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%gamma), JUSLIN_MAX_PAIRS, &
         c_loc(this%db%ngamma), &
         CSTR("gamma"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%c), JUSLIN_MAX_PAIRS, c_loc(this%db%nc), &
         CSTR("c"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%d), JUSLIN_MAX_PAIRS, c_loc(this%db%nd), &
         CSTR("d"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%h), JUSLIN_MAX_PAIRS, c_loc(this%db%nh), &
         CSTR("h"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%n), JUSLIN_MAX_PAIRS, c_loc(this%db%nn), &
         CSTR("n"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%alpha), JUSLIN_MAX_EL**3, c_loc(this%db%nalpha), &
         CSTR("alpha"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%omega), JUSLIN_MAX_EL**3, c_loc(this%db%nomega), &
         CSTR("omega"), CSTR("See functional form."))
    call ptrdict_register_integer_list_property(m, &
         c_loc1(this%db%m), JUSLIN_MAX_EL**3, c_loc(this%db%nm), &
         CSTR("m"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r1), JUSLIN_MAX_PAIRS, c_loc(this%db%nr1), &
         CSTR("r1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%r2), JUSLIN_MAX_PAIRS, c_loc(this%db%nr2), &
         CSTR("r2"), CSTR("See functional form."))
#ifdef SCREENING
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or1), JUSLIN_MAX_PAIRS, c_loc(this%db%nor1), &
         CSTR("or1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%or2), JUSLIN_MAX_PAIRS, c_loc(this%db%nor2), &
         CSTR("or2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor1),JUSLIN_MAX_PAIRS,c_loc(this%db%nbor1), &
         CSTR("bor1"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%bor2),JUSLIN_MAX_PAIRS,c_loc(this%db%nbor2), &
         CSTR("bor2"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmin),JUSLIN_MAX_PAIRS,c_loc(this%db%nCmin), &
         CSTR("Cmin"), CSTR("See functional form."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Cmax),JUSLIN_MAX_PAIRS,c_loc(this%db%nCmax), &
         CSTR("Cmax"), CSTR("See functional form."))
#endif

  endsubroutine REGISTER_FUNC
