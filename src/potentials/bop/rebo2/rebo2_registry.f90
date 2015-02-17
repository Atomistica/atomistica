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

!    this = default_parameters

    call rebo2_default_Fcc_table(this%in_Fcc, this%in_dFdi, this%in_dFdj, this%in_dFdk)
    call rebo2_default_Fch_table(this%in_Fch)
    call rebo2_default_Fhh_table(this%in_Fhh)
    call rebo2_default_Pcc_table(this%in_Pcc)
    call rebo2_default_Pch_table(this%in_Pch)
    call rebo2_default_Tcc_table(this%in_Tcc)

#ifdef SCREENING
    m = ptrdict_register_section(cfg, CSTR(BOP_STR), &
         CSTR("The screened 2nd generation REBO (Brenner 2002) potential."))
#else
    m = ptrdict_register_section(cfg, CSTR(BOP_STR), &
         CSTR("The 2nd generation REBO (Brenner 2002) potential."))
#endif

    call ptrdict_register_string_property(m, c_locs(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements for which to use this potential (default: C,H)."))

!! ================= Added Some Parameters H-H =================

    call ptrdict_register_real_property(m, c_loc(this%hh_Q), &
         CSTR("HH_Q"), CSTR("Q for H-H interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%hh_A), &
         CSTR("HH_A"), CSTR("A for H-H interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%hh_alpha), &
         CSTR("HH_alpha"), CSTR("alpha for H-H interaction (inner)."))

    call ptrdict_register_real_property(m, c_loc(this%hh_B1), &
         CSTR("HH_B1"), CSTR("B1 for H-H interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%hh_beta1), &
         CSTR("HH_beta1"), CSTR("beta1 for H-H interaction (inner)."))

!! ================= C-H Parameters ====================

    call ptrdict_register_real_property(m, c_loc(this%ch_Q), &
         CSTR("CH_Q"), CSTR("Q for C-H interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%ch_A), &
         CSTR("CH_A"), CSTR("A for C-H interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%ch_alpha), &
         CSTR("CH_alpha"), CSTR("alpha for C-H interaction (inner)."))

    call ptrdict_register_real_property(m, c_loc(this%ch_B1), &
         CSTR("CH_B1"), CSTR("B1 for C-H interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%ch_beta1), &
         CSTR("CH_beta1"), CSTR("beta1 for C-H interaction (inner)."))

!! ================= C-C Parameters ====================

    call ptrdict_register_real_property(m, c_loc(this%cc_Q), &
         CSTR("CC_Q"), CSTR("Q for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_A), &
         CSTR("CC_A"), CSTR("A for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_alpha), &
         CSTR("CC_alpha"), CSTR("alpha for C-C interaction (inner)."))

    call ptrdict_register_real_property(m, c_loc(this%cc_B1), &
         CSTR("CC_B1"), CSTR("B1 for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_B2), &
         CSTR("CC_B2"), CSTR("B2 for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_B3), &
         CSTR("CC_B3"), CSTR("B3 for C-C interaction (inner)."))

    call ptrdict_register_real_property(m, c_loc(this%cc_beta1), &
         CSTR("CC_beta1"), CSTR("beta1 for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_beta2), &
         CSTR("CC_beta2"), CSTR("beta2 for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_beta3), &
         CSTR("CC_beta3"), CSTR("beta3 for C-C interaction (inner)."))

#ifdef SCREENING
    call ptrdict_register_real_property(m, c_loc(this%Cmin), CSTR("Cmin"), &
         CSTR("Lower screening cut-off (should be >= 1)."))
    call ptrdict_register_real_property(m, c_loc(this%Cmax), CSTR("Cmax"), &
         CSTR("Upper screening cut-off (should be <= 3)."))
#endif

    call ptrdict_register_real_property(m, c_loc(this%cc_in_r1), &
         CSTR("CC_in_r1"), CSTR("r1 for C-C interaction (inner)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_in_r2), &
         CSTR("CC_in_r2"), CSTR("r2 for C-C interaction (inner)."))
#ifdef SCREENING
    call ptrdict_register_real_property(m, c_loc(this%cc_ar_r1), &
         CSTR("CC_ar_r1"), &
         CSTR("r1 for C-C interaction (attractive-repulsive)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_ar_r2), &
         CSTR("CC_ar_r2"), &
         CSTR("r2 for C-C interaction (attractive-repulsive)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_bo_r1), &
         CSTR("CC_bo_r1"), CSTR("r1 for C-C interaction (bond-order)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_bo_r2), &
         CSTR("CC_bo_r2"), CSTR("r2 for C-C interaction (bond-order)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_nc_r1), &
         CSTR("CC_nc_r1"), &
         CSTR("r1 for C-C interaction (neighbor-conjugation)."))
    call ptrdict_register_real_property(m, c_loc(this%cc_nc_r2), &
         CSTR("CC_nc_r2"), &
         CSTR("r2 for C-C interaction (neighbor-conjugation)."))
#endif
    call ptrdict_register_real_property(m, c_loc(this%ch_r1), CSTR("CH_r1"), &
         CSTR("r1 for C-H interaction."))
    call ptrdict_register_real_property(m, c_loc(this%ch_r2), CSTR("CH_r2"), &
         CSTR("r2 for C-H interaction."))
    call ptrdict_register_real_property(m, c_loc(this%hh_r1), CSTR("HH_r1"), &
         CSTR("r1 for H-H interaction."))
    call ptrdict_register_real_property(m, c_loc(this%hh_r2), CSTR("HH_r2"), &
         CSTR("r2 for H-H interaction."))

    call ptrdict_register_boolean_property(m, c_loc(this%with_dihedral), &
         CSTR("dihedral"), CSTR("Include the dihedral term?"))

    call ptrdict_register_boolean_property(m, c_loc(this%zero_tables), &
         CSTR("zero_tables"), CSTR("Initialize all tables to zero."))

    call ptrdict_register_array3d_property(m, c_loc111(this%in_Fcc), 5, 5, 10, &
         CSTR("Fcc"), CSTR("Fcc-table"))
    call ptrdict_register_array3d_property(m, c_loc111(this%in_dFdi), 5, 5, 10, &
         CSTR("dFdi"), CSTR("dFdi"))
    call ptrdict_register_array3d_property(m, c_loc111(this%in_dFdj), 5, 5, 10, &
         CSTR("dFdj"), CSTR("dFdj"))
    call ptrdict_register_array3d_property(m, c_loc111(this%in_dFdk), 5, 5, 10, &
         CSTR("dFdk"), CSTR("dFdk"))
    call ptrdict_register_array3d_property(m, c_loc111(this%in_Fch), 5, 5, 10, &
         CSTR("Fch"), CSTR("Fch-table"))
    call ptrdict_register_array3d_property(m, c_loc111(this%in_Fhh), 5, 5, 10, &
         CSTR("Fhh"), CSTR("Fhh-table"))
    call ptrdict_register_array2d_property(m, c_loc11(this%in_Pcc), 6, 6, &
         CSTR("Pcc"), CSTR("Pcc-table"))
    call ptrdict_register_array2d_property(m, c_loc11(this%in_Pch), 6, 6, &
         CSTR("Pch"), CSTR("Pch-table"))
    call ptrdict_register_array3d_property(m, c_loc111(this%in_Tcc), 5, 5, 10, &
         CSTR("Tcc"), CSTR("Tcc-table"))

  endsubroutine REGISTER_FUNC
