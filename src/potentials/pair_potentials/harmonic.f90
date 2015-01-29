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

! @meta
!   shared
!   classtype:harmonic_t classname:Harmonic interface:potentials
!   features:per_at
! @endmeta

!>
!! Harmonic interaction potential (i.e. springs)
!<

#include "macros.inc"
#include "filter.inc"

module harmonic
  use libAtoms_module

  use ptrdict

  use logging
  use timer

  use neighbors
  use particles
  use filter

  implicit none

  private

  public :: harmonic_t
  type harmonic_t

     !
     ! Element on which to apply the force
     !

     character(MAX_EL_STR) :: element1 = "*"
     character(MAX_EL_STR) :: element2 = "*"
     integer               :: el1
     integer               :: el2

     !
     ! constants
     !
     
     real(DP)      :: k = 1.0_DP
     real(DP)      :: r0 = 1.0_DP
     real(DP)      :: cutoff = 1.5_DP
     logical(BOOL) :: shift = .false.

     !
     ! derived parameters
     !

     real(DP) :: offset

  endtype harmonic_t


  public :: bind_to
  interface bind_to
     module procedure harmonic_bind_to
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure harmonic_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure harmonic_register
  endinterface

contains

  !>
  !! Initialization
  !!
  !! Initialization
  !<
  subroutine harmonic_bind_to(this, p, nl, ierror)
    implicit none

    type(harmonic_t),  intent(inout) :: this
    type(particles_t), intent(in)    :: p
    type(neighbors_t), intent(inout) :: nl
    integer, optional, intent(inout) :: ierror

    ! ---

    integer :: i, j, k

    ! ---

    this%el1 = filter_from_string(this%element1, p)
    this%el2 = filter_from_string(this%element2, p)

    write (ilog, '(A)')            "- harmonic_init -"
    call filter_prlog(this%el1, p, indent=5)
    call filter_prlog(this%el2, p, indent=5)
    write (ilog, '(5X,A,F20.10)')  "k       = ", this%k
    write (ilog, '(5X,A,F20.10)')  "r0      = ", this%r0
    write (ilog, '(5X,A,F20.10)')  "cutoff  = ", this%cutoff
    write (ilog, '(5X,A,L)')       "shift   = ", this%shift

    do i = 1, p%nel
       do j = 1, p%nel
          if (IS_EL2(this%el1, i) .and. IS_EL2(this%el2, j)) then
             call request_interaction_range(nl, this%cutoff, i, j)
          endif
       enddo
    enddo

    this%offset = 0.0_DP
    if (this%shift) then
       this%offset = 0.5_DP*this%k*(this%cutoff-this%r0)**2
    endif

    write (ilog, *)

  endsubroutine harmonic_bind_to


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine harmonic_energy_and_forces(this, p, nl, epot, f, wpot, &
       epot_per_at, wpot_per_at, ierror)
    implicit none

    type(harmonic_t),   intent(inout) :: this
    type(particles_t),  intent(in)    :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(inout) :: epot
    real(DP),           intent(inout) :: f(3, p%maxnatloc)  !< forces
    real(DP),           intent(inout) :: wpot(3, 3)
    real(DP), optional, intent(inout) :: epot_per_at(p%maxnatloc)
#ifdef LAMMPS
    real(DP), optional, intent(inout) :: wpot_per_at(6, p%maxnatloc)
#else
    real(DP), optional, intent(inout) :: wpot_per_at(3, 3, p%maxnatloc)
#endif
    integer,  optional, intent(inout) :: ierror

    ! ---

    integer             :: i, j
    integer(NEIGHPTR_T) :: jn

    real(DP)  :: dr(3), df(3), dw(3, 3)
    real(DP)  :: cut_sq, abs_dr, for, en, fac12, fac6

    ! ---

    call timer_start("harmonic_force")

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    cut_sq = this%cutoff**2

    do i = 1, p%nat
       do jn = nl%seed(i), nl%last(i)
          j = GET_NEIGHBOR(nl, jn)

          if (i > j) then
             if ( ( IS_EL(this%el1, p, i) .and. IS_EL(this%el2, p, j) ) .or. &
                  ( IS_EL(this%el2, p, i) .and. IS_EL(this%el1, p, j) ) ) then

                DIST_SQ(p, nl, i, jn, dr, abs_dr)

                if (abs_dr < cut_sq) then
                   abs_dr = sqrt(abs_dr)

                   for = this%k*(this%r0-abs_dr)
                   en  = 0.5_DP*for*(this%r0-abs_dr) - this%offset

                   epot = epot + en
                   df   = for * dr/abs_dr

                   VEC3(f, i) = VEC3(f, i) + df
                   VEC3(f, j) = VEC3(f, j) - df

                   dw    = -outer_product(dr, df)
                   wpot  = wpot + dw

                   if (present(epot_per_at)) then
                      en = en/2
                      epot_per_at(i) = epot_per_at(i) + en
                      epot_per_at(j) = epot_per_at(j) + en
                   endif

                   if (present(wpot_per_at)) then
                      dw = dw/2
                      SUM_VIRIAL(wpot_per_at, i, dw)
                      SUM_VIRIAL(wpot_per_at, j, dw)
                   endif

                endif
             endif
          endif
       enddo
    enddo

    call timer_stop("harmonic_force")

  endsubroutine harmonic_energy_and_forces


  subroutine harmonic_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(harmonic_t), target      :: this
    type(c_ptr),      intent(in)  :: cfg
    type(c_ptr),      intent(out) :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("Harmonic"), &
         CSTR("Harmonic springs"))

    call ptrdict_register_string_property(m, c_locs(this%element1), &
         MAX_EL_STR, CSTR("el1"), CSTR("First element."))
    call ptrdict_register_string_property(m, c_locs(this%element2), &
         MAX_EL_STR, CSTR("el2"), CSTR("Second element."))

    call ptrdict_register_real_property(m, c_loc(this%k), &
         CSTR("k"), CSTR("Spring constant."))
    call ptrdict_register_real_property(m, c_loc(this%r0), &
         CSTR("r0"), CSTR("Equilibrium length."))
    call ptrdict_register_real_property(m, c_loc(this%cutoff), CSTR("cutoff"), &
         CSTR("Cutoff length."))
    call ptrdict_register_boolean_property(m, c_loc(this%shift), &
         CSTR("shift"), CSTR("Shift potential to zero energy at cutoff."))

  endsubroutine harmonic_register

endmodule harmonic
