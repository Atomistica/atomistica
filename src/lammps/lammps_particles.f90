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

#include "macros.inc"

!>
!! Particle information
!!
!! Position and cell information are stored in the data structures
!! of this module. This is a LAMMPS compatibility structure that keeps
!! pointers to the LAMMPS data structures.
!<
module particles
  use, intrinsic :: iso_c_binding

  use supplib

  implicit none

  private

  !>
  !! Highest element number stored in the periodic table module
  !<
  integer, parameter :: MAX_Z = ubound(ElementName, 1)
  public :: MAX_Z

  !
  ! This stores the static information, i.e. the *positions*
  !

  public :: particles_t
  type particles_t

     !>
     !! Number of particles in system (including ghost particles)
     !<
     integer                     :: nat = 0

     !>
     !! Number of particles on this processor (excluding ghost particles)
     !<
     integer                     :: natloc = 0

     !>
     !! Length of the position array
     !<
     integer                     :: maxnatloc = 0

     !
     ! All particel data is managed by LAMMPS. The fields below are pointers to
     ! LAMMPS data structures.
     !

     !>
     !! Unique atom id
     !<
     integer(C_INT), pointer     :: tag(:) => NULL()

     !>
     !! Internal element numbers
     !<
     integer(C_INT), pointer     :: el(:) => NULL()

     !>
     !! These positions are always local and may be outside the global box.
     !<
     real(C_DOUBLE), pointer     :: r_non_cyc(:, :) => NULL()

     !>
     !! Mapping of internal element numbers to real elements
     !<
     integer                     :: nel           !> number of distinct elements
     integer                     :: el2Z(MAX_Z)   !> id - i.e. from 1 to nel

     !>
     !! Communication border
     !<
     real(C_DOUBLE)              :: border

     !>
     !! Maximum range of interaction
     !<
     real(C_DOUBLE), allocatable :: interaction_range(:, :)

  endtype particles_t

  public :: request_border
  interface request_border
     module procedure particles_request_border
  endinterface request_border

  public :: set_interaction_range
  interface set_interaction_range
     module procedure particles_set_interaction_range
  endinterface

contains

  !>
  !! Create a new instance
  !<
  subroutine particles_new(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(inout) :: this_cptr

    ! ---

    type(particles_t), pointer :: this_fptr

    ! ---

    allocate(this_fptr)
    this_cptr = c_loc(this_fptr)

  endsubroutine particles_new


  !>
  !! Destroy the instance
  !<
  subroutine particles_free(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr

    ! ---

    type(particles_t), pointer :: this_fptr

    ! ---

    call c_f_pointer(this_cptr, this_fptr)
    deallocate(this_fptr)

  endsubroutine particles_free


  !>
  !! Constructor
  !<
  subroutine particles_init(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr

    ! ---

    type(particles_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    this%nat = 0
    this%natloc = 0
    this%maxnatloc = 0

    this%nel = 0

    this%border = 0.0_DP

    this%tag => NULL()
    this%el => NULL()
    this%r_non_cyc => NULL()

    allocate(this%interaction_range(MAX_Z, MAX_Z))
    this%interaction_range = 0.0_DP

  endsubroutine particles_init


  !>
  !! Destructore
  !<
  subroutine particles_del(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr

    ! ---

    type(particles_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    deallocate(this%interaction_range)

  endsubroutine particles_del


  !>
  !! Associate an internal element number with a chemical element
  !<
  subroutine particles_set_element(this_cptr, el_str_cptr, nel, el_no, Z, error) &
       bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),              value       :: this_cptr
    type(C_PTR),              value       :: el_str_cptr
    integer(C_INT),           value       :: nel
    integer(C_INT),           value       :: el_no
    integer(C_INT),           intent(out) :: Z
#ifdef NO_BIND_C_OPTIONAL
    integer(C_INT),           intent(out) :: error
#else
    integer(C_INT), optional, intent(out) :: error
#endif

    ! ---

    type(particles_t), pointer :: this

    ! ---

#ifdef NO_BIND_C_OPTIONAL
    error = ERROR_NONE
#else
    INIT_ERROR(error)
#endif
    call c_f_pointer(this_cptr, this)

    Z = atomic_number(a2s(c_f_string(el_str_cptr)))
    if (Z <= 0) then
#ifdef NO_BIND_C_OPTIONAL
       call push_error_with_info("Cannot find element '" // a2s(c_f_string(el_str_cptr)) // "'.", __FILE__, __LINE__)
       error = ERROR_UNSPECIFIED
       return
#else
       RAISE_ERROR("Cannot find element '" // a2s(c_f_string(el_str_cptr)) // "'.", error)
#endif
    endif
    this%nel = nel
    this%el2Z(el_no) = Z

  endsubroutine particles_set_element


  !>
  !! Assign pointers to data
  !>
  subroutine particles_set_pointers(this_cptr, nat, natloc, maxnatloc, tag, el, r) &
       bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),    value :: this_cptr
    integer(C_INT), value :: nat, natloc, maxnatloc
    type(C_PTR),    value :: tag, el, r

    ! ---

    type(particles_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    this%nat = nat
    this%natloc = natloc
    this%maxnatloc = maxnatloc
    call c_f_pointer(tag, this%tag, [ maxnatloc ])
    call c_f_pointer(el, this%el, [ maxnatloc ])
    call c_f_pointer(r, this%r_non_cyc, [ 3, maxnatloc ])

  endsubroutine particles_set_pointers


  !>
  !! Request a certain size for the communication border
  !>
  subroutine particles_request_border(this, border)
    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: border

    ! ---

    this%border = max(border, this%border)

  endsubroutine particles_request_border


  !>
  !! Report the true interaction range for this model (which can differ from
  !! the neighbor list cutoff)
  !<
  subroutine particles_set_interaction_range(this, range, el1, el2)
    type(particles_t), intent(inout) :: this
    real(DP),          intent(in)    :: range
    integer, optional, intent(in)    :: el1
    integer, optional, intent(in)    :: el2

    ! ---

    integer :: rel2

    ! ---

    if (present(el1)) then
       if (present(el2)) then
          rel2 = el2
       else
          rel2 = el1
       endif
       this%interaction_range(el1, rel2) = &
          max(this%interaction_range(el1, rel2), range)
       if (el1 /= rel2)  then
          this%interaction_range(rel2, el1) = &
             max(this%interaction_range(rel2, el1), range)
       endif
    else
       this%interaction_range = max(this%interaction_range, range)
    endif

  endsubroutine particles_set_interaction_range


  !>
  !! Return cutoffs for element combination el1, el2
  !<
  subroutine particles_get_interaction_range(this_cptr, el1, el2, range) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),    value       :: this_cptr
    integer(C_INT), value       :: el1, el2
    real(DP),       intent(out) :: range

    ! ---

    type(particles_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    range = this%interaction_range(el1, el2)

  endsubroutine particles_get_interaction_range


  !>
  !! Get the value of the border
  !>
  subroutine particles_get_border(this_cptr, border) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value       :: this_cptr
    real(DP),    intent(out) :: border

    ! ---

    type(particles_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    border = this%border

  endsubroutine particles_get_border


  ! --- Auxiliary stuff ---

  !>
  !! Return error string
  !<
  subroutine get_full_error_string(str) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    character(kind=c_char, len=1), intent(inout)  :: str(*)

    ! ---

    integer          :: i, l
    character(1000)  :: errstr

    ! ---
    
    errstr = get_error_string_and_clear()
    l = len_trim(errstr)
    do i = 1, l
       str(i) = errstr(i:i)
    enddo
    str(l+1) = C_NULL_CHAR

  endsubroutine get_full_error_string

endmodule particles
