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
!! Binning and neighbor list module
!<
module neighbors
  use, intrinsic :: iso_c_binding

  use supplib

  use particles

  implicit none

  private

  public :: neighptr_t
  integer, parameter :: NEIGHPTR_T = C_INTPTR_T

  public :: neighbors_t
  type neighbors_t

     !>
     !! Seed for the neighbor list for the first set of neighbors
     !<
     integer(NEIGHPTR_T), pointer :: seed(:)

     !>
     !! End type neighbors_tthe neighbor list for the first set of neighbors
     !<
     integer(NEIGHPTR_T), pointer :: last(:)

     !>
     !! Size of the neighbor list
     !<
     integer                      :: neighbors_size

     !>
     !! Neighbor list for the second set of neighbors
     !<
     integer(C_INT),      pointer :: neighbors(:)

     !>
     !! Neighbor list cutoffs
     !<
     real(DP),        allocatable :: cutoff(:, :)

  endtype neighbors_t

  public :: request_interaction_range
  interface request_interaction_range
     module procedure neighbors_request_interaction_range
  endinterface

  public :: update
  interface update
     module procedure neighbors_update
  endinterface

contains
  
  !>
  !! Create a new instance
  !<
  subroutine neighbors_new(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(inout) :: this_cptr

    ! ---

    type(neighbors_t), pointer :: this_fptr

    ! ---

    allocate(this_fptr)
    this_cptr = c_loc(this_fptr)

  endsubroutine neighbors_new


  !>
  !! Destroy the instance
  !<
  subroutine neighbors_free(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr

    ! ---

    type(neighbors_t), pointer :: this_fptr

    ! ---

    call c_f_pointer(this_cptr, this_fptr)
    deallocate(this_fptr)

  endsubroutine neighbors_free


  !>
  !! Constructor
  !<
  subroutine neighbors_init(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr

    ! ---

    type(neighbors_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    allocate(this%cutoff(MAX_Z, MAX_Z))
    this%cutoff = 0.0_DP

  endsubroutine neighbors_init


  !>
  !! Destructor
  !<
  subroutine neighbors_del(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr

    ! ---

    type(neighbors_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    deallocate(this%cutoff)

  endsubroutine neighbors_del


  !>
  !! Assign pointers to data
  !>
  subroutine neighbors_set_pointers(this_cptr, nat, seed, last, &
       neighbors_size, neighbors) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),    value :: this_cptr
    integer(C_INT), value :: nat
    type(C_PTR),    value :: seed, last
    integer(C_INT), value :: neighbors_size
    type(C_PTR),    value :: neighbors

    ! ---

    type(neighbors_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    call c_f_pointer(seed, this%seed, [nat])
    call c_f_pointer(last, this%last, [nat])
    this%neighbors_size = neighbors_size
    call c_f_pointer(neighbors, this%neighbors, [neighbors_size])

  endsubroutine neighbors_set_pointers


  !>
  !! Dummy subroutine. LAMMPS should ensure the neighbor list has been build
  !<
  subroutine neighbors_update(this, p, error)
    implicit none

    type(neighbors_t), intent(inout) :: this
    type(particles_t), intent(in)    :: p
    integer, optional, intent(inout) :: error

    ! ---

  endsubroutine neighbors_update


  !>
  !! Request and interaction range
  !<
  subroutine neighbors_request_interaction_range(this, cutoff, el1, el2)
    implicit none

    type(neighbors_t), intent(inout) :: this
    real(DP),          intent(in)    :: cutoff
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
       this%cutoff(el1, rel2) = max(this%cutoff(el1, rel2), cutoff)
       if (el1 /= rel2)  then
          this%cutoff(rel2, el1) = max(this%cutoff(rel2, el1), cutoff)
       endif
    else
       this%cutoff = max(this%cutoff, cutoff)
    endif

  endsubroutine neighbors_request_interaction_range


  !>
  !! Return cutoffs for element combination el1, el2
  !<
  subroutine neighbors_get_cutoff(this_cptr, el1, el2, cutoff) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),    value       :: this_cptr
    integer(C_INT), value       :: el1, el2
    real(DP),       intent(out) :: cutoff

    ! ---

    type(neighbors_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    cutoff = this%cutoff(el1, el2)

  endsubroutine neighbors_get_cutoff


  !>
  !! Dump cutoffs to the log file
  !<
  subroutine neighbors_dump_cutoffs(this_cptr, p_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value :: this_cptr, p_cptr

    ! ---

    type(neighbors_t), pointer :: this
    type(particles_t), pointer :: p

    integer :: i, j

    ! ---

    call c_f_pointer(this_cptr, this)
    call c_f_pointer(p_cptr, p)

    call prlog("- neighbors_dump_cutoffs -")
    call prlog("     Communication border is " // p%border)

    do i = 1, p%nel
       do j = i, p%nel
          call prlog("     " // i // "-" // j // " (Z = " // p%el2Z(i) // &
               "-" // p%el2Z(j) // "), cutoff = " // this%cutoff(i, j))
       enddo
    enddo

    call prlog

  endsubroutine neighbors_dump_cutoffs

endmodule neighbors
