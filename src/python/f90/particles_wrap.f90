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

module particles_wrap
  use supplib

  use particles

  implicit none

contains

  subroutine f_particles_new(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), intent(out)  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---
   
    allocate(this)
    this_cptr = c_loc(this)

  endsubroutine f_particles_new


  subroutine f_particles_free(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    deallocate(this)

  endsubroutine f_particles_free


  subroutine f_particles_init(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call init(this)

  endsubroutine f_particles_init

  subroutine f_particles_allocate(this_cptr, nat, error) &
       bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),     value          :: this_cptr
    integer(c_int),  value          :: nat
    integer(c_int),  intent(inout)  :: error

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call allocate(this, nat, error=error)

  endsubroutine f_particles_allocate


  subroutine f_particles_del(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call del(this)

  endsubroutine f_particles_del


  subroutine f_particles_update_elements(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call update_elements(this)

  endsubroutine f_particles_update_elements


  subroutine f_particles_set_cell(this_cptr, cell, pbc, error) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),     value          :: this_cptr
    real(c_double)                  :: cell(3, 3)
    logical(c_bool)                 :: pbc(3)
    integer(c_int),  intent(inout)  :: error

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call set_cell(this, cell, logical(pbc), error=error)

  endsubroutine f_particles_set_cell


  subroutine f_particles_set_lees_edwards(this_cptr, dx, dv, error) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),     value          :: this_cptr
    real(c_double)                  :: dx(3)
    real(c_double)                  :: dv(3)
    integer(c_int),  intent(inout)  :: error

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call set_lees_edwards(this, dx, dv, error=error)

  endsubroutine f_particles_set_lees_edwards


  subroutine f_particles_inbox(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call inbox(this)

  endsubroutine f_particles_inbox


  subroutine f_particles_I_changed_positions(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value  :: this_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call I_changed_positions(this)

  endsubroutine f_particles_I_changed_positions


  !>
  !! Get data structure elements of a particle_t.
  !<
  subroutine f_particles_get_data(this_cptr, data_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value        :: this_cptr
    type(c_ptr), intent(out)  :: data_cptr

    ! ---

    type(particles_t), pointer  :: this

    ! ---
   
    call c_f_pointer(this_cptr, this)
    data_cptr = c_loc(this%data)

  endsubroutine f_particles_get_data
  
  
  !>
  !! Set the tag that is stored in the particles_t
  !<
  subroutine f_particles_set_tag(this_cptr, tag) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value  :: this_cptr
    type(C_PTR), value  :: tag

    ! ---

    type(particles_t), pointer  :: this

    ! ---
   
    call c_f_pointer(this_cptr, this)
    this%tag = tag

  endsubroutine f_particles_set_tag

 
  !>
  !! Get the tag that is stored in the particles_t
  !<
  subroutine f_particles_get_tag(this_cptr, tag) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), value        :: this_cptr
    type(C_PTR), intent(out)  :: tag

    ! ---

    type(particles_t), pointer  :: this

    ! ---
   
    call c_f_pointer(this_cptr, this)
    tag = this%tag

  endsubroutine f_particles_get_tag


  !>
  !! Return number of elements of this particles object
  !<
  function f_particles_get_nel(this_cptr) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr), value :: this_cptr

    integer(c_int)     :: f_particles_get_nel

    ! ---

    type(particles_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    f_particles_get_nel  = this%nel

  endfunction f_particles_get_nel

endmodule particles_wrap
