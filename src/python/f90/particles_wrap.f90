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
module particles_wrap
  use libAtoms_module

  use data
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

  subroutine f_particles_allocate(this_cptr, nat, ierror) &
       bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),     value          :: this_cptr
    integer(c_int),  value          :: nat
    integer(c_int),  intent(inout)  :: ierror

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call allocate(this, nat, ierror=ierror)

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


  subroutine f_particles_set_cell(this_cptr, cell, pbc, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),     value          :: this_cptr
    real(c_double)                  :: cell(3, 3)
    logical(c_bool)                 :: pbc(3)
    integer(c_int),  intent(inout)  :: ierror

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call set_cell(this, cell, logical(pbc), ierror=ierror)

  endsubroutine f_particles_set_cell


  subroutine f_particles_set_lees_edwards(this_cptr, dx, dv, ierror) bind(C)
    use, intrinsic :: iso_c_binding

    implicit none

    type(c_ptr),     value          :: this_cptr
    real(c_double)                  :: dx(3)
    real(c_double)                  :: dv(3)
    integer(c_int),  intent(inout)  :: ierror

    ! ---

    type(particles_t), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    call set_lees_edwards(this, dx, dv, ierror=ierror)

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

endmodule particles_wrap
