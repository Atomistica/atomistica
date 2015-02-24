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

!<
!! Callables dispatch module.
!!
!! Callables are intented as plugins for the standalone code, and will be
!! called after the second Verlet step.
!! 
!>
module callables
  use libAtoms_module

  use data
  use particles
  use neighbors
  use dynamics
  use filter

  use {classname}

  implicit none

  private

  public :: callables_t
  type callables_t

     type({classtype}), allocatable :: {classname}(:)

  endtype callables_t

  ! Note: callables_t is hidden. Everything is passed as type(C_PTR) to hide
  ! the complexity of callables_t from the compiler. This speeds up compile
  ! times and avoids nasty compiler crashes. However, this invalidates Fortran
  ! interfaces since the compiler can't match a generic call to datatype.

  public :: callables_alloc, callables_free, callables_register_data
  public :: callables_init, callables_del, callables_bind_to
  public :: callables_invoke

contains

  !>
  !! Allocator
  !!
  !! Allocate memory for new callables instance
  !<
  subroutine callables_alloc(this_cptr)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(out) :: this_cptr

    ! ---

    type(callables_t), pointer :: this

    ! ---

    allocate(this)
    this_cptr = c_loc(this)

  endsubroutine callables_alloc


  !>
  !! Free memory
  !!
  !! Free memory occupied by a callables instance
  !<
  subroutine callables_free(this_cptr)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(callables_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)
    deallocate(this)

  endsubroutine callables_free


  !>
  !! Register any field with a particles object
  !!
  !! Call the register_data of all callables.
  !<
  subroutine callables_register_data(this_cptr, p)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(inout) :: p

    ! ---

    type(callables_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    call c_f_pointer(this_cptr, this)

#define REGISTER_DATA(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call register_data(this%x(i), p) ; enddo ; endif

    REGISTER_DATA({classname})

#undef REGISTER_DATA

  endsubroutine callables_register_data


  !>
  !! Constructor
  !!
  !! Call the constructors of all callables, and
  !! removes the respective lists from memory.
  !<
  subroutine callables_init(this_cptr)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(callables_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    call c_f_pointer(this_cptr, this)

#define INIT(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call init(this%x(i)) ; enddo ; endif

    INIT({classname})

#undef INIT

  endsubroutine callables_init


  !>
  !! Destructor
  !!
  !! Call the destructors of all callables, and
  !! removes the respective lists from memory.
  !<
  subroutine callables_del(this_cptr)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(callables_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    call c_f_pointer(this_cptr, this)

#define DEL(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call del(this%x(i)) ; enddo ; deallocate(this%x) ; endif

    DEL({classname})

#undef DEL

  endsubroutine callables_del


  !>
  !! Bind the callables to a certain Particles and Neighbors object
  !!
  !! Bind the callables to a certain Particles and Neighbors object. This will
  !! tell the potential to initialize its internal buffers according to the 
  !! array sizes used by the Particles and Neighbors object. All subsequent
  !! calls to energy and forces *must* be carried out with the same Particles
  !! and Neighbors object. If either Particles or Neighbors object changes,
  !! bind_to will need to be called again.
  !<
  subroutine callables_bind_to(this_cptr, p, nl, pots, coul, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    type(C_PTR),        intent(in)    :: pots
    type(C_PTR),        intent(in)    :: coul
    integer,  optional, intent(out)   :: ierror

    ! ---

    type(callables_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(this_cptr, this)

#define BIND_TO(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call bind_to(this%x(i), p, nl, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    BIND_TO({classname})

#undef BIND_TO

#define BIND_TO_WITH_POTS(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call bind_to_with_pots(this%x(i), p, nl, pots, ierror=ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    BIND_TO_WITH_POTS({classname})

#undef BIND_TO_WITH_POTS

  endsubroutine callables_bind_to


  !>
  !! Invoke the respective callables
  !!
  !! Calls the invoke method of all callables.
  !<
  subroutine callables_invoke(this_cptr, dyn, nl, pots, coul, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(dynamics_t),   intent(inout) :: dyn
    type(neighbors_t),  intent(inout) :: nl
    type(C_PTR),        intent(in)    :: pots
    type(C_PTR),        intent(in)    :: coul
    integer,  optional, intent(out)   :: ierror

    ! ---

    type(callables_t), pointer :: this

    ! ---

    integer  :: i

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(this_cptr, this)


#define INVOKE(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call invoke(this%x(i), dyn, nl, ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    INVOKE({classname})

#undef INVOKE

#define INVOKE_WITH_POTS(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call invoke_with_pots(this%x(i), dyn, nl, pots, ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    INVOKE_WITH_POTS({classname})

#undef INVOKE_WITH_POTS

#define INVOKE_WITH_POTS_AND_COUL(x)  if (allocated(this%x)) then ; do i = lbound(this%x, 1), ubound(this%x, 1) ; call invoke_with_pots_and_coul(this%x(i), dyn, nl, pots, coul, ierror) ; PASS_ERROR(ierror) ; enddo ; endif

    INVOKE_WITH_POTS_AND_COUL({classname})

#undef INVOKE_WITH_POTS_AND_COUL

  endsubroutine callables_invoke

endmodule callables

