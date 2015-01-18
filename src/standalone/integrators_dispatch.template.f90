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
!>
!! Integrator dispatch module
!<

#include "macros.inc"

module integrators
  use libAtoms_module

  use logging

  use particles
  use dynamics

  use {classname}  

  implicit none
  
  private

  !
  ! Dispatch type
  !

  public :: integrators_t
  type integrators_t

     type({classtype}), allocatable :: {classname}

  endtype integrators_t

  !
  ! Integrators
  !

  public :: integrators_init, integrators_del, integrators_step1
  public :: integrators_step2

contains

  !>
  !! Constructor
  !<
  subroutine integrators_init(this_cptr, p)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(inout) :: p

    ! ---

    type(integrators_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    if (.not. any( (/ &
         allocated(this%{classname}), &
         .false. /) ) &
         ) then

       !
       ! Use the velocity verlet integrator by default
       !

       call prlog("- integrators_init -")
       call prlog("     Using default (Verlet) integrator.")
       call prlog

       allocate(this%verlet)
    endif

#define INIT(x)  if (allocated(this%x)) then ; call init(this%x, p) ; endif

    INIT({classname})

#undef INIT
    
  endsubroutine integrators_init


  !>
  !! Destructor
  !<
  subroutine integrators_del(this_cptr)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR), intent(in) :: this_cptr

    ! ---

    type(integrators_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

#define DEL(x)  if (allocated(this%x)) then ; call del(this%x) ; endif

    DEL({classname})

#undef DEL

  endsubroutine integrators_del


  !>
  !! Position update and velocity estimation
  !<
  subroutine integrators_update(this_cptr, p)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(particles_t), intent(inout) :: p

    ! ---

    type(integrators_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

#define UDPATE(x)  if(allocated(this%x)) then ; call update(this%x, p) ; endif

    UPDATE({classname})

#undef UPDATE

  endsubroutine integrators_update


  !>
  !! Step1, called before force calculation
  !!
  !! Step1, called before force calculation. This is the position update and
  !! the velocity estimation in the usual velocity-Verlet algorithm.
  !<
  subroutine integrators_step1(this_cptr, pots_cptr, dyn, max_dt, max_dr, &
       max_dr_sq_out, ierror)
    use, intrinsic :: iso_c_binding

    use potentials

    implicit none

    type(C_PTR),        intent(in)    :: this_cptr
    type(C_PTR),        intent(in)    :: pots_cptr
    type(dynamics_t),   intent(inout) :: dyn
    real(DP), optional, intent(in)    :: max_dt
    real(DP), optional, intent(in)    :: max_dr
    real(DP), optional, intent(inout) :: max_dr_sq_out
    integer,  optional, intent(inout) :: ierror

    ! ---

    type(integrators_t), pointer :: this
    type(potentials_t), pointer :: pots

    ! ---

    call c_f_pointer(this_cptr, this)
    call c_f_pointer(pots_cptr, pots)

#define STEP1(x)  if (allocated(this%x)) then ; call step1(this%x, dyn%p, dyn%v, dyn%f, dyn%dt, max_dt, max_dr, max_dr_sq_out) ; endif

    STEP1({classname})

#undef STEP1

#define STEP1_WITH_DYN(x)  if (allocated(this%x)) then ; call step1_with_dyn(this%x, dyn, max_dt, max_dr, max_dr_sq_out) ; endif

    STEP1_WITH_DYN({classname})

#undef STEP1_WITH_DYN

#define STEP1_WITH_BAROSTAT(x)  if (allocated(this%x)) then ; call step1_with_barostat(this%x, pots%sliding_p, dyn%p, dyn%v, dyn%f, dyn%dt, max_dt, max_dr, max_dr_sq_out) ; endif

    STEP1_WITH_BAROSTAT({classname})

#undef STEP1_WITH_BAROSTAT

    if (allocated(this%sliding_t)) then
       if (.not. allocated(pots%sliding_p)) then
          RAISE_ERROR("SlidingT can only be used in conjunction with SlidingP.", ierror)
       else
          if (size(pots%sliding_p) > 1) then
             RAISE_ERROR("There can only be a single SlidingP object.", ierror)
          endif
       endif
    endif

  endsubroutine integrators_step1


  !>
  !! Step 2, called after force calculation (Velocity correction
  !!
  !! Step 2, called after force calculation. This is the velocity correction
  !! in the usual velocity-Verlet algorithm.
  !<
  subroutine integrators_step2(this_cptr, dyn, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),       intent(in)    :: this_cptr
    type(dynamics_t),  intent(inout) :: dyn
    integer, optional, intent(inout) :: ierror

    ! ---

    type(integrators_t), pointer :: this

    ! ---

    call c_f_pointer(this_cptr, this)

#define STEP2(x)  if (allocated(this%x)) then ; call step2(this%x, dyn%p, dyn%v, dyn%f, dyn%dt) ; endif

    STEP2({classname})

#undef STEP2_WITH_DYN

#define STEP2_WITH_DYN(x)  if (allocated(this%x)) then ; call step2_with_dyn(this%x, dyn) ; endif

    STEP2_WITH_DYN({classname})

#undef STEP2_WITH_DYN
    
#define STEP2_WITH_WPOT(x)  if (allocated(this%x)) then ; call step2(this%x, dyn%p, dyn%v, dyn%f, dyn%wpot, dyn%dt) ; endif

    STEP2_WITH_WPOT({classname})

#undef STEP2_WITH_WPOT

  endsubroutine integrators_step2

endmodule integrators
