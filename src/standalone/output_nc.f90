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
!   classtype:output_nc_t classname:OutputNC interface:callables
! @endmeta

!>
!! The NetCDF output module (AMBER-style)
!!
!! The NetCDF output module (AMBER-style).
!! Note: This is a thin layer that has the "callables" interface and hooks into
!! the standalone code. The output modules from src/io are used.
!<

#include "macros.inc"

module output_nc
  use supplib

  use particles
  use dynamics
  use neighbors

  use nc

  implicit none

  private

  character(MAX_NAME_STR), parameter  :: TI_STR = "output_nc.time"
  character(MAX_NAME_STR), parameter  :: EPOT_PER_AT_STR = "potential_energy"

  public :: output_nc_t
  type output_nc_t

     type(dynamics_t), pointer  :: dyn

     !
     ! Interval in which to output another trajectory step
     !

     real(DP)           :: freq  = -1.0_DP
     logical(BOOL)      :: epot_per_at = .false.

     !
     ! Time
     !

     real(DP), pointer  :: ti

     !
     ! NetCDF stuff
     !

     type(nc_t)         :: nc

  endtype output_nc_t


  public :: init
  interface init
     module procedure output_nc_init
  endinterface

  public :: del
  interface del
     module procedure output_nc_del
  endinterface

  public :: register_data
  interface register_data
     module procedure output_nc_register_data
  endinterface

  public :: bind_to_with_pots
  interface bind_to_with_pots
     module procedure output_nc_bind_to
  endinterface

  public :: invoke
  interface invoke
     module procedure output_nc_invoke
  endinterface

  public :: register
  interface register
    module procedure output_nc_register
  endinterface

!--- Internal

  interface set_dynamics
     module procedure output_nc_set_dynamics
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine output_nc_init(this)
    implicit none

    type(output_nc_t), intent(inout)  :: this

    ! ---

    nullify(this%dyn)
    nullify(this%ti)

  endsubroutine output_nc_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine output_nc_del(this)
    implicit none

    type(output_nc_t), intent(inout)  :: this

    ! ---

    if (associated(this%dyn)) then
       call write_frame(this%nc, this%dyn%ti, this%dyn%p)
       call close(this%nc)
    endif

  endsubroutine output_nc_del


  !>
  !! Register any additonal fields with the particles object
  !!
  !! Register any additonal fields with the particles object
  !<
  subroutine output_nc_register_data(this, p)
    implicit none

    type(output_nc_t), intent(inout)  :: this
    type(particles_t), intent(inout)  :: p

    ! ---

    call add_real_attr(p%data, TI_STR)
    if (this%epot_per_at) then
      call add_real(p%data, EPOT_PER_AT_STR, F_TO_TRAJ)
    endif

  endsubroutine output_nc_register_data


  !>
  !! Notify output object of particle and neighbor list objects
  !!
  !! Notify output object of particle and neighbor list objects
  !<
  subroutine output_nc_bind_to(this, p, nl, pots_cptr, ierror)
    use, intrinsic :: iso_c_binding
    use potentials

    implicit none

    type(output_nc_t),  intent(inout) :: this
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    type(C_PTR),        intent(in)    :: pots_cptr
    integer,  optional, intent(out)   :: ierror

    ! ---

    type(potentials_t), pointer :: pots

    ! ---

    INIT_ERROR(ierror)
    call c_f_pointer(pots_cptr, pots)

    if (this%epot_per_at) then
      if (associated(pots%epot_per_at)) then
        RAISE_ERROR("Another module allocated per-atom potential energy array.", ierror)
      else
        call ptr_by_name(p%data, EPOT_PER_AT_STR, pots%epot_per_at)
      endif
    endif

  endsubroutine output_nc_bind_to


  !>
  !! Open output file, etc.
  !!
  !! Open output file, etc.
  !<
  subroutine output_nc_set_dynamics(this, dyn)
    implicit none

    type(output_nc_t), intent(inout)  :: this
    type(dynamics_t),  target         :: dyn

    ! ---

    this%dyn  => dyn

    call attr_by_name(dyn%p%data, TI_STR, this%ti)

    call create(this%nc, dyn%p, "traj.nc")

#ifndef _MP
    call write_prmtop(dyn%p, "traj.prmtop")
#endif
    call write_frame(this%nc, dyn%ti, dyn%p)

  endsubroutine output_nc_set_dynamics


  !>
  !! Output a snapshot
  !!
  !! Output a snapshot
  !<
  subroutine output_nc_invoke(this, dyn, nl, ierror)
    implicit none

    type(output_nc_t), intent(inout)  :: this
    type(dynamics_t),  target         :: dyn
    type(neighbors_t), target         :: nl
    integer, optional, intent(out)    :: ierror

    ! ---

    INIT_ERROR(ierror)

    if (.not. associated(this%dyn, dyn)) then
       call set_dynamics(this, dyn)
    endif

    this%ti = this%ti + dyn%dt

    if (this%ti >= this%freq) then
       call write_frame(this%nc, dyn%ti, dyn%p)
       this%ti = 0.0_DP
    endif

  endsubroutine output_nc_invoke


  subroutine output_nc_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(output_nc_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)          :: cfg
    type(c_ptr), intent(out)         :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("OutputNC"), &
         CSTR("NetCDF output module."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Output interval."))
    call ptrdict_register_boolean_property(m, c_loc(this%epot_per_at), &
         CSTR("epot"), CSTR("Output per-atom potential energies."))

  endsubroutine output_nc_register

endmodule output_nc
