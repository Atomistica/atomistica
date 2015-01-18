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
!   classtype:output_xyz_t classname:OutputXYZ interface:callables
! @endmeta

!>
!! The XYZ output module
!!
!! The XYZ output module.
!! Note: This is a thin layer that has the "callables" interface and hooks into
!! the standalone code. The output modules from src/io are used.
!<

#include "macros.inc"

module output_xyz
  use libAtoms_module

  use particles
  use dynamics
  use neighbors

  use xyz

  implicit none

  private

  public :: output_xyz_t
  type output_xyz_t
     
     !
     ! Interval in which to output another trajectory step
     !
		 
     real(DP) :: freq = 1.0_DP

     !
     ! Time
     !

     real(DP) :: ti
		
     !
     ! Output file
     !
		
     integer  :: un
	
  endtype output_xyz_t


  public :: init
  interface init
     module procedure output_xyz_init
  endinterface

  public :: del
  interface del
     module procedure output_xyz_del
  endinterface

! NOT REQUIRED
!  interface register_data
!     module procedure output_xyz_register_data
!  endinterface

  public :: invoke
  interface invoke
     module procedure output_xyz_invoke
  endinterface

  public :: register
  interface register
    module procedure output_xyz_register
  endinterface

contains 

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine output_xyz_init(this)
    implicit none

    type(output_xyz_t), intent(inout)  :: this

    ! ---
    
    this%un  = fopen("traj.xyz", F_WRITE)
   
    this%ti  = 0.0_DP
    
  endsubroutine output_xyz_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine output_xyz_del(this)
    implicit none

    type(output_xyz_t), intent(inout)  :: this

    ! ---

    call fclose(this%un)

  endsubroutine output_xyz_del


  !>
  !! Output a snapshot
  !!
  !! Output a snapshot
  !<
  subroutine output_xyz_invoke(this, dyn, nl, ierror)
    implicit none

    type(output_xyz_t), intent(inout)  :: this
    type(dynamics_t), intent(in)       :: dyn
    type(neighbors_t), intent(in)      :: nl
    integer, intent(inout), optional   :: ierror

    ! ---

    this%ti = this%ti + dyn%dt
    
    if (this%ti >= this%freq) then
       call write_xyz(this%un, dyn%p, length_to_A)
       this%ti = 0.0_DP
    endif

  endsubroutine output_xyz_invoke


  subroutine output_xyz_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(output_xyz_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)           :: cfg
    type(c_ptr), intent(out)          :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("OutputXYZ"), &
         CSTR("XYZ output module."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Output interval."))

  endsubroutine output_xyz_register

endmodule output_xyz
