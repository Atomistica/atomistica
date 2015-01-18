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
!   classtype:output_cfg_t classname:OutputCFG interface:callables
! @endmeta

!>
!! The CFG output module (AtomEye)
!!
!! The CFG output module (AtomEye).
!! Note: This is a thin layer that has the "callables" interface and hooks into
!! the standalone code. The output modules from src/io are used.
!<

#include "macros.inc"

module output_cfg
  use libAtoms_module

  use timer

  use particles
  use dynamics
  use neighbors

  use cfg

  implicit none

  private

  public :: output_cfg_t
  type output_cfg_t
	
     !
     ! Configuration variables
     !
		 
     real(DP)  :: freq  = -1.0_DP

     real(DP)  :: ti
     integer   :: n
		
  endtype output_cfg_t


  public :: init
  interface init
     module procedure output_cfg_init
  endinterface

! NOT REQUIRED
!  interface del
!     module procedure output_cfg_del
!  endinterface

! NOT REQUIRED
!  interface register_data
!     module procedure output_cfg_register_data
!  endinterface

  public :: invoke
  interface invoke
     module procedure output_cfg_invoke
  endinterface

  public :: register
  interface register
    module procedure output_cfg_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine output_cfg_init(this)
    implicit none

    type(output_cfg_t), intent(inout)  :: this

    ! ---

    this%ti  = 0.0_DP

    this%n   = 0
    
  endsubroutine output_cfg_init


  !>
  !! Output a snapshot
  !!
  !! Output a snapshot
  !<
  subroutine output_cfg_invoke(this, dyn, nl, ierror)
    implicit none

    type(output_cfg_t), intent(inout)  :: this
    type(dynamics_t), intent(in)       :: dyn
    type(neighbors_t), intent(in)      :: nl
    integer, intent(inout), optional   :: ierror

    ! ---

    character(9)  :: fn

    ! ---

    call timer_start("output_cfg_invoke")

    this%ti  = this%ti + dyn%dt

    if (this%ti > this%freq) then

       this%n  = this%n + 1

       write (fn, '(I5.5,A4)')  this%n, ".cfg"

       call write_cfg(fn, dyn%p)

       this%ti  = 0.0_DP

    endif

    call timer_stop("output_cfg_invoke")

  endsubroutine output_cfg_invoke


  subroutine output_cfg_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(output_cfg_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)           :: cfg
    type(c_ptr), intent(out)          :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("OutputCFG"), &
         CSTR("CFG output module (AtomEye extended format)."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq" ), &
         CSTR("Output interval."))

  endsubroutine output_cfg_register

endmodule output_cfg
