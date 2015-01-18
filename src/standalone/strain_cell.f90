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
!   classtype:strain_cell_t classname:StrainCell interface:callables
! @endmeta

!>
!! Continuously deform the unit cell
!!
!! Continuously deform the unit cell using a constant strain rate
!<

#include "macros.inc"

module strain_cell
  use libAtoms_module

  use particles
  use neighbors
  use dynamics
  
  implicit none

  private

  public :: strain_cell_t
  type strain_cell_t

     !
     ! Hydrostatic pressure components
     !

     real(DP)  :: gamma(3)  = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
     real(DP)  :: dcell(3)  = (/ 0.0_DP, 0.0_DP, 0.0_DP /)

  endtype strain_cell_t


  public :: init
  interface init
     module procedure strain_cell_init
  endinterface

  public :: del
  interface del
     module procedure strain_cell_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure strain_cell_bind_to
  endinterface bind_to

  public :: invoke
  interface invoke
     module procedure strain_cell_invoke
  endinterface

  public :: register
  interface register
    module procedure strain_cell_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a DeformCell object
  !<
  subroutine strain_cell_init(this, gamma, dcell)
    implicit none

    type(strain_cell_t), intent(inout)  :: this
    real(DP),  optional, intent(in)     :: gamma(3)
    real(DP),  optional, intent(in)     :: dcell(3)

    ! ---

    call prlog("- strain_cell_init -")
    call prlog("     $Id$")

    ASSIGN_PROPERTY(gamma)
    ASSIGN_PROPERTY(dcell)

    call prlog("     Straining unit cell with parameters")
    call prlog("     gamma  = " // this%gamma)
    call prlog("     dcell  = " // this%dcell)

  endsubroutine strain_cell_init


  !>
  !! Destructor
  !!
  !! Destroy a DeformCell object
  !<
  subroutine strain_cell_del(this)
    implicit none

    type(strain_cell_t), intent(inout)   :: this

    ! ---

  end subroutine strain_cell_del


  !>
  !! Notify the ContactArea object of the Particles and Neighbors objects
  !<
  subroutine strain_cell_bind_to(this, p, nl, ierror)
    implicit none

    type(strain_cell_t), intent(inout)  :: this
    type(particles_t),   intent(inout)  :: p
    type(neighbors_t),   intent(inout)  :: nl
    integer,   optional, intent(out)    :: ierror

    ! ---

    INIT_ERROR(ierror)

    call require_orthorhombic_cell(p, error=ierror)
    PASS_ERROR(ierror)

  endsubroutine strain_cell_bind_to


  !>
  !! Change the box size
  !!
  !! Change the box size
  !<
  subroutine strain_cell_invoke(this, dyn, nl, error)
    implicit none

    type(strain_cell_t), intent(inout)  :: this
    type(dynamics_t),    intent(inout)  :: dyn
    type(neighbors_t),   intent(in)     :: nl
    integer,   optional, intent(inout)  :: error

    ! ---

    real(DP)  :: dc(3)

    ! ---

    dc = ( matmul(dyn%p%Abox, this%gamma) + this%dcell )*dyn%dt

    call set_cell( &
         dyn%p, &
         (/ dyn%p%Abox(1, 1) + dc(1), &
            dyn%p%Abox(2, 2) + dc(2), &
            dyn%p%Abox(3, 3) + dc(3) /), &
         error=error)

  endsubroutine strain_cell_invoke


  subroutine strain_cell_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(strain_cell_t), target, intent(inout)  :: this
    type(c_ptr),                 intent(in)     :: cfg
    type(c_ptr),                 intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("StrainCell"), &
         CSTR("Continuously deform the unit cell."))

    call ptrdict_register_point_property(m, c_loc(this%gamma(1)), &
         CSTR("gamma"), &
         CSTR("Strain rate (i.e. dcell = gamma*cell)."))
    call ptrdict_register_point_property(m, c_loc(this%dcell(1)), &
         CSTR("dcell"), &
         CSTR("Absolute change in cell size (per time unit)."))

  endsubroutine strain_cell_register

endmodule strain_cell
