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
!   classtype:output_cell_t classname:OutputCell interface:callables
! @endmeta

!>
!! Continuously output the cell shape
!!
!! Continuously output the cell shape
!<

#include "macros.inc"

module output_cell
  use libAtoms_module

  use io
  use logging

  use particles
  use dynamics
  use neighbors

  implicit none

  private

  public :: output_cell_t
  type output_cell_t

     real(DP)  :: freq  = -1.0_DP

     integer   :: un

     !
     ! Averaging
     !

     real(DP)  :: t

     real(DP)  :: s(3)
     real(DP)  :: dx(3)

  endtype output_cell_t


  public :: init
  interface init
     module procedure output_cell_init
  endinterface

  public :: del
  interface del
     module procedure output_cell_del
  endinterface

  public :: invoke
  interface invoke
     module procedure output_cell_invoke
  endinterface

  public :: register
  interface register
    module procedure output_cell_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine output_cell_init(this)
    implicit none

    type(output_cell_t), intent(inout)  :: this

    ! ---

#ifdef _MP
    if (mpi_id() == 0) then
#endif

    write (ilog, '(A)')            "- output_cell_init -"
    write (ilog, '(5X,A,F20.10)')  "freq = ", this%freq

    this%un = fopen("cell.out", F_WRITE)

    write (this%un, '(A1,1X,A8,7A20)')  "#", "1:it", "2:time", &
         "3:sx", "4:sy", "5:sz", "6:dx", "7:dy", "8:dz"

    this%t   = 0.0_DP
    this%s   = 0.0_DP
    this%dx  = 0.0_DP

    write (ilog, *)

#ifdef _MP
    endif
#endif

  endsubroutine output_cell_init


  !>
  !! Destructor
  !!
  !! Delete a output_cell object
  !<
  subroutine output_cell_del(this)
    implicit none

    type(output_cell_t), intent(inout)  :: this

    ! ---

#ifdef _MP
    if (mpi_id() == 0) then
#endif

    call fclose(this%un)

#ifdef _MP
    endif
#endif

  endsubroutine output_cell_del


  !>
  !! Output the cells dimensions
  !!
  !! Output the cells dimensions
  !<
  subroutine output_cell_invoke(this, dyn, nl, ierror)
    implicit none

    type(output_cell_t), intent(inout)  :: this
    type(dynamics_t), intent(in)        :: dyn
    type(neighbors_t), intent(in)       :: nl
    integer, intent(inout), optional    :: ierror

    ! ---

    this%t = this%t + dyn%dt

    this%s   = this%s   + (/ dyn%p%Abox(1, 1), dyn%p%Abox(2, 2), dyn%p%Abox(3, 3) /) * dyn%dt
    this%dx  = this%dx  + dyn%p%shear_dx * dyn%dt

    if (this%freq < 0 .or. this%t >= this%freq) then

       this%s   = this%s / this%t
       this%dx  = this%dx / this%t

       write (this%un, '(I9,X,7ES20.10)')  dyn%it, dyn%ti, this%s, this%dx

       this%t   = 0.0_DP
       this%s   = 0.0_DP
       this%dx  = 0.0_DP

    endif

  endsubroutine output_cell_invoke


  subroutine output_cell_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(output_cell_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)             :: cfg
    type(c_ptr), intent(out)            :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("OutputCell"), &
         CSTR("Output the cell size and Lees-Edwards displacement 'cell.out'."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Output frequency (-1 means output every time step)."))

  endsubroutine output_cell_register

endmodule output_cell
