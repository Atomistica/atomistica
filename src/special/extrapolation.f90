!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
!! and others. See the AUTHORS file in the top-level Atomistica directory.
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
!! Extrapolation
!!
!! Extrapolation of properties from past time steps to the current time step.
!! The extrapolation scheme is the one of Alfe, Comp. Phys. Comp. 118, 31
!! (1999), originally desgined for the extrapolation charge densities in
!! density functional theory calculations.
!!
!<

#include "macros.inc"

module extrapolation
  use supplib

  use particles

  implicit none

  private

  public :: extrapolation_t
  type extrapolation_t

     integer               :: extrapolation_memory = 3  !< Number of past steps to keep

     integer               :: history_counter = 0
     real(DP), allocatable :: r(:, :, :)
     real(DP), allocatable :: q(:, :)

  endtype extrapolation_t


  public :: init
  interface init
     module procedure extrapolation_init
  endinterface

  public :: del
  interface del
     module procedure extrapolation_del
  endinterface

  public :: extrapolate
  interface extrapolate
     module procedure extrapolation_extrapolate
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine extrapolation_init(this, p, extrapolation_memory)
    implicit none

    type(extrapolation_t), intent(inout) :: this
    type(particles_t),     intent(in)    :: p
    integer,               intent(in)    :: extrapolation_memory

    ! ---

    this%extrapolation_memory = extrapolation_memory

    this%history_counter = 0
    if (allocated(this%r)) then
       deallocate(this%r)
    endif
    if (allocated(this%q)) then
       deallocate(this%q)
    endif
    if (this%extrapolation_memory >= 2) then
       allocate(this%r(3, p%maxnatloc, this%extrapolation_memory))
       allocate(this%q(p%maxnatloc, this%extrapolation_memory))
    endif
    
  endsubroutine extrapolation_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine extrapolation_del(this)
    implicit none

    type(extrapolation_t), intent(inout)  :: this

    ! ---

    if (allocated(this%r))  deallocate(this%r)
    if (allocated(this%q))  deallocate(this%q)
    this%extrapolation_memory = 0
    
  endsubroutine extrapolation_del


  !>
  !! Extrapolate a scalar quantity from past time steps
  !<
  subroutine extrapolation_extrapolate(this, p, q, error)
    implicit none

    type(extrapolation_t), intent(inout) :: this      !< extrapolation object
    type(particles_t),     intent(in)    :: p         !< Particles
    real(DP),              intent(inout) :: q(p%nat)  !< Scalar quantity
    integer,     optional, intent(out)   :: error     !< Error status

    ! ---

    integer :: i, k, l

    real(DP) :: a(this%extrapolation_memory-1, this%extrapolation_memory-1)
    real(DP) :: b(this%extrapolation_memory-1)
    real(DP) :: alpha(this%extrapolation_memory-1)
    real(DP) :: q0(p%nat)

    real(DP) :: drk(3, p%natloc), drl(3, p%natloc)

    ! ---

    INIT_ERROR(error)

    if (this%extrapolation_memory < 2) then
       return
    endif

    q0 = q

    if (this%history_counter >= this%extrapolation_memory) then
       do k = 1, this%extrapolation_memory-1
          drk = this%r(1:3, 1:p%natloc, modulo(this%history_counter-k, this%extrapolation_memory)+1) - &
                this%r(1:3, 1:p%natloc, modulo(this%history_counter-k-1, this%extrapolation_memory)+1)
          do l = 1, this%extrapolation_memory-1
             drl = this%r(1:3, 1:p%natloc, modulo(this%history_counter-l, this%extrapolation_memory)+1) - &
                   this%r(1:3, 1:p%natloc, modulo(this%history_counter-l-1, this%extrapolation_memory)+1)
             a(k, l) = dot_product(reshape(drk, [3*p%natloc]), reshape(drl, [3*p%natloc]))
          enddo
          b(k) = dot_product( &
              reshape(PCN3(p, 1:p%natloc) - &
                      this%r(1:3, 1:p%natloc, modulo(this%history_counter-1, this%extrapolation_memory)+1), &
                      [3*p%natloc]), &
              reshape(this%r(1:3, 1:p%natloc, modulo(this%history_counter-k, this%extrapolation_memory)+1) - &
                      this%r(1:3, 1:p%natloc, modulo(this%history_counter-k-1, this%extrapolation_memory)+1), &
                      [3*p%natloc]) &
              )
       enddo

       alpha = matmul(inverse(a, error=error), b)
       if (error == ERROR_NONE) then
          !                q(t) - q(t-dt)
          q = q0 + alpha(1)*(q0 - this%q(1:p%nat, modulo(this%history_counter-1, this%extrapolation_memory)+1))
          do i = 2, this%extrapolation_memory-1
             q = q + alpha(i)*(this%q(1:p%nat, modulo(this%history_counter-i+1, this%extrapolation_memory)+1) - &
                               this%q(1:p%nat, modulo(this%history_counter-i, this%extrapolation_memory)+1))
          enddo
       else
          call prlog("Warning: Unable to extrapolate quantity. Resetting history.")
          call clear_error(error)
          this%history_counter = 0
       endif
    endif

    this%history_counter = this%history_counter+1
    i = modulo(this%history_counter-1, this%extrapolation_memory)+1

    ! This is current r(t+dt)
    this%r(1:3, 1:p%nat, i) = PCN3(p, 1:p%nat)

    ! This is last q(t)
    this%q(1:p%nat, i) = q0(1:p%nat)
    
  endsubroutine extrapolation_extrapolate

endmodule extrapolation
