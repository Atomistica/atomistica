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
!   shared:directory
! @endmeta

!>
!! Dispatch module for the eigenvalue solvers
!!
!! Dispatch module for the eigensystem solver. Note that the solver computes
!! the density matrix, not the eigenvectors!
!<

#include "macros.inc"

module dense_solver
  use, intrinsic :: iso_c_binding

  use supplib

  use particles

  use materials

  use dense_hamiltonian_type

  use dense_solver_cp
  use dense_solver_lapack

  implicit none

  private

  public :: dense_solver_t
  type dense_solver_t

     type(dense_solver_cp_t),       pointer :: cp       => NULL()
     type(dense_solver_lapack_t),   pointer :: lapack   => NULL()

  endtype dense_solver_t

  
  !
  ! Interface definition
  !

  public :: init
  interface init
     module procedure dense_solver_init
  endinterface

  public :: del
  interface del
     module procedure dense_solver_del
  endinterface

  public :: diag_start
  interface diag_start
     module procedure dense_solver_diag_start
  endinterface

  public :: diag_stop
  interface diag_stop
     module procedure dense_solver_diag_stop
  endinterface

  public :: diag_HS
  interface diag_HS
     module procedure dense_solver_diag
  endinterface

  public :: e_bs
  interface e_bs
     module procedure dense_solver_e_bs
  endinterface

  public :: mulliken
  interface mulliken
     module procedure dense_solver_mulliken
  endinterface

  public :: get_dict
  interface get_dict
     module procedure dense_solver_get_dict
  endinterface

  public :: register
  interface register
     module procedure dense_solver_register
  endinterface

contains

  integer function l2i(l)
    use, intrinsic :: iso_c_binding
    logical(C_BOOL) :: l
    if (l) then
       l2i = 1
    else
       l2i = 0
    endif
  endfunction l2i

  !**********************************************************************
  ! Initialize the solver
  !**********************************************************************
  subroutine dense_solver_init(this, error)
    implicit none

    type(dense_solver_t),           intent(inout)  :: this
    integer,              optional, intent(out)    :: error

    ! ---

    integer :: n

    ! ---

    INIT_ERROR(error)

    if (associated(this%cp) .and. &
         associated(this%lapack)) then
       n = sum( [ l2i(this%cp%enabled), &
                  l2i(this%lapack%enabled) ] )
       if (n > 1) then
          RAISE_ERROR("Please specify only a single solver for the dense NOTB module.", error)
       endif

       ! If no solver is given, select LAPACK as the default
       if (n == 0) then
          this%lapack%enabled = .true.
       endif

       if (.not. this%cp%enabled)        deallocate(this%cp)
       if (.not. this%lapack%enabled)    deallocate(this%lapack)
    endif

    if (associated(this%cp)) then
       call init(this%cp, error=error)
       PASS_ERROR(error)
    endif
    if (associated(this%lapack)) then
       call init(this%lapack, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_solver_init


  !**********************************************************************
  ! Delete the solver
  !**********************************************************************
  subroutine dense_solver_del(this)
    implicit none

    type(dense_solver_t), intent(inout)  :: this

    ! ---

    if (associated(this%cp)) then
       call del(this%cp)
       deallocate(this%cp)
    endif
    if (associated(this%lapack)) then
       call del(this%lapack)
       deallocate(this%lapack)
    endif

  endsubroutine dense_solver_del


  !**********************************************************************
  ! Diagonalize the system and determine the density matrix
  !**********************************************************************
  subroutine dense_solver_diag_start(this, tb, error)
    implicit none

    type(dense_solver_t),                intent(inout)  :: this
    type(dense_hamiltonian_t), target,   intent(inout)  :: tb
    integer,                   optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    if (associated(this%cp)) then
       call diag_start(this%cp, tb, error=error)
       PASS_ERROR(error)
    endif
    if (associated(this%lapack)) then
       call diag_start(this%lapack, tb, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_solver_diag_start


  !**********************************************************************
  ! Diagonalize the system and determine the density matrix
  !**********************************************************************
  subroutine dense_solver_diag_stop(this, tb, error)
    implicit none

    type(dense_solver_t),                intent(inout)  :: this
    type(dense_hamiltonian_t), target,   intent(inout)  :: tb
    integer,                   optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    if (associated(this%cp)) then
       call diag_stop(this%cp, tb, error=error)
       PASS_ERROR(error)
    endif
    if (associated(this%lapack)) then
       call diag_stop(this%lapack, tb, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_solver_diag_stop


  !**********************************************************************
  ! Diagonalize the system and determine the density matrix
  !**********************************************************************
  subroutine dense_solver_diag(this, tb, noc, phi, error)
    implicit none

    type(dense_solver_t),      intent(inout) :: this
    type(dense_hamiltonian_t), intent(inout) :: tb
    real(DP),                  intent(in)    :: noc
    real(DP),        optional, target        :: phi(*)
    integer,         optional, intent(out)   :: error

    ! ---

    INIT_ERROR(error)

    if (associated(this%cp)) then
       call diag_HS(this%cp, tb, noc, phi, error=error)
       PASS_ERROR(error)
    endif
    if (associated(this%lapack)) then
       call diag_HS(this%lapack, tb, noc, phi, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_solver_diag


  !>
  !! Mulliken charges
  !!
  !! Mulliken charges. Updates the charges of atoms that
  !! are treated by tight-binding.
  !<
  subroutine dense_solver_mulliken(this, tb, q, error)
    implicit none

    type(dense_solver_t),                 intent(inout)  :: this    !< Solver object
    type(dense_hamiltonian_t),            intent(in)     :: tb
    real(DP),                             intent(out)    :: q(:)    !< Charges
    integer,                    optional, intent(out)    :: error


    ! ---

    INIT_ERROR(error)

    if (associated(this%cp)) then
       call mulliken(this%cp, tb, q, error=error)
       PASS_ERROR(error)
    endif
    if (associated(this%lapack)) then
       call mulliken(this%lapack, tb, q, error=error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_solver_mulliken


  !**********************************************************************
  ! The band-structure energy
  !**********************************************************************
  function dense_solver_e_bs(this, tb, error)
    implicit none

    type(dense_solver_t),      intent(inout) :: this
    type(dense_hamiltonian_t), intent(in)    :: tb
    integer,         optional, intent(out)   :: error

    real(DP)                                 :: dense_solver_e_bs

    ! ---

    INIT_ERROR(error)

    if (associated(this%cp)) then
       dense_solver_e_bs = e_bs(this%cp, tb, error=error)
       PASS_ERROR(error)
    endif
    if (associated(this%lapack)) then
       dense_solver_e_bs = e_bs(this%lapack, tb, error=error)
       PASS_ERROR(error)
    endif

  endfunction dense_solver_e_bs


  !>
  !! Return dictionary object containing pointers to internal data
  !<
  subroutine dense_solver_get_dict(this, dict, error)
    implicit none

    type(dense_solver_t), intent(inout) :: this        !< NOTB object
    type(ptrdict_t),      intent(inout) :: dict
    integer,    optional, intent(out)   :: error       !< Error signals

    ! ---

    INIT_ERROR(error)

    if (associated(this%lapack)) then
       call get_dict(this%lapack, dict, error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_solver_get_dict


  !>
  !! Register solver dispatch module and all solvers
  !<
  subroutine dense_solver_register(this, cfg)
    implicit none

    type(dense_solver_t), target     :: this
    type(c_ptr),          intent(in) :: cfg

    ! ---

    allocate(this%cp)
    allocate(this%lapack)

    call register(this%cp, cfg)
    call register(this%lapack, cfg)

  endsubroutine dense_solver_register

endmodule dense_solver
