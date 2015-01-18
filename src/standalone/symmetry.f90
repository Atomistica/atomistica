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
!! Treatment of crystal symmetries
!<

#include "macros.inc"

module symmetry
  use libAtoms_module

  use logging

  use particles
  use cyclic

  implicit none

  private

  public :: symmetry_t
  type symmetry_t

     integer                :: n            ! Number of operations
     integer, allocatable   :: sym(:)       ! Pointer to the actual operation
     real(DP), allocatable  :: t(:, :, :)   ! Generator of the operation in real space
     real(DP), allocatable  :: u(:, :, :)   ! Generator of the operation in reciprocal space

  endtype symmetry_t

  integer, parameter  :: SYM_IDENTITY   = 1
  integer, parameter  :: SYM_INVERSION  = 13

  public :: SYM_IDENTITY, SYM_INVERSION

  integer, save                    :: n            ! Number of operations
  character(8), save, allocatable  :: id(:)        ! Associated names
  real(DP), save, allocatable      :: t(:, :, :)   ! Generator


  public :: del
  interface del
     module procedure symmetry_del
  endinterface

  real(DP), parameter :: EPS = 1d-10

  public :: symmetry_init, symmetry_analysis, symmetry_check_reciprocal

contains

  !>
  !! Generate a list of symmetry operations
  !! Adopted from crystsym.f of the MBPP code
  !<
  subroutine symmetry_init
    implicit none

    real(DP)  :: f
    integer   :: i, j

    ! ---

    ! For now, we only cover the hexagonal group

    if (.not. allocated(id)) then

       n = 24

       allocate(id(n), t(3, 3, n))

       t(1:3,1:3,1:n) = 0.D0

       f = 0.5D0*sqrt(3.D0)

       t(1,1,2) = 0.5D0
       t(1,2,2) = -f
       t(2,1,2) = f
       t(2,2,2) = 0.5D0
       t(1,1,7) = -0.5D0
       t(1,2,7) = -f
       t(2,1,7) = -f
       t(2,2,7) = 0.5D0

       t(3,3,1: 6) = 1.D0
       t(3,3,7:12) = -1.D0

       do i = 1, 2
          t(i,i,1) = 1.D0
          do j = 1, 2
             t(i,j, 6) = t(j,i, 2)
             t(i,j, 3) = sum(t(i,1:2,2)*t(1:2,j,2))
             t(i,j, 8) = sum(t(i,1:2,2)*t(1:2,j,7))
             t(i,j,12) = sum(t(i,1:2,7)*t(1:2,j,2))
          enddo
       enddo
       do i = 1, 2
          do j = 1, 2
             t(i,j, 5) = t(j,i, 3)
             t(i,j, 4) = sum(t(i,1:2, 2)*t(1:2,j,3))
             t(i,j, 9) = sum(t(i,1:2, 2)*t(1:2,j,8))
             t(i,j,10) = sum(t(i,1:2,12)*t(1:2,j,3))
             t(i,j,11) = sum(t(i,1:2,12)*t(1:2,j,2))
          enddo
       enddo

       t(1:3,1:3,13:24) = -t(1:3,1:3,1:12)

       id(1:12) = &
            (/'E       ','+C6z    ','+C3z    ','C2z     ', &
            '-C3z    ','-C6z    ','C2      ','C2      ', &
            'C2x     ','C2      ','C2      ','C2y     '/)

       id(13:24)(1:1) = 'I'
       id(13:24)(2:8) = id(1:12)(1:7)

    endif

  endsubroutine symmetry_init

  
  !>
  !! Analyse the symmetry of this crystal structure
  !<
  subroutine symmetry_analysis(s, p, A, B)
    implicit none

    type(symmetry_t), intent(inout)  :: s
    type(particles_t), intent(in)    :: p
    real(DP), intent(in)             :: A(3, 3)
    real(DP), intent(in)             :: B(3, 3)

    ! ---

    integer   :: i, j, k

    integer   :: nsymlat
    integer   :: symlat(n)
    real(DP)  :: u(3, 3, n)

    real(DP)  :: r(p%nat, 3), sym_r(p%nat, 3)

    real(DP)  :: help(3, 3)
    integer   :: ipiv(3), info

    ! ---

    call del(s)

    allocate(s%sym(n), s%t(3, 3, n), s%u(3, 3, n))

    ! Determine the point group of the LATTICE

    nsymlat = 0
    do i = 1, n
       u(:, :, i) = matmul(B, matmul(t(:, :, i), A))

       if (all(abs(u(:, :, i) - nint(u(:, :, i))) < EPS)) then
          ! This operation is commensurable with our lattice
          nsymlat          = nsymlat + 1
          symlat(nsymlat)  = i
       endif
    enddo

    ! --- Output ---
    write (ilog, '(1X, A)')  "The lattice is commensurable with the following operations:"
    write (ilog, *)  (id(symlat(i)), i = 1, nsymlat)
    write (ilog, *)
    ! --------------


    ! Determine the point group of the CRYSTAL

    s%n = 0
    do i = 1, nsymlat
       do j = 1, p%nat
          r(j, :)      = in_bounds(p, POS3(p, j))
          sym_r(j, :)  = in_bounds(p, matmul(t(:, :, symlat(i)), r(j, :)))
       enddo

       if (identical(p, p%nat, r, sym_r)) then
          s%n         = s%n+1
          s%sym(s%n)  = symlat(i)

          s%t(:, :, s%n)  = t(:, :, symlat(i))
          help            = t(:, :, symlat(i))
          s%u(:, :, s%n)  = 0
          do k = 1, 3
             s%u(k, k, s%n) = 1
          enddo
          call dgesv(3, 3, help, 3, ipiv, s%u(:, :, s%n), 3, info)

          if (info /= 0) then
             stop '[symmetry_analysis] Error inverting symmetry operation.'
          endif
       endif
    enddo

    ! --- Output ---
    write (ilog, '(1X, A)')  "The crystal is commensurable with the following operations:"
    write (ilog, *)  (id(s%sym(i)), i = 1, s%n)
    write (ilog, *)
    ! --------------

    s%n = 2
    s%sym(2) = s%sym(SYM_INVERSION)
    s%u(:, :, 2) = s%u(:, :, SYM_INVERSION)
    s%n = 1

  endsubroutine symmetry_analysis


  !>
  !! Destructor
  !! 
  !! Free memory
  !<
  subroutine symmetry_del(this)
    implicit none

    type(symmetry_t)  :: this

    ! ---

    if (allocated(this%sym))  deallocate(this%sym)
    if (allocated(this%t))    deallocate(this%t)
    if (allocated(this%u))    deallocate(this%u)
    
  endsubroutine symmetry_del


  !>
  !! Check if the list of points does have the given symmetry properties
  !! with respect to the real space
  !<
  function symmetry_check_real(s, p, n, x)
    implicit none

    type(symmetry_t), intent(in)   :: s
    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: n
    real(DP), intent(in)           :: x(n, 3)

    logical                        :: symmetry_check_real

    ! ---

    integer   :: i, j
    real(DP)  :: org_x(n, 3), sym_x(n, 3)

    ! ---

    symmetry_check_real = .true.
    do i = 1, s%n
       do j = 1, n
          org_x(j, :)  = in_bounds(p, x(j, :))
          sym_x(j, :)  = in_bounds(p, matmul(t(:, :, s%sym(i)), org_x(j, :)))
       enddo

       if (.not. identical(p, n, org_x, sym_x)) then
          symmetry_check_real = .false.
       endif
    enddo
    
  endfunction symmetry_check_real


  !>
  !! Check if the list of points does have the given symmetry properties
  !! with respect to the reciprocal space
  !<
  function symmetry_check_reciprocal(s, p, n, x)
    implicit none

    type(symmetry_t), intent(in)   :: s
    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: n
    real(DP), intent(in)           :: x(n, 3)

    logical                        :: symmetry_check_reciprocal

    ! ---

    integer   :: i, j
    real(DP)  :: org_x(n, 3), sym_x(n, 3)

    ! ---

    symmetry_check_reciprocal = .true.
    do i = 1, s%n
       do j = 1, n
          org_x(j, :)  = cyclic_in_reciprocal_bounds(p, x(j, :))
          sym_x(j, :)  = cyclic_in_reciprocal_bounds(p, matmul(s%u(:, :, i), org_x(j, :)))
       enddo

       if (.not. identical_reciprocal(p, n, org_x, sym_x)) then
          write (*, *)  id(s%sym(i))
          write (*, '(3F10.5)')  s%u(:, :, i)
          write (*, *)  '---'
          write (*, '(3F10.5)')  (org_x(j, :), j = 1, n)
          write (*, *)  '---'
          write (*, '(3F10.5)')  (sym_x(j, :), j = 1, n)
          write (*, *)  '---'

          symmetry_check_reciprocal = .false.
       endif
    enddo
    
  endfunction symmetry_check_reciprocal


  !>
  !! Check if the two structures are the same
  !<
  function identical(p, nat, r1, r2)
    implicit none

    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: nat
    real(DP), intent(in)           :: r1(nat, 3)
    real(DP), intent(in)           :: r2(nat, 3)

    logical                        :: identical

    ! ---

    integer   :: i, j
    real(DP)  :: n, d, mindist(nat)
    real(DP)  :: dr(3)

    ! ---

    do i = 1, nat
       d = 10*EPS
       do j = 1, nat
          dr = in_bounds(p, r1(i, :) - r2(j, :))
          n = sqrt(dot_product(dr, dr))

!          write (*, *)  r1(i, :)
!          write (*, *)  r2(j, :)
!          write (*, *)  n
!          write (*, *)  '...'          

          if (n < d) then
             d = n
          endif
       enddo

       mindist(i) = d
    enddo

    identical = all(mindist < EPS)

  endfunction identical


  !>
  !! Check if the two structures are the same
  !<
  function identical_reciprocal(p, nat, r1, r2) result(identical)
    implicit none

    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: nat
    real(DP), intent(in)           :: r1(nat, 3)
    real(DP), intent(in)           :: r2(nat, 3)

    logical                        :: identical

    ! ---

    integer   :: i, j
    real(DP)  :: n, d, mindist(nat)
    real(DP)  :: dr(3)

    ! ---

    do i = 1, nat
       d = 10*EPS
       do j = 1, nat
          dr = cyclic_in_reciprocal_bounds(p, r1(i, :) - r2(j, :))
          n = sqrt(dot_product(dr, dr))

!         write (*, *)  r1(i, :)
!         write (*, *)  r2(j, :)
!         write (*, *)  n
!         write (*, *)  '...'          

          if (n < d) then
             d = n
          endif
       enddo

       mindist(i) = d
    enddo

    identical = all(mindist < EPS)

  endfunction identical_reciprocal

endmodule symmetry
