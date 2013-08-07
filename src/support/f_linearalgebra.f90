!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
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

!>
!! Matrix functions and helpers
!!
!! Matrix functions and helpers
!<
module linearalgebra
  use system_module

  implicit none

  private

  public :: tr
  interface tr
     module procedure dtr
     module procedure ztr
  endinterface

  public :: multr
  interface multr
     module procedure dmultr
     module procedure zmultr
  endinterface

  public :: identity
  interface identity
     module procedure didentity
     module procedure zidentity
  endinterface

  public :: normsq
  interface normsq
     module procedure dnormsq
  endinterface normsq

  public :: norm
  interface norm
     module procedure dnorm
  endinterface norm

  !% Return the trace of a matrix.
  public :: trace
  interface trace
     module procedure matrix_trace
  endinterface

  public :: ev_bounds
  interface ev_bounds
     subroutine dev_bounds(n, H, l, u) bind(C)
       use, intrinsic :: iso_c_binding
       integer(C_INT), value       :: n
       real(C_DOUBLE), intent(in)  :: H(n, n)
       real(C_DOUBLE), intent(out) :: l
       real(C_DOUBLE), intent(out) :: u
     endsubroutine dev_bounds
     subroutine zev_bounds(n, H, l, u) bind(C)
       use, intrinsic :: iso_c_binding
       integer(C_INT),    value       :: n
       complex(C_DOUBLE), intent(in)  :: H(n, n)
       real(C_DOUBLE),    intent(out) :: l
       real(C_DOUBLE),    intent(out) :: u
     endsubroutine zev_bounds
  endinterface

  public :: iterative_matrix_inverse
  interface
     subroutine iterative_matrix_inverse(mat, invmat, n, prev, epsilon, work1, &
       work2, error, cublas_handle, nit) bind(C)
       use, intrinsic :: iso_c_binding
       real(C_DOUBLE),           intent(in)     :: mat(n, n)     !< Matrix to be inverted
       real(C_DOUBLE),           intent(inout)  :: invmat(n, n)  !< Inverse
       integer(C_INT),           value          :: n             !< Matrix size (n,n)
       logical(C_BOOL),          value          :: prev          !< Previous inverse supplied in invM?
       real(C_DOUBLE),           value          :: epsilon       !< Converge criterion on inverse matrix elements
#ifdef NO_BIND_C_OPTIONAL
       real(C_DOUBLE),           target         :: work1(n, n)   !< Workspace matrix
       real(C_DOUBLE),           target         :: work2(n, n)   !< Workspace matrix
       integer(C_INT),           intent(out)    :: error
       type(C_PTR),              value          :: cublas_handle
       integer(C_INT),           intent(out)    :: nit
#else
       real(C_DOUBLE), optional, target         :: work1(n, n)   !< Workspace matrix
       real(C_DOUBLE), optional, target         :: work2(n, n)   !< Workspace matrix
       integer(C_INT), optional, intent(out)    :: error
       type(C_PTR),    optional, value          :: cublas_handle
       integer(C_INT), optional, intent(out)    :: nit
#endif
     endsubroutine iterative_matrix_inverse
  endinterface

  public :: cross_product

contains

  !>
  !! Trace of a real matrix
  !!
  !! Trace of a real matrix
  !<
  function dtr(n, mat)
    implicit none

    integer, intent(in)   :: n
    real(DP), intent(in)  :: mat(n, n)
    real(DP)              :: dtr

    ! ---

    integer   :: i
    real(DP)  :: t

    ! ---

    t = 0.0
    do i = 1, n
       t = t + mat(i, i)
    enddo
    dtr = t

  endfunction dtr


  !>
  !! Trace of a complex matrix
  !!
  !! Trace of a complex matrix
  !<
  function ztr(n, mat)
    implicit none

    integer, intent(in)      :: n
    complex(DP), intent(in)  :: mat(n, n)
    complex(DP)              :: ztr

    ! ---

    integer      :: i
    complex(DP)  :: t

    ! ---

    t = 0.0
    do i = 1, n
       t = t + mat(i, i)
    enddo
    ztr = t

  endfunction ztr


  !>
  !! Trace of a matrix product
  !!
  !! Take the trace of the product of two matrices [O(N^2) operation]
  !<
  function dmultr(n, mat1, mat2)
    implicit none

    integer, intent(in)   :: n
    real(DP), intent(in)  :: mat1(n, n)
    real(DP), intent(in)  :: mat2(n, n)
    real(DP)              :: dmultr

    ! ---

    dmultr  = sum(transpose(mat1)*mat2)

  endfunction dmultr


  !>
  !! Trace of a matrix product
  !!
  !! Take the trace of the product of two matrices [O(N^2) operation]
  !<
  function zmultr(n, mat1, mat2)
    implicit none

    integer, intent(in)      :: n
    complex(DP), intent(in)  :: mat1(n, n)
    complex(DP), intent(in)  :: mat2(n, n)
    complex(DP)              :: zmultr

    ! ---

    integer      :: i, j
    complex(DP)  :: t

    ! ---

    t = 0.0
    do i = 1, n
       do j = 1, n
          t = t + mat1(i, j)*mat2(j, i)
       enddo
    enddo
    zmultr = t

  endfunction zmultr


  !>
  !! Identity matrix
  !!
  !! Identity matrix
  !<
  subroutine didentity(n, mat)
    implicit none

    integer, intent(in)      :: n
    real(DP), intent(inout)  :: mat(n, n)

    ! ---

    integer   :: i

    ! ---

    mat = 0.0_DP
    do i = 1, n
       mat(i, i) = 1.0_DP
    enddo

  endsubroutine didentity


  !>
  !! Identity matrix
  !!
  !! Identity matrix
  !<
  subroutine zidentity(n, mat)
    implicit none

    integer, intent(in)         :: n
    complex(DP), intent(inout)  :: mat(n, n)

    ! ---

    integer   :: i

    ! ---

    mat = 0.0_DP
    do i = 1, n
       mat(i, i) = 1.0_DP
    enddo

  endsubroutine zidentity


  !>
  !! Add a scalar matrix
  !!
  !! Add a scalar to (the diagonal elements of) a matrix
  !<
  subroutine add_scalar(n, s, mat)
    implicit none

    integer, intent(in)      :: n
    real(DP), intent(in)     :: s
    real(DP), intent(inout)  :: mat(n, n)

    ! ---

    integer   :: i

    ! ---

    do i = 1, n
       mat(i, i) = mat(i, i) + s
    enddo

  endsubroutine add_scalar


  ! normsq()
  ! returns (X.dot.X)
  pure function dnormsq(vector) result(normsq) 

    real(dp), intent(in), dimension(:) :: vector
    real(dp)             :: normsq
   
    normsq = dot_product(vector,vector)

  end function dnormsq

  ! norm()
  ! returns SQRT((X.dot.X))
  pure function dnorm(vector) result(norm)

    real(dp), intent(in),dimension(:) :: vector
    real(dp)             :: norm

    norm = sqrt(dot_product(vector,vector))

  end function dnorm


  pure function cross_product(x,y) ! x ^ y

    real(dp), dimension(3), intent(in):: x,y
    real(dp), dimension(3)            :: cross_product

    cross_product(1) = + ( x(2)*y(3) - x(3)*y(2) )
    cross_product(2) = - ( x(1)*y(3) - x(3)*y(1) )
    cross_product(3) = + ( x(1)*y(2) - x(2)*y(1) )

  end function cross_product


  !returns trace of a matrix
  function matrix_trace(matrix) result(tr)
    real(dp),intent(in), dimension(:,:) ::matrix
    real(dp)::tr
    integer::i,N

    N = min(size(matrix,1),size(matrix,2))
   
    tr=0.0_dp
    do i=1,N
       tr=tr+matrix(i,i)
    end do

  end function matrix_trace

endmodule linearalgebra
