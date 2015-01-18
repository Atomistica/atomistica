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
!**********************************************************************
! Thread local storage: Every thread has a set of pointers for local
! reduction operations. These are shared among different modules
! to save memory overhead.
!
! Thanks Jim Dempsey at the Intel support forum for pointing out
! this strategy.
!**********************************************************************

#include "macros.inc"

module tls
  use error_module
  use system_module

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  private

  public :: tls_sca1, tls_vec1, tls_mat1, tls_mat2
  real(DP), allocatable, save  :: tls_sca1(:)
  real(DP), allocatable, save  :: tls_vec1(:, :)
  real(DP), allocatable, save  :: tls_mat1(:, :)
  real(DP), allocatable, save  :: tls_mat2(:, :)
  !$omp threadprivate(tls_sca1, tls_vec1, tls_mat1, tls_mat2)

  public :: tls_init, tls_del, tls_reduce

contains

  !**********************************************************************
  ! Allocate TLS
  !**********************************************************************
  subroutine tls_init(n, sca, vec, mat, ierror)
    implicit none

    integer, intent(in)               :: n
    integer, intent(in), optional     :: sca
    integer, intent(in), optional     :: vec
    integer, intent(in), optional     :: mat
    integer, intent(inout), optional  :: ierror

    ! ---

    if (present(sca)) then

       if (sca <= 1) then

          if (allocated(tls_sca1) .and. size(tls_sca1) /= n) then
             deallocate(tls_sca1)
          endif

          if (.not. allocated(tls_sca1)) then
             allocate(tls_sca1(n))
          endif

          tls_sca1  = 0.0_DP
          
       else

          RAISE_ERROR("Only up to 1 scalar array supported right now.", ierror)

       endif

    endif


    if (present(vec)) then
    
       if (vec <= 1) then

          if (allocated(tls_vec1) .and. size(tls_vec1) /= 3*n) then
             deallocate(tls_vec1)
          endif

          if (.not. allocated(tls_vec1)) then
             allocate(tls_vec1(3, n))
          endif
    
          tls_vec1  = 0.0_DP

       else

          RAISE_ERROR("Only up to 1 vector array supported right now.", ierror)

       endif

    endif


    if (present(mat)) then
    
       if (mat <= 2) then

          if (allocated(tls_mat1) .and. size(tls_mat1) /= n*n) then
             deallocate(tls_mat1)
          endif

          if (.not. allocated(tls_mat1)) then
             allocate(tls_mat1(n, n))
          endif
    
          tls_mat1  = 0.0_DP

          if (mat == 2) then

             if (allocated(tls_mat2) .and. size(tls_mat2) /= n*n) then
                deallocate(tls_mat2)
             endif

             if (.not. allocated(tls_mat2)) then
                allocate(tls_mat2(n, n))
             endif

             tls_mat2  = 0.0_DP

          endif

       else

          RAISE_ERROR("Only up to 2 matrices supported right now.", ierror)

       endif

    endif

  endsubroutine tls_init


  !**********************************************************************
  ! Deallocate TLS
  !**********************************************************************
  subroutine tls_del
    implicit none

    ! ---

    if (allocated(tls_sca1)) then
       deallocate(tls_sca1)
    endif

    if (allocated(tls_vec1)) then
       deallocate(tls_vec1)
    endif

    if (allocated(tls_mat1)) then
       deallocate(tls_mat1)
    endif

    if (allocated(tls_mat2)) then
       deallocate(tls_mat2)
    endif
    
  endsubroutine tls_del


  !**********************************************************************
  ! Perform a reduction operation on something stored in a thread
  ! local storage.
  !**********************************************************************
  subroutine tls_reduce(n, sca1, vec1, mat1, mat2)
    implicit none

    integer, intent(in)                :: n
    real(DP), intent(inout), optional  :: sca1(:)
    real(DP), intent(inout), optional  :: vec1(:, :)
    real(DP), intent(inout), optional  :: mat1(:, :)
    real(DP), intent(inout), optional  :: mat2(:, :)

    ! ---

#ifdef _OPENMP

    integer  :: i, j, j1, j2, dn, tnum, numt

    ! ---

    numt  = omp_get_num_threads()
    tnum  = omp_get_thread_num()

    dn    = n/numt+1

    if (present(sca1)) then
       do i = 0, numt-1
          !$omp barrier
          j   = mod(tnum + i, numt)
          j1  = j*dn+1
          j2  = min((j+1)*dn, n)
          sca1(j1:j2)  = sca1(j1:j2) + tls_sca1(j1:j2)
       enddo
    endif

    if (present(vec1)) then
       do i = 0, numt-1
          !$omp barrier
          j   = mod(tnum + i, numt)
          j1  = j*dn+1
          j2  = min((j+1)*dn, n)
          VEC3(vec1, j1:j2)  = VEC3(vec1, j1:j2) + VEC3(tls_vec1, j1:j2)
       enddo
    endif

    if (present(mat1)) then
       do i = 0, numt-1
          !$omp barrier
          j   = mod(tnum + i, numt)
          j1  = j*dn+1
          j2  = min((j+1)*dn, n)
          mat1(1:n, j1:j2)  = mat1(1:n, j1:j2) + tls_mat1(1:n, j1:j2)
       enddo
    endif

    if (present(mat2)) then
       do i = 0, numt-1
          !$omp barrier
          j   = mod(tnum + i, numt)
          j1  = j*dn+1
          j2  = min((j+1)*dn, n)
          mat2(1:n, j1:j2)  = mat2(1:n, j1:j2) + tls_mat2(1:n, j1:j2)
       enddo
    endif

#else

    if (present(sca1)) then
       sca1(1:n)  = sca1(1:n) + tls_sca1(1:n)
    endif

    if (present(vec1)) then
       VEC3(vec1, 1:n)  = VEC3(vec1, 1:n) + VEC3(tls_vec1, 1:n)
    endif

    if (present(mat1)) then
       mat1(1:n, 1:n)  = mat1(1:n, 1:n) + tls_mat1(1:n, 1:n)
    endif

    if (present(mat2)) then
       mat2(1:n, 1:n)  = mat2(1:n, 1:n) + tls_mat2(1:n, 1:n)
    endif

#endif

    !$omp barrier

  endsubroutine tls_reduce

endmodule tls
