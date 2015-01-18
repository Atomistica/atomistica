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
! Module for export to VTK (Visualization ToolKit, www.vtk.org)
! The files intended to be viewed with ParaView (www.paraview.org).
!**********************************************************************

module vtk
  use supplib

  use particles

  implicit none

contains

  !**********************************************************************
  ! Open a file for writing
  !**********************************************************************
  subroutine vtk_header(un, nat, r)
    implicit none

    integer, intent(in)          :: un
    integer, intent(in)          :: nat
    real(DP), intent(in)          :: r(nat, 3)    
    
    ! ---

    integer                      :: i
    
    ! ---

    write (un, '(A)') '# vtk DataFile Version 2.0'
    write (un, '(A)') 'MDCore'
    ! For now, we only support ASCII VTKs
    write (un, '(A)') 'ASCII'
    ! And, of course, only point clouds
    write (un, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (un, '(A, I10, A)') 'POINTS ', nat, ' double'
    do i = 1, nat
       write (un, '(3F20.10)') r(i, :)
    enddo

!    write (un, '(A, I10, I10)') 'CELLS ', np, 3*np
!    do i = 1, np
!       ! This is the VTK_LINE
!       write (un, '(A, I10, I10)') '2 ', p(i)%i-1, p(i)%j-1
!    enddo
    
!    write (un, '(A, I10)') 'CELL_TYPES', np
!    do i = 1, np
!       ! 3 = VTK_LINE
!       write (un, '(A)')  '3'
!    enddo

!    write (un, *)

  endsubroutine vtk_header


  !**********************************************************************
  ! Open a file for writing
  !**********************************************************************
  subroutine vtk_header_with_cells(un, nat, r, nnmax, seed, last, ne)
    implicit none

    integer, intent(in)          :: un
    integer, intent(in)          :: nat
    real(DP), intent(in)          :: r(nat, 3)
    integer, intent(in)          :: nnmax
    integer, intent(in)          :: seed(nat+1)
    integer, intent(in)          :: last(nat+1)
    integer, intent(in)          :: ne(nnmax)
    
    ! ---

    integer                      :: i, j, ncells
    
    ! ---

    write (un, '(A)') '# vtk DataFile Version 2.0'
    write (un, '(A)') 'MDCore'
    ! For now, we only support ASCII VTKs
    write (un, '(A)') 'ASCII'
    ! And, of course, only point clouds
    write (un, '(A)') 'DATASET UNSTRUCTURED_GRID'

    write (un, '(A, I10, A)') 'POINTS ', nat, ' double'
    do i = 1, nat
       write (un, '(3F20.10)') r(i, :)
    enddo

    ncells = 0
    do i = 1, nat
       do j = seed(i), last(i)
          if (ne(j) > i) then
             ncells = ncells+1
          endif
       enddo
    enddo

    write (un, '(A, I10, I10)') 'CELLS ', ncells, 3*ncells
    do i = 1, nat
       do j = seed(i), last(i)
          if (ne(j) > i) then
          ! This is the VTK_LINE
             write (un, '(A, I10, I10)') '2 ', i-1, ne(j)-1
          endif
       enddo
    enddo
    
    write (un, '(A, I10)') 'CELL_TYPES', ncells
    do i = 1, ncells
       ! 3 = VTK_LINE
       write (un, '(A)')  '3'
    enddo

!    write (un, *)

  endsubroutine vtk_header_with_cells


  !**********************************************************************
  ! Prepare for output of point data (i.e., velocities, ...)
  !**********************************************************************
  subroutine vtk_start_point_data(un, nat)
    implicit none

    integer, intent(in)  :: un
    integer, intent(in)  :: nat
    
    ! ---

    write (un, '(A, I10)') 'POINT_DATA', nat

  endsubroutine vtk_start_point_data


  !**********************************************************************
  ! Prepare for output of cell data (i.e., bonding information, ...)
  !**********************************************************************  
  subroutine vtk_start_cell_data(un, nat, nnmax, seed, last, ne)
    implicit none

    integer, intent(in)   :: un
    integer, intent(in)   :: nat
    integer, intent(in)   :: nnmax
    integer, intent(in)   :: seed(nat+1)
    integer, intent(in)   :: last(nat+1)
    integer, intent(in)   :: ne(nnmax)
    
    ! ---

    integer  :: i, j, ncells

    ! ---

    ncells = 0
    do i = 1, nat
       do j = seed(i), last(i)
          if (ne(j) > i) then
             ncells = ncells+1
          endif
       enddo
    enddo

    write (un, '(A, I10)') 'CELL_DATA', ncells

  endsubroutine vtk_start_cell_data


  !**********************************************************************
  ! Write a list of scalars
  !**********************************************************************
  subroutine vtk_write_scalars(un, name, nat, l)
    implicit none

    integer, intent(in)       :: un
    character(*), intent(in)  :: name
    integer, intent(in)       :: nat
    real(DP), intent(in)       :: l(nat)
    
    ! ---

    integer                   :: i
    
    ! ---

    write (un, '(A, A, A)') 'SCALARS ', name, ' double 1'
    write (un, '(A)') 'LOOKUP_TABLE default'

    do i = 1, nat
       write (un, '(F20.10)') l(i)
    enddo

!    write (un, *)

  endsubroutine vtk_write_scalars


  !**********************************************************************
  ! Write a list of scalars
  !**********************************************************************
  subroutine vtk_write_scalars_for_cell(un, name, nat, nn, seed, last, ne, l)
    implicit none

    integer, intent(in)       :: un
    character(*), intent(in)  :: name
    integer, intent(in)       :: nat
    integer, intent(in)       :: nn
    integer, intent(in)       :: seed(nat+1)
    integer, intent(in)       :: last(nat+1)
    integer, intent(in)       :: ne(nn)
    real(DP), intent(in)      :: l(nn)
    
    ! ---

    integer                   :: i, j
    
    ! ---

    write (un, '(A, A, A)') 'SCALARS ', name, ' double 1'
    write (un, '(A)') 'LOOKUP_TABLE default'

    do i = 1, nat
       do j = seed(i), last(i)
          if (ne(j) > i) then
             write (un, '(F20.10)') l(j)
          endif
       enddo
    enddo

!    write (un, *)

  endsubroutine vtk_write_scalars_for_cell


  !**********************************************************************
  ! Write a list of vectors
  !**********************************************************************
  subroutine vtk_write_vectors(un, name, nat, l)
    implicit none

    integer, intent(in)       :: un
    character(*), intent(in)  :: name
    integer, intent(in)       :: nat
    real(DP), intent(in)      :: l(nat, 3)
    
    ! ---

    integer                   :: i
    
    ! ---

    write (un, '(A, A, A)') 'VECTORS ', name, ' double'

    do i = 1, nat
       write (un, '(3F20.10)') l(i, :)
    enddo

!    write (un, *)

  endsubroutine vtk_write_vectors

endmodule vtk
