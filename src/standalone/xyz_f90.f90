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
#include "macros.inc"

!>
!! XYZ output module
!!
!! The XYZ output module.
!<
module xyz
  use libAtoms_module

  use io

  use particles

  interface read_xyz
     module procedure read_xyz_un, read_xyz_fn
  endinterface

  interface write_xyz
     module procedure write_xyz_un, write_xyz_fn
  endinterface

contains

  !>
  !! Reax XYZ file
  !!
  !! Reax XYZ file from unit.
  !<
  subroutine read_xyz_un(p, un, error)
    implicit none

    type(particles_t), intent(inout)  :: p
    integer, intent(in)               :: un
    integer, intent(inout), optional  :: error

    ! ---

    integer          :: i, nat, Z

    ! ---

    if (.not. initialized(p)) then
       call init(p)
    endif

    read (un, *)  nat
    read (un, *)

    call allocate(p, nat)

    do i = 1, p%nat
       read (un, *)  p%sym(i), PNC3(p, i)
#ifndef IMPLICIT_R
       POS3(p, i)  = PNC3(p, i)
#endif

       p%index(i) = i

       Z = atomic_number(p%sym(i))
       if (Z > 0 .and. Z <= MAX_Z) then
          p%Z(i)  = Z
          p%m(i)  = ElementMass(Z)
       else
          RAISE_ERROR("Unknown chemical symbol '" // trim(p%sym(i)) // "' encountered.", error)
       endif

    enddo

    call update_elements(p)

  endsubroutine read_xyz_un


  !>
  !! Reax XYZ file
  !!
  !! Reax XYZ file from named file.
  !<
  subroutine read_xyz_fn(p, fn, error)
    implicit none

    type(particles_t), intent(inout)  :: p
    character(*), intent(in)          :: fn
    integer, intent(inout), optional  :: error

    ! ---

    integer  :: un

    ! ---

    un  = fopen(fn, F_READ)
    call read_xyz_un(p, un, error)
    call fclose(un)
    PASS_ERROR_WITH_INFO("Filename '" // trim(fn) // "'.", error)

  endsubroutine read_xyz_fn


  !>
  !! Write an XYZ file.
  !!
  !! Write an XYZ file.
  !<
  subroutine write_xyz_un(un, p, conv, q)
    implicit none

    integer, intent(in)             :: un               !< Unit number
    type(particles_t), intent(in)   :: p                !< Particles object
    real(DP), intent(in), optional  :: conv             !< Scaling factor for lengths
    real(DP), intent(in), dimension(:), optional :: q   !< Additional array (5th column of xyz-file)

    ! ---

    integer   :: i
    real(DP)  :: conv_loc

    ! ---

    if(present(conv)) then
       conv_loc = conv
    else
       conv_loc = 1.0_DP
    endif


    write (un,*) p%nat
    write (un,*) 'Cluster'
    if(.not. present(q)) then
       do i = 1, p%nat
          write (un,'(a4,3f11.3)')  p%sym(i), POS3(p, i)*conv_loc
       enddo
    else
       do i = 1, p%nat
          write (un,'(a4,4f11.3)')  p%sym(i), POS3(p, i)*conv_loc, q(i)
       enddo
    end if

  endsubroutine write_xyz_un


  !>
  !! Write an XYZ file.
  !!
  !! Write an XYZ file.
  !<
  subroutine write_xyz_fn(fn, p, conv, q)
    implicit none

    character(*), intent(in)       :: fn                !< File name
    type(particles_t), intent(in)  :: p                 !< Particles object
    real(DP), intent(in), optional :: conv              !< Scaling factor for lengths
    real(DP), intent(in), dimension(:), optional :: q   !< Additional array (5th column of xyz-file)

    ! ---

    integer  :: un

    ! ---

    un  = fopen(fn, F_WRITE)
    call write_xyz_un(un, p, conv, q)
    call fclose(un)

  endsubroutine write_xyz_fn
  
endmodule xyz
