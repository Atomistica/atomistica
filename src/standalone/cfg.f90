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
!! CFG input/output module (AtomEye)
!!
!! Write particle configurations into a CFG (AtomEye) format. This format
!! supports additional arbitrary atomic properties which enables output
!! of, e.g., charges.
!!
!! See
!! http://mt.seas.upenn.edu/Archive/Graphics/A/
!<
module cfg
  use libAtoms_module

  use io
  use particles

  interface read_cfg
     module procedure read_cfg_un, read_cfg_fn
  endinterface

  interface write_cfg
     module procedure write_cfg_un, write_cfg_fn
  endinterface

contains

  !>
  !! Read CFG from file unit
  !!
  !! Reads a CFG file provided an already open unit \param un is provided.
  !<
  subroutine read_cfg_un(p, un, error)
    implicit none

    type(particles_t), intent(inout)  :: p       !> Particle configuration
    integer, intent(in)               :: un      !> File unit
    integer, intent(inout), optional  :: error  !> Error passing

    ! ---

    integer          :: i, j, nat, Z
    character(80)    :: dummy1, dummy2, dummy3, dummy4
    real(DP)         :: x(3)

    ! ---

    if (.not. initialized(p)) then
       call init(p)
    endif

    read (un, *)  dummy1, dummy2, dummy3, dummy4, nat
    read (un, *)

    call allocate(p, nat)

    do i = 1, 3
       do j = 1, 3
          read (un, *)  dummy1, dummy2, dummy3, p%Abox(i, j)
       enddo
    enddo

    call set_cell(p, p%Abox(:, :), error=error)
    PASS_ERROR(error)

    read (un, *) 
    read (un, *)  dummy1, dummy2, i
    do j = 1, i-3
       read (un, *)
    enddo

    do i = 1, p%nat
       read (un, *)  p%m(i)
       read (un, *)  p%sym(i)
       read (un, *)  x(:)
       PNC3(p, i)         = matmul(p%Abox(:, :), x(:))
#ifndef IMPLICIT_R
       POS3(p, i)         = PNC3(p, i)
#endif
       VEC3(p%r_cont, i)  = PNC3(p, i)
       p%g(i) = 0

       p%index(i) = i
       Z = atomic_number(p%sym(i))
       if (Z > 0 .and. Z < MAX_Z) then
          p%Z(j) = Z
       endif

    enddo

    call update_elements(p)

  endsubroutine read_cfg_un


  !>
  !! Read CFG from file name
  !!
  !! Reads a CFG file from a given file name \param fn.
  !<
  subroutine read_cfg_fn(p, fn, error)
    implicit none

    type(particles_t), intent(inout)  :: p       !> Particle configuration
    character(*), intent(in)          :: fn      !> File name
    integer, intent(inout), optional  :: error  !> Error passing

    ! ---

    integer  :: un

    ! ---

    un  = fopen(fn, F_READ)
    call read_cfg_un(p, un, error)
    call fclose(un)
    PASS_ERROR_WITH_INFO("Filename '" // trim(fn) // "'.", error)

  endsubroutine read_cfg_fn


  !>
  !! Write CFG to file unit
  !!
  !! Write a CFG file provided an already open unit \param un (in write mode) is provided.
  !<
  subroutine write_cfg_un(un, p, conv_in, error)
    implicit none

    integer, intent(in)               :: un       !> File unit
    type(particles_t), intent(in)     :: p        !> Particle configuration
    real(DP), intent(in), optional    :: conv_in  !> Length conversion factor
    integer, intent(inout), optional  :: error   !> Error passing

    ! ---

    integer, parameter  :: MAX_AUX = 100

    integer             :: i, j, k, naux
    character(80)       :: fmt

    real(DP)            :: conv
    real(DP)            :: aux(MAX_AUX)

    ! ---

    conv  = 1.0_DP
    if (present(conv_in)) then
       conv  = conv_in
    endif

    write (un, '(A,I10)')  "Number of particles = ", p%nat
    write (un, '(A)')      "A = 1.0 Angstrom"

    do i = 1, 3
       do j = 1, 3
          write (un, '(A,I1.1,A,I1.1,A,ES16.9,A)')  "H0(", i, ",", j, ") = ", p%Abox(j, i)*conv, " A"
       enddo
    enddo

    naux = p%data%n_real + p%data%n_integer + 3*p%data%n_real3

    if (naux > MAX_AUX) then
       RAISE_ERROR("Too many auxiliary properties.", error)
    endif

    write (un, '(A)')     ".NO_VELOCITY."
    write (un, '(A,I5)')  "entry_count = ", 3+naux
    do i = 1, p%data%n_real
       write (un, '(A,I2.2,A,A)')     "auxiliary[", i-1, "] = ", trim(p%data%name_real(i))
    enddo

    k  = p%data%n_real
    do i = 1, p%data%n_integer
       write (un, '(A,I2.2,A,A)')     "auxiliary[", k+i-1, "] = ", trim(p%data%name_integer(i))
    enddo

    k  = p%data%n_real + p%data%n_integer
    do i = 1, p%data%n_real3
       write (un, '(A,I2.2,A,A)')     "auxiliary[", k+3*(i-1), "] = ", trim(p%data%name_real3(i)) // "_x"
       write (un, '(A,I2.2,A,A)')     "auxiliary[", k+3*(i-1)+1, "] = ", trim(p%data%name_real3(i)) // "_y"
       write (un, '(A,I2.2,A,A)')     "auxiliary[", k+3*(i-1)+2, "] = ", trim(p%data%name_real3(i)) // "_z"
    enddo

    write (fmt, '(A,I2.2,A)')  "(", 3+naux, "(ES16.9,1X))"

    do i = 1, p%nat
       write (un, '(ES16.9)')   p%m(i)
       write (un, '(A)')        p%sym(i)
       if (naux > 0) then
          do j = 1, p%data%n_real
             aux(j)  = p%data%data_real(i, j)
          enddo

          k  = p%data%n_real
          do j = 1, p%data%n_integer
             aux(j+k)  = p%data%data_integer(i, j)
          enddo

          k  = p%data%n_real + p%data%n_integer
          do j = 1, p%data%n_real3
#ifdef SEP_XYZ
             aux(3*(j-1)+k+1)  = p%data%data_real3(i, 1, j)
             aux(3*(j-1)+k+2)  = p%data%data_real3(i, 2, j)
             aux(3*(j-1)+k+3)  = p%data%data_real3(i, 3, j)
#else
             aux(3*(j-1)+k+1)  = p%data%data_real3(1, i, j)
             aux(3*(j-1)+k+2)  = p%data%data_real3(2, i, j)
             aux(3*(j-1)+k+3)  = p%data%data_real3(3, i, j)
#endif
          enddo

          write (un, fmt) &
               matmul(p%Bbox, POS3(p, i)), &
               conv*aux(1:naux)
       else
          write (un, fmt) &
               matmul(p%Bbox, POS3(p, i))
       endif
    enddo

  endsubroutine write_cfg_un



  !>
  !! Write CFG to file name
  !!
  !! Writes a CFG file to a file given by the name \param fn.
  !<
  subroutine write_cfg_fn(fn, p, conv_in, error)
    implicit none

    character(*), intent(in)          :: fn       !> File name
    type(particles_t), intent(in)     :: p        !> Particle configuration
    real(DP), intent(in), optional    :: conv_in  !> Length conversion factor
    integer, intent(inout), optional  :: error   !> Error passing

    ! ---

    integer  :: un

    ! ---

    un  = fopen(fn, F_WRITE)
    call write_cfg_un(un, p, conv_in, error=error)
    PASS_ERROR(error)
    call fclose(un)

  endsubroutine write_cfg_fn

endmodule cfg
