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
!****************************************************************
! The PDB output module
!****************************************************************

#include "macros.inc"

module pdb
  use libAtoms_module

  use io
  use particles

  interface write_pdb
     module procedure write_pdb_un, write_pdb_fn
  endinterface

contains

  !**********************************************************************
  ! Append a configuration to a PDB trajectory file
  !**********************************************************************    
  subroutine write_pdb_un(un, p, beta, conv, ierror)
    implicit none

    integer, intent(in)               :: un
    type(particles_t), intent(in)     :: p
    real(DP), intent(in), optional    :: beta(:)
    real(DP), intent(in), optional    :: conv
    integer, intent(inout), optional  :: ierror

    ! ---

    character(80)    :: fmt_pdb

    character*6      :: rec_name       ! record name
    character*4      :: atom_name      ! atom name
    character*1      :: altloc         ! alternate location indicator
    character*3      :: res_name       ! residue name
    character*1      :: chainid       ! Chain identifier

    integer          :: resseq             ! Residue sequence number
    character*1      :: icode         ! Code for insertion of residues
    real*8           :: occupancy           ! Occupancy
    character*4      :: segid          ! Segment identifier, left-justified
    character*2      :: element        ! Elementsymbol, right-justified
    character*2      :: charge         ! Charge on the atom
    character*6      :: endtoken

    real*8           :: add
    integer          :: i, k, Z

    real*8           :: r(3), c

    ! ---

    fmt_pdb = '(A6,I5,1X,A4,A1,A3,A1,1X,I4,A1,3X,3(F8.3),2(F6.2),6X,A4,A2,A2)'

    rec_name ='ATOM'
    atom_name = 'C'
    altloc = ''
    res_name = ''
    chainid = ''
    resseq = 1
    icode = ''

    occupancy = 1.0
    segid = ''
    element = ''
    charge = ''
    endtoken ='END'

    resseq = 0

    k = 0

    if (present(conv)) then
       write (un, '(A6,3F9.3,3F7.2,A10,I3)')  &
            "CRYST1", p%Abox(1, 1)*conv, p%Abox(2, 2)*conv, p%Abox(3, 3)*conv, 90.0d0, 90.0d0, 90.0d0, "P", 1
    else
       write (un, '(A6,3F9.3,3F7.2,A10,I3)')  &
            "CRYST1", p%Abox(1, 1), p%Abox(2, 2), p%Abox(3, 3), 90.0d0, 90.0d0, 90.0d0, "P", 1
    endif

    do i = 1, p%natloc

       add = 0

       if (p%Z(i) > 0 .and. p%Z(i) <= MAX_Z) then
          atom_name = ElementName(p%Z(i))

          if (present(conv)) then
             r = POS3(p, i)*conv
          else
             r = POS3(p, i)
          endif

          if (present(beta)) then
             c = beta(i)
          else
             c = p%g(i)
          endif

          write(un, fmt_pdb) rec_name, i, atom_name, &
               altloc, res_name, chainid, resseq, icode, &
               r, &
               occupancy, c, segid, element, charge

          k = k + 1
       else
          RAISE_ERROR("Unknown atomic number encountered.", ierror)
       endif

    enddo

    write(un,'(A3)') endtoken

  endsubroutine write_pdb_un


  !**********************************************************************
  ! Write a configuration to a PDB file
  !**********************************************************************    
  subroutine write_pdb_fn(fn, p, beta, conv, ierror)
    implicit none

    character(*), intent(in)          :: fn
    type(particles_t), intent(in)     :: p
    real(DP), intent(in), optional    :: beta(:)
    real(DP), intent(in), optional    :: conv
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: un

    ! ---

    un  = fopen(fn, F_WRITE)
    call write_pdb_un(un, p, beta, conv, ierror)
    call fclose(un)
    PASS_ERROR_WITH_INFO("Filename '" // trim(fn) // "'.", ierror)

  endsubroutine write_pdb_fn
  
endmodule pdb
