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
!! Write LAMMPS native atom data
!!
!! Write LAMMPS native atom data
!<
module lammps_data
  use libAtoms_module

  use io

  use particles
  use molecules

  interface write_lammps_data
     module procedure write_lammps_data_un, write_lammps_data_fn
  endinterface

contains

  !>
  !! Write LAMMPS atom data file with unit number
  !!
  !! Write LAMMPS atom data file with unit number. At the moment,
  !! a hacked version that works for me, now. (Tommi Jaervi, 18.11.2009)
  !! (XXX marks some hardcoded features.)
  !!
  !! Features: Rectangular box assumed.
  !!
  !! Mode       Corresponding lammps atom_style (see lammps documentation for command read_data)
  !! ---------------------------
  !! 0          atomic
  !! 1          charge
  !! 2          full (+ identify water molecules to be output as molecules, output bonds and angles)
  !<
  subroutine write_lammps_data_un(un, p, ntypes, types, outtypes, mode, mol, q, ierror)
    implicit none

    integer, intent(in)                    :: un        !< File unit
    type(particles_t), intent(in)          :: p         !< Particles object
    integer, intent(in)                    :: ntypes    !< Number of types
    integer, dimension(:), intent(in)      :: types     !< Atom types, atomic number Z (mdcore)
    integer, dimension(:), intent(in)      :: outtypes  !< Atom type number (LAMMPS)
    integer, intent(in)                    :: mode      !< Mode
    type(molecules_t), intent(inout), optional   :: mol       !< Molecules object
    real(DP), dimension(:), intent(in), optional :: q         !< Atom charges
    integer, intent(out), optional         :: ierror    !< Error signals

    ! ---

    integer  :: i, j, k, l                              ! Loops etc.
    integer, dimension(100) :: z_to_outtype             ! Conversion from Z to LAMMPS type

    integer  :: nangles, nbonds                         ! Number of angles and bonds (mode 2)
    integer  :: natmol                                  ! Number of atoms in this molecule (mode 2)
    integer, dimension(3) :: atmol                      ! Atoms in this molecule (mode 2)
    integer  :: nh, no                                  ! To identify water (mode 2)

    ! --- Checks

    INIT_ERROR(ierror)

    if(mode < 0 .or. mode > 2) then
       RAISE_ERROR("write_lammps_data_un: Unknown mode " // mode // ".", ierror)
    end if
    if (mode==1) then
       if (.not. present(q)) then
          RAISE_ERROR("write_lammps_data_un: Mode 1 requires charge array.", ierror)
       endif
    endif
    if(mode==2) then
       if (present(mol)) then
          if (.not. mol%use_imol) then
             RAISE_ERROR("write_lammps_data_un: Mode 2 requires imol-array in molecules!", ierror)
          endif
       else
          RAISE_ERROR("write_lammps_data_un: Mode 2 requires molecules object!", ierror)
       endif
    end if

    ! --- Construct stuff needed below

    ! Molecule id's
    if(mode==2) then
       call molecules_update_head(mol, p, ierror)
       PASS_ERROR(ierror)
    end if

    ! Z -> LAMMPS type conversion
    z_to_outtype = 0
    do i = 1, ntypes
       z_to_outtype(types(i)) = outtypes(i)
    end do

    ! Mode 2: Count the number of H2O molecules
    if(mode==2) then
       nbonds = 0  ! number of H2O's
       ! loop over molecules
       do i = 1, mol%n_molec
          ! collect atoms
          natmol = 0
          no = 0
          nh = 0
          j = mol%head(i)
          do while(j > 0)
             natmol = natmol + 1
             atmol(natmol) = j

             ! mark O and count O, H
             if(p%Z(j)==8) then
                no = no + 1
             end if
             if(p%Z(j)==1) then
                nh = nh + 1
             end if

             j = mol%next(j)
          end do

          ! atoms collected, check if it's H2O
          if(natmol==3 .and. no==1 .and. nh==2) then
             nbonds = nbonds + 1
          end if
       end do  ! end of loop over molecules
    end if


    ! --- Start

    write (un,*), "# LAMMPS data written by MDCORE, write_lammps_data_un()"
    write (un,*), ""
    write (un,*), p%nat, "atoms"
    write (un,*), ntypes, "atom types"
    write (un,*), ""
    write (un,*), 0.0, p%Abox(1,1), "xlo xhi"
    write (un,*), 0.0, p%Abox(2,2), "ylo yhi"
    write (un,*), 0.0, p%Abox(3,3), "zlo zhi"
    write (un,*), ""
    if(mode==2) then
       write (un,*), "1 bond types" ! XXX
       write (un,*), "1 angle types" ! XXX
       write (un,*), nbonds*2, "bonds"  ! XXX
       write (un,*), nbonds, "angles"  ! XXX
       write (un,*), ""
    end if
    write (un,*), "Masses"
    write (un,*), ""
    do i = 1, ntypes
       if (types(i) > 0 .and. types(i) <= MAX_Z) then
          write (un,*), outtypes(i), ElementMass(types(i))
       else
          RAISE_ERROR("Unknown atomic number " // types(i) // ".", ierror)
       endif
    end do
    write (un,*), ""
    write (un,*), "Atoms"
    write (un,*), ""
    do i = 1, p%nat
       if (mode == 0) then
          write (un,'(2I10,3ES20.10)'), i, z_to_outtype(p%Z(i)), POS3(p, i)
       else if(mode==1) then
          write (un,'(2I10,4ES20.10)'), i, z_to_outtype(p%Z(i)), q(i), POS3(p, i)
       else if(mode==2) then
          write (un,'(3I10,4ES20.10)'), i, mol%imol(i), z_to_outtype(p%Z(i)), q(i), POS3(p, i)
       end if
    end do

    ! --- XXX: H2O specific bonds
    if(mode==2) then
       write (un,*), ""
       write (un,*), "Bonds"
       write (un,*), ""
       nbonds = 0
       ! loop over molecules
       do i = 1, mol%n_molec
          ! collect atoms
          natmol = 0
          no = 0
          nh = 0
          j = mol%head(i)
          do while(j > 0)
             natmol = natmol + 1
             atmol(natmol) = j

             ! mark O and count O, H
             if(p%Z(j)==8) then
                k = natmol
                no = no + 1
             end if
             if(p%Z(j)==1) then
                nh = nh + 1
             end if

             j = mol%next(j)
          end do

          ! atoms collected, write H2O data
          if(natmol==3 .and. no==1 .and. nh==2) then
             ! bonds
             do l = 1, natmol
                if(l /= k) then
                   ! bonds from H to O, that's why O omitted
                   nbonds = nbonds + 1
                   write (un,*), nbonds, 1, atmol(k), atmol(l)  ! XXX: bond type 1
                end if
             end do
          end if
       end do  ! end of loop over molecules
    end if

    ! --- XXX: H2O specific angles
    if(mode==2) then
       write (un,*), ""
       write (un,*), "Angles"
       write (un,*), ""
       nangles = 0
       ! loop over molecules
       do i = 1, mol%n_molec
          ! collect atoms
          natmol = 0
          no = 0
          nh = 0
          j = mol%head(i)
          do while(j > 0)
             natmol = natmol + 1
             atmol(natmol) = j

             ! mark O and count O, H
             if(p%Z(j)==8) then
                k = natmol
                no = no + 1
             end if
             if(p%Z(j)==1) then
                nh = nh + 1
             end if

             j = mol%next(j)
          end do

          ! atoms collected, write H2O data
          if(natmol==3 .and. no==1 .and. nh==2) then
             ! angles
             nangles = nangles + 1
             if(k==1) then
                write (un,*), nangles, 1, atmol(2), atmol(1), atmol(3)  ! XXX: bond type 1
             else if(k==2) then
                write (un,*), nangles, 1, atmol(1), atmol(2), atmol(3)  ! XXX: bond type 1
             else
                write (un,*), nangles, 1, atmol(1), atmol(3), atmol(2)  ! XXX: bond type 1
             end if
          end if
       end do  ! end of loop over molecules
    end if

  end subroutine write_lammps_data_un


  !>
  !! Write LAMMPS atom data file with file name
  !!
  !! Write LAMMPS atom data file with file name
  !<
  subroutine write_lammps_data_fn(fn, p, ntypes, types, outtypes, mode, mol, q, ierror)
    implicit none

    character(*), intent(in)               :: fn        !< File name
    type(particles_t), intent(in)          :: p         !< Particles object
    integer, intent(in)                    :: ntypes    !< Number of types
    integer, dimension(:), intent(in)      :: types     !< Atom types, atomic number Z (mdcore)
    integer, dimension(:), intent(in)      :: outtypes  !< Atom type number (LAMMPS)
    integer, intent(in)                    :: mode      !< Mode
    type(molecules_t), intent(inout), optional       :: mol       !< Molecules object
    real(DP), dimension(:), intent(in), optional     :: q         !< Atom charges
    integer, intent(out), optional         :: ierror    !< Error signals

    ! ---

    integer  :: un

    ! ---

    INIT_ERROR(ierror)

    un  = fopen(fn, F_WRITE)
    call write_lammps_data_un(un, p, ntypes, types, outtypes, mode, mol, q, ierror)
    PASS_ERROR(ierror)
    call fclose(un)

  end subroutine write_lammps_data_fn

end module lammps_data
