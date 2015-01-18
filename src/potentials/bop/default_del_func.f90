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
  !! Remove the object
  !<
  subroutine DEL_FUNC(this)
    implicit none

    type(BOP_TYPE), intent(inout) :: this

    ! ---

    if (this%neighbor_list_allocated) then
       deallocate(this%neb)
       deallocate(this%nbb)
#ifndef LAMMPS
       deallocate(this%dcell)
#endif
       deallocate(this%bndtyp)
       deallocate(this%bndlen)
       deallocate(this%bndnm)
       deallocate(this%cutfcnar)
       deallocate(this%cutdrvar)

#ifdef SCREENING
       deallocate(this%cutfcnbo)
       deallocate(this%cutdrvbo)
       deallocate(this%sneb_seed)
       deallocate(this%sneb_last)
       deallocate(this%sneb)
       deallocate(this%sbnd)
       deallocate(this%sfacbo)
       deallocate(this%cutdrarik)
       deallocate(this%cutdrarjk)
       deallocate(this%cutdrboik)
       deallocate(this%cutdrbojk)
#endif

       this%neighbor_list_allocated = .false.
    endif

  endsubroutine DEL_FUNC
