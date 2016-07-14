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
!! General tight-binding methods
!<

#include "macros.inc"
#include "filter.inc"

module dense_hamiltonian
  use supplib

  use particles
  use materials

  use dense_hamiltonian_type

  implicit none

  private

  public :: dense_hamiltonian_t

  !
  ! Interfaces
  !

  public :: init
  interface init
     module procedure dense_hamiltonian_init
  endinterface

  public :: del
  interface del
     module procedure dense_hamiltonian_del
  endinterface

  public :: update_orbitals
  interface update_orbitals
     module procedure dense_hamiltonian_update_orbitals
  endinterface

  public :: e_atomic
  interface e_atomic
     module procedure dense_hamiltonian_e_atomic
  endinterface

  public :: get_dict
  interface get_dict
     module procedure dense_hamiltonian_get_dict
  endinterface

contains

  !**********************************************************************
  ! Constructor
  !**********************************************************************
  subroutine dense_hamiltonian_init(this, mat, p, f, error)
    implicit none

    type(dense_hamiltonian_t),           intent(inout) :: this
    type(materials_t), target,           intent(in)    :: mat
    type(particles_t), target, optional, intent(in)    :: p
    integer,                   optional, intent(in)    :: f
    integer,                   optional, intent(inout) :: error

    ! ---

    INIT_ERROR(error)

    this%mat = c_loc(mat)

    this%f = 0
    if (present(f)) then
       this%f = f
    endif

    if (present(p)) then
       call update_orbitals(this, p, error)
    endif

  endsubroutine dense_hamiltonian_init


  !**********************************************************************
  ! Destructor
  !**********************************************************************
  subroutine dense_hamiltonian_del(this)
    implicit none

    type(dense_hamiltonian_t), target :: this

    ! ---

    call dense_hamiltonian_deallocate(c_loc(this))

  endsubroutine dense_hamiltonian_del


  !**********************************************************************
  ! Set the particles object
  !**********************************************************************
  subroutine dense_hamiltonian_update_orbitals(this, p, error)
    implicit none

    type(dense_hamiltonian_t), target     :: this
    type(particles_t), target, intent(in) :: p
    integer, intent(inout), optional      :: error

    ! ---

    integer   :: i, j, enr, enrj, ia
    real(DP)  :: c
#ifdef LAMMPS
    logical   :: found
#endif

    type(materials_t),    pointer :: this_mat
    type(notb_element_t), pointer :: this_at(:)

    ! ---

    INIT_ERROR(error)

    call timer_start("dense_hamiltonian_update_orbitals")

    this%p = c_loc(p)
    call c_f_pointer(this%mat, this_mat)

    this%el = c_loc(p%el(1))

    !
    ! Determine the total number of orbitals
    !

    c = 0.0
    this%norb = 0
    this%norbloc = 0
    do i = 1, p%nat

       if (IS_EL(this%f, p, i)) then

          if (.not. element_by_Z(this_mat, p%el2Z(p%el(i)), enr=enr)) then
              RAISE_ERROR_AND_STOP_TIMER("Could not find Slater-Koster tables for element '"//trim(ElementName(p%el2Z(p%el(i))))//"'.", "dense_hamiltonian_update_orbitals", error)
          endif

          this%norb = this%norb + this_mat%e(enr)%no
          if (i <= p%natloc) then
             this%norbloc = this%norbloc + this_mat%e(enr)%no
          endif

          do j = 1, i
             
             if (IS_EL(this%f, p, j)) then

                if (.not. element_by_Z(this_mat, p%el2Z(p%el(j)), enr=enrj)) then
                   RAISE_ERROR_AND_STOP_TIMER("Could not find Slater-Koster tables for element '"//trim(ElementName(p%el2Z(p%el(j))))//"'.", "dense_hamiltonian_update_orbitals", error)
                endif

                c = max(c, this_mat%cut(enr, enrj))
                c = max(c, this_mat%R(enr, enrj)%cut)
                
             endif

          enddo

       endif

    enddo

    !write (*, *)  "update_orbitals, this%norb = ", this%norb, p%nat, p%natloc

    this%cutoff = c

    call dense_hamiltonian_allocate(c_loc(this), p%nat, this%norb)
    call c_f_pointer(this%at, this_at, [this%nat])

    ia   = 1
    do i = 1, p%nat

       if (IS_EL(this%f, p, i)) then

#ifdef LAMMPS
          if (i <= p%natloc) then
#endif

             if (.not. element_by_Z(this_mat, p%el2Z(p%el(i)), enr=enr)) then
                RAISE_ERROR_AND_STOP_TIMER("Could not find Slater-Koster tables for element '"//trim(ElementName(p%el2Z(p%el(i))))//"'.", "dense_hamiltonian_update_orbitals", error)
             endif
             this_at(i)     = this_mat%e(enr)
             this_at(i)%o1  = ia
             ia             = ia + this_at(i)%no

!             write (*, *) i//" is "//this_at(i)%enr

#ifdef LAMMPS
          else
             ! FIXME!!! Slow, but should not be the time relevant step in TB
             found = .false.
             do j = 1, p%natloc
                if (p%tag(i) == p%tag(j)) then
                   this_at(i) = this_at(j)
!                   write (*, *) i//"->"//j//"; "//p%tag(i)//"->"//p%tag(j)//"; is "//this_at(i)%enr
                   found = .true.
                endif
             enddo

             if (.not. found) then
                RAISE_ERROR_AND_STOP_TIMER("Could not find tag "//p%tag(i)//" of atom "//i//" in simulation.", "dense_hamiltonian_update_orbitals", error)
             endif
          endif
#endif

       endif

    enddo

    call timer_stop("dense_hamiltonian_update_orbitals")

  endsubroutine dense_hamiltonian_update_orbitals


  !>
  !! Atomic energies
  !!
  !! Return the energy of the system decomposed into charge neutral
  !! isolated atoms.
  !!
  !<
  real(DP) function dense_hamiltonian_e_atomic(this, p, error) result(e)
    implicit none

    type(dense_hamiltonian_t), intent(in)     :: this   !< Hamiltonian object
    type(particles_t),         intent(in)     :: p
    integer,         optional, intent(inout)  :: error  !< Error signals

    ! ---

    integer  :: i, a0, q0

    type(notb_element_t), pointer  :: at(:)

    ! ---

    INIT_ERROR(error)
    
    call c_f_pointer(this%at, at, [this%nat])

    e = 0.0_DP
    do i = 1, p%natloc
       q0 = int(at(i)%q0)

       do a0 = 1, q0/2
          e = e + 2*at(i)%e(get_orbital(at(i)%no, a0))
       enddo
       if (mod(q0, 2) /= 0) then
          e = e + at(i)%e(get_orbital(at(i)%no, q0/2+1))
       endif
       ! If atom has fractional charge...
       e = e + (at(i)%q0-q0)*at(i)%e(get_orbital(at(i)%no, q0/2+1))
    enddo

  endfunction dense_hamiltonian_e_atomic


  !>
  !! Return dictionary object containing pointers to internal data
  !<
  subroutine dense_hamiltonian_get_dict(this, dict, error)
    implicit none

    type(dense_hamiltonian_t), intent(inout) :: this        !< NOTB object
    type(ptrdict_t),           intent(inout) :: dict
    integer,         optional, intent(out)   :: error       !< Error signals

    ! ---

    INIT_ERROR(error)

    if (this%nk == 1) then
       if (c_associated(this%H)) then
          call ptrdict_register_array2d_property(dict%ptrdict, this%H, &
                                                 this%norb, this%norb, &
                                                 CSTR("Hamiltonian_matrix"), &
                                                 CSTR("N/A"))
       endif 
       if (c_associated(this%S)) then
          call ptrdict_register_array2d_property(dict%ptrdict, this%S, &
                                                 this%norb, this%norb, &
                                                 CSTR("overlap_matrix"), &
                                                 CSTR("N/A"))
       endif
       if (c_associated(this%rho)) then
          call ptrdict_register_array2d_property(dict%ptrdict, this%rho, &
                                                 this%norb, this%norb, &
                                                 CSTR("density_matrix"), &
                                                 CSTR("N/A"))
       endif
    else
       if (c_associated(this%H)) then
          call ptrdict_register_array3d_property(dict%ptrdict, this%H, &
                                                 this%norb, this%norb, this%nk, &
                                                 CSTR("Hamiltonian_matrix"), &
                                                 CSTR("N/A"))
       endif
       if (c_associated(this%S)) then
          call ptrdict_register_array3d_property(dict%ptrdict, this%S, &
                                                 this%norb, this%norb, this%nk, &
                                                 CSTR("overlap_matrix"), &
                                                 CSTR("N/A"))
       endif
       if (c_associated(this%rho)) then
          call ptrdict_register_array3d_property(dict%ptrdict, this%rho, &
                                                 this%norb, this%norb, this%nk, &
                                                 CSTR("density_matrix"), &
                                                 CSTR("N/A"))
       endif
    endif

  endsubroutine dense_hamiltonian_get_dict

endmodule dense_hamiltonian
