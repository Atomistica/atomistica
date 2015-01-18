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

module input_trajectory
  use libAtoms_module

  use particles

  use native_io
  use cfg
  use xyz

#ifdef HAVE_NETCDF
  use nc 
#endif

  implicit none

  integer, parameter  :: IN_ATOMS  = 0
  integer, parameter  :: IN_CFG    = 1
  integer, parameter  :: IN_NC     = 2
  integer, parameter  :: IN_XYZ    = 3


  integer, parameter  :: EXTS_LEN  = 5

  character(4), parameter :: EXTS(EXTS_LEN) = &
       (/ ".dat", ".out", ".cfg", ".nc ", ".xyz" /)
  integer, parameter :: KINDS(EXTS_LEN) = &
       (/ IN_ATOMS, IN_ATOMS, IN_CFG, IN_NC, IN_XYZ /)

  type input_t

     integer        :: nframes

     integer        :: kind

     character(80)  :: fn

     !
     ! Modules
     !

#ifdef HAVE_NETCDF
     type(nc_t)     :: nc
#endif

  endtype input_t


  interface open
     module procedure input_open
  endinterface

  interface close
     module procedure input_close
  endinterface

  interface get_frame
     module procedure input_get_frame
  endinterface

contains

  !****************************************************************
  ! Open a file
  !****************************************************************
  subroutine input_open(this, fn, p, error)
    implicit none

    type(input_t), intent(inout)      :: this
    character*(*)                     :: fn
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: error

    ! ---

    integer        :: i, j
    character(80)  :: ext

    ! ---

    this%kind = -1

    i = 0
    j = index(fn, ".")
    do while (j > 0)
       i = i + j
       j = index(fn(i+1:), ".")
    enddo

    ext = ""
    if (i > 0) then

       ext = fn(i:len(fn))

       do i = 1, EXTS_LEN
          if (trim(ext) == trim(EXTS(i))) then
             this%kind = KINDS(i)
          endif
       enddo

    endif

    selectcase(this%kind)
    case(IN_ATOMS)
       call read_atoms(p, fn)
!       call cyclic_from_cyc_dat(p, "cyc.dat")
       this%nframes = 1
    case(IN_CFG)
       call read_cfg(p, fn, error=error)
       PASS_ERROR(error)
       this%nframes = 1
#ifdef HAVE_NETCDF
    case(IN_NC)
       call open(this%nc, p, fn)
       this%nframes = this%nc%nframes
#endif
    case(IN_XYZ)
       ! Only single frames are supported for now
       call read_xyz(p, fn, error=error)
       PASS_ERROR(error)
       this%nframes = 1
    case default
       RAISE_ERROR("Unknown file extension encountered: '" // trim(ext) // "'.", error)
    endselect

  endsubroutine input_open


  !****************************************************************
  ! Read a frame from a file
  !****************************************************************
  subroutine input_get_frame(this, it, ti, p)
    implicit none

    type(input_t), intent(in)         :: this
    integer, intent(in)               :: it
    real(DP), intent(out)             :: ti
    type(particles_t), intent(inout)  :: p

    ! ---

    ti  = 0.0_DP

#ifdef HAVE_NETCDF
    if (this%kind == IN_NC) then
       call read_frame(this%nc, it, ti, p)
    endif
#endif

  endsubroutine input_get_frame


  !****************************************************************
  ! Close a NetCDF file
  !****************************************************************
  subroutine input_close(this)
    implicit none

    type(input_t), intent(inout)  :: this

    ! ---

#ifdef HAVE_NETCDF
    if (this%kind == IN_NC) then
       call close(this%nc)
    endif
#endif

  endsubroutine input_close

endmodule input_trajectory
