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
!! Read, write and modify NetCDF trajectory files
!!
!! This module reads, writes and modifies NetCDF trajectory files. The file format
!! follows the AMBER conventions which VMD can read. Additional fields can also be
!! added which can be displayed with AtomEye.
!!
!! See
!!
!! http://ambermd.org/formats.html
!!
!! http://www.ks.uiuc.edu/Research/vmd/
!!
!! http://www.csanyi.net and go to "group" -> "group software" -> "AtomEye"
!!
!! This is included for compile option -DHAVE_NETCDF. Otherwise a
!! dummy module with error messages is used. Don't forget to add a
!! corresponding subroutine to the dummy module when adding
!! functionality.
!<

#ifdef HAVE_NETCDF

module nc
  use supplib

  use particles

#ifdef _MP

  use mpi
  use pnetcdf
  use communicator

#define nf90_close nf90mpi_close
#define nf90_def_dim nf90mpi_def_dim
#define nf90_def_var nf90mpi_def_var
#define nf90_enddef nf90mpi_enddef
#define nf90_get_var nf90mpi_get_var
#define nf90_inq_dimid nf90mpi_inq_dimid
#define nf90_inq_varid nf90mpi_inq_varid
#define nf90_inquire_dimension nf90mpi_inquire_dimension
#define nf90_inquire_variable nf90mpi_inquire_variable
#define nf90_put_att nf90mpi_put_att
#define nf90_put_var nf90mpi_put_var
#define nf90_redef nf90mpi_redef
#define nf90_strerror nf90mpi_strerror
#define nf90_sync nf90mpi_sync

#else

  use netcdf

#endif

  use versioninfo

#ifndef _MP
  use iso_fortran_env
#endif

  implicit none

#ifndef _MP
  integer, parameter :: MPI_OFFSET_KIND = ATOMIC_INT_KIND
#endif

  character(*), parameter, private  :: MODULE_STR = "NC"

  character(*), parameter, private  :: NC_FRAME_STR         = "frame"
  character(*), parameter, private  :: NC_SPATIAL_STR       = "spatial"
  character(*), parameter, private  :: NC_ATOM_STR          = "atom"
  character(*), parameter, private  :: NC_CELL_SPATIAL_STR  = "cell_spatial"
  character(*), parameter, private  :: NC_CELL_ANGULAR_STR  = "cell_angular"
  character(*), parameter, private  :: NC_LABEL_STR         = "label"

  character(*), parameter, private  :: NC_TIME_STR          = "time"
  character(*), parameter, private  :: NC_CELL_ORIGIN_STR   = "cell_origin"
  character(*), parameter, private  :: NC_CELL_LENGTHS_STR  = "cell_lengths"
  character(*), parameter, private  :: NC_CELL_ANGLES_STR   = "cell_angles"

  character(*), parameter, private  :: NC_SHEAR_DX_STR      = "shear_dx"

  character(*), parameter, private  :: NC_UNITS_STR         = "units"
  character(*), parameter, private  :: NC_SCALE_FACTOR_STR  = "scale_factor"

  type nc_t

     !
     ! Mode (read/write) and NetCDF file handle
     !

     integer               :: mode
     integer               :: ncid

     !
     ! Total number of frames in file and current
     ! frame for consecutive writes
     !

     integer               :: nframes
     integer               :: frame_no

     !
     ! Amber convention
     !

     integer               :: frame_dim
     integer               :: spatial_dim
     integer               :: atom_dim
     integer               :: cell_spatial_dim
     integer               :: cell_angular_dim
     integer               :: label_dim

     integer               :: spatial_var
     integer               :: cell_spatial_var
     integer               :: cell_angular_var

     integer               :: time_var
     integer               :: cell_origin_var
     integer               :: cell_lengths_var
     integer               :: cell_angles_var

     !
     ! MDCore convention
     !

     integer               :: Z_var
     integer               :: shear_dx_var

     !
     ! Dynamic fields
     !

     integer, pointer      :: real_attr_var(:)
     integer, pointer      :: real_attr_ndims(:)
     integer, pointer      :: integer_attr_var(:)
     integer, pointer      :: integer_attr_ndims(:)
     integer, pointer      :: real3_attr_var(:)
     integer, pointer      :: real3_attr_ndims(:)
     integer, pointer      :: real3x3_attr_var(:)
     integer, pointer      :: real3x3_attr_ndims(:)

     integer, pointer      :: real_var(:)
     integer, pointer      :: real_ndims(:)
     integer, pointer      :: integer_var(:)
     integer, pointer      :: integer_ndims(:)
     integer, pointer      :: real3_var(:)
     integer, pointer      :: real3_ndims(:)
     integer, pointer      :: real3x3_var(:)
     integer, pointer      :: real3x3_ndims(:)

     !
     ! Temporary buffers
     !

#ifndef _MP
     real(DP), pointer      :: tmp_real(:)
     integer, pointer       :: tmp_integer(:)
     real(DP), pointer      :: tmp_real3(:, :)
     real(DP), pointer      :: tmp_real3x3(:, :, :)
#endif

  endtype nc_t

  interface create
     module procedure nc_create
  endinterface

  interface open
     module procedure nc_open
  endinterface

  interface close
     module procedure nc_close
  endinterface

  interface get_time
     module procedure nc_get_time
  endinterface

  interface find_frame
     module procedure nc_find_frame
  endinterface

  interface read_frame
     module procedure nc_read_frame
  endinterface

  interface write_constant
     module procedure nc_write_constant
  endinterface

  interface write_frame
     module procedure nc_write_frame
  endinterface

  interface write_field
     module procedure nc_write_field
  endinterface

contains

#define CHECK_NETCDF_ERROR(x, ierror)  if (x /= NF90_NOERR) then ; RAISE_ERROR("NetCDF error: " // trim(nf90_strerror(x)), ierror) ; endif
#define CHECK_NETCDF_ERROR_WITH_INFO(x, info, ierror)  if (x /= NF90_NOERR) then ; RAISE_ERROR(info // trim(nf90_strerror(x)), ierror) ; endif

  !>
  !! Write the prmtop (topology) file
  !!
  !! Write the prmtop (topology) file. Required for VMD only.
  !<
  subroutine write_prmtop(p, fn, ierror)
    implicit none

    type(particles_t), intent(in)  :: p
    character(*), intent(in)       :: fn
    integer, intent(inout), optional :: ierror

    ! ---

    integer  :: un, i

    ! ---

    un = fopen(fn, F_WRITE)

    write (un, '(A)')          "%VERSION  MDCore"
    write (un, '(A)')          "%FLAG TITLE"
    write (un, '(A)')          "%FORMAT(20a4)"
    write (un, '(A4,75X,A1)')  "NASN", " "

    write (un, '(A)')          "%FLAG POINTERS"
    write (un, '(A)')          "%FORMAT(10I8)"
    write (un, '(10I8)')       p%nat, (0, i = 1, 11)
    write (un, '(10I8)')       (0, i = 1, 12)
    write (un, '(6I8)')        (0, i = 1, 6)

    write (un, '(A)')          "%FLAG ATOM_NAME"
    write (un, '(A)')          "%FORMAT(20a4)"
    write (un, '(20A4)')       (p%sym(p%global2local(i)), i = 1, p%nat)

    write (un, '(A)')          "%FLAG CHARGE"
    write (un, '(A)')          "%FORMAT(5E16.5)"
    write (un, '(5E16.5)')     (0.0_DP, i = 1, p%nat)

    write (un, '(A)')          "%FLAG MASS"
    write (un, '(A)')          "%FORMAT(5E16.5)"
    write (un, '(5E16.5)')     (p%m(p%global2local(i)), i = 1, p%nat)

    !    write (un, '(12I6)')    (0, i = 1, p%nat)

    !    write (un, '(12I6)')    (0, i = 1, p%nat)

    call fclose(un)

  endsubroutine write_prmtop


  !>
  !! Create a new NetCDF trajectory file
  !!
  !! Create a new NetCDF trajectory file
  !<
  subroutine nc_create(this, p, fn, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    character(*), intent(in)          :: fn
    integer, intent(inout), optional  :: ierror

    ! ---

    character(1000) :: versionstr

    integer                       :: i
    integer(kind=MPI_OFFSET_KIND) :: totnat

    ! ---

    this%mode            = F_WRITE
    this%nframes         = 0
    this%cell_origin_var = -1

#ifdef _MP
    CHECK_NETCDF_ERROR( nf90mpi_create(mod_communicator%mpi%communicator, fn, NF90_64BIT_OFFSET, MPI_INFO_NULL, this%ncid), ierror )
#else
    CHECK_NETCDF_ERROR( nf90_create(fn, NF90_64BIT_OFFSET, this%ncid), ierror )
#endif

    !
    ! Dimensions
    !

    CHECK_NETCDF_ERROR( nf90_def_dim(this%ncid, NC_FRAME_STR, nf90_unlimited, this%frame_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_def_dim(this%ncid, NC_SPATIAL_STR, 3, this%spatial_dim), ierror )
    totnat = p%totnat
    CHECK_NETCDF_ERROR( nf90_def_dim(this%ncid, NC_ATOM_STR, totnat, this%atom_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_def_dim(this%ncid, NC_CELL_SPATIAL_STR, 3, this%cell_spatial_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_def_dim(this%ncid, NC_CELL_ANGULAR_STR, 3, this%cell_angular_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_def_dim(this%ncid, NC_LABEL_STR, 10, this%label_dim), ierror )

    !
    ! Variables
    !

    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_SPATIAL_STR, NF90_CHAR, (/ this%spatial_dim /), this%spatial_var), ierror )
    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_CELL_SPATIAL_STR, NF90_CHAR, (/ this%spatial_dim /), this%cell_spatial_var), ierror )
    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_CELL_ANGULAR_STR, NF90_CHAR, (/ this%label_dim, this%spatial_dim /), this%cell_angular_var), ierror )

    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_TIME_STR, NF90_FLOAT, (/ this%frame_dim /), this%time_var), ierror )
    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_CELL_LENGTHS_STR, NF90_FLOAT, (/ this%cell_spatial_dim, this%frame_dim /), this%cell_lengths_var), ierror )
    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_CELL_ANGLES_STR, NF90_FLOAT, (/ this%cell_angular_dim, this%frame_dim /), this%cell_angles_var), ierror )

!    CHECK_NETCDF_ERROR( nf90_def_var(this%ncid, NC_SHEAR_DX_STR, NF90_FLOAT, (/ this%spatial_dim, this%frame_dim /), this%shear_dx_var), ierror )
    ! Shear is read, but not written. (True cell is written. Reading support
    ! is for backwards compatibility.)
    this%shear_dx_var = -1

    !
    ! Dynamic variables
    !

    this%Z_var  = -1

    !
    ! Attributes
    !

    this%real_attr_var    => NULL()
    this%real_attr_ndims  => NULL()

    if (p%data%n_real_attr > 0) then
       allocate(this%real_attr_var(p%data%n_real_attr))
       allocate(this%real_attr_ndims(p%data%n_real_attr))

       this%real_attr_var    = -1
       this%real_attr_ndims  = -1

       do i = 1, p%data%n_real_attr
          if (iand(p%data%tag_real_attr(i), F_TO_TRAJ) /= 0) then
             CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real_attr(i), NF90_FLOAT, (/ this%frame_dim /), this%real_attr_var(i) ), ierror )
             this%real_attr_ndims(i)  = 1
          endif
       enddo
    endif

    this%integer_attr_var    => NULL()
    this%integer_attr_ndims  => NULL()

    if (p%data%n_integer_attr > 0) then
       allocate(this%integer_attr_var(p%data%n_integer_attr))
       allocate(this%integer_attr_ndims(p%data%n_integer_attr))

       this%integer_attr_var    = -1
       this%integer_attr_ndims  = -1

       do i = 1, p%data%n_integer_attr
          if (iand(p%data%tag_integer_attr(i), F_TO_TRAJ) /= 0) then
             CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_integer_attr(i), NF90_INT, (/ this%frame_dim /), this%integer_attr_var(i) ), ierror )
             this%integer_attr_ndims(i)  = 1
          endif
       enddo
    endif

    this%real3_attr_var    => NULL()
    this%real3_attr_ndims  => NULL()

    if (p%data%n_real3_attr > 0) then
       allocate(this%real3_attr_var(p%data%n_real3_attr))
       allocate(this%real3_attr_ndims(p%data%n_real3_attr))

       this%real3_attr_var    = -1
       this%real3_attr_ndims  = -1

       do i = 1, p%data%n_real3_attr
          if (iand(p%data%tag_real3_attr(i), F_TO_TRAJ) /= 0) then
             CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3_attr(i), NF90_FLOAT, (/ this%spatial_dim, this%frame_dim /), this%real3_attr_var(i) ), ierror )
             this%real3_attr_ndims(i)  = 2
          endif
       enddo
    endif

    this%real3x3_attr_var    => NULL()
    this%real3x3_attr_ndims  => NULL()

    if (p%data%n_real3x3_attr > 0) then
       allocate(this%real3x3_attr_var(p%data%n_real3x3_attr))
       allocate(this%real3x3_attr_ndims(p%data%n_real3x3_attr))

       this%real3x3_attr_var    = -1
       this%real3x3_attr_ndims  = -1

       do i = 1, p%data%n_real3x3_attr
          if (iand(p%data%tag_real3x3_attr(i), F_TO_TRAJ) /= 0) then
             CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3x3_attr(i), NF90_FLOAT, (/ this%spatial_dim, this%spatial_dim, this%frame_dim /), this%real3x3_attr_var(i) ), ierror )
             this%real3x3_attr_ndims(i)  = 3
          endif
       enddo
    endif

    !
    ! Fields
    !

    this%real_var    => NULL()
    this%real_ndims  => NULL()

    if (p%data%n_real > 0) then
       allocate(this%real_var(p%data%n_real))
       allocate(this%real_ndims(p%data%n_real))

       this%real_var    = -1
       this%real_ndims  = -1

       do i = 1, p%data%n_real
          if (iand(p%data%tag_real(i), F_TO_TRAJ) /= 0) then
             if (iand(p%data%tag_real(i), F_CONSTANT) /= 0) then
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real(i), NF90_FLOAT, (/ this%atom_dim /), this%real_var(i) ), ierror )
                this%real_ndims(i)  = 1
             else
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real(i), NF90_FLOAT, (/ this%atom_dim, this%frame_dim /), this%real_var(i) ), ierror )
                this%real_ndims(i)  = 2
             endif

             CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real_var(i), NC_UNITS_STR, p%data%unit_real(i)), ierror )
             CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real_var(i), NC_SCALE_FACTOR_STR, p%data%conv_real(i)), ierror )
          endif
       enddo
    endif

    this%integer_var    => NULL()
    this%integer_ndims  => NULL()

    if (p%data%n_integer > 0) then
       allocate(this%integer_var(p%data%n_integer))
       allocate(this%integer_ndims(p%data%n_integer))

       this%integer_var    = -1
       this%integer_ndims  = -1

       do i = 1, p%data%n_integer
          if (iand(p%data%tag_integer(i), F_TO_TRAJ) /= 0) then
             if (iand(p%data%tag_integer(i), F_CONSTANT) /= 0) then
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_integer(i), NF90_INT, (/ this%atom_dim /), this%integer_var(i) ), ierror )

                if (trim(p%data%name_integer(i)) == trim(Z_STR)) then
                   this%Z_var  = this%integer_var(i)
                endif

                this%integer_ndims(i)  = 1
             else
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_integer(i), NF90_INT, (/ this%atom_dim, this%frame_dim /), this%integer_var(i) ), ierror )
                this%integer_ndims(i)  = 2
             endif
          endif
       enddo
    endif

    this%real3_var    => NULL()
    this%real3_ndims  => NULL()

    if (p%data%n_real3 > 0) then
       allocate(this%real3_var(p%data%n_real3))
       allocate(this%real3_ndims(p%data%n_real3))

       this%real3_var    = -1
       this%real3_ndims  = -1

       do i = 1, p%data%n_real3
          if (iand(p%data%tag_real3(i), F_TO_TRAJ) /= 0) then
             if (iand(p%data%tag_real3(i), F_CONSTANT) /= 0) then
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3(i), NF90_FLOAT, (/ this%spatial_dim, this%atom_dim /), this%real3_var(i) ), ierror )
                this%real3_ndims(i)  = 2
             else
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3(i), NF90_FLOAT, (/ this%spatial_dim, this%atom_dim, this%frame_dim /), this%real3_var(i) ), ierror )
                this%real3_ndims(i)  = 3
             endif

             CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real3_var(i), NC_UNITS_STR, p%data%unit_real3(i)), ierror )
             CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real3_var(i), NC_SCALE_FACTOR_STR, p%data%conv_real3(i)), ierror )
          endif
       enddo
    endif

    this%real3x3_var    => NULL()
    this%real3x3_ndims  => NULL()

    if (p%data%n_real3x3 > 0) then
       allocate(this%real3x3_var(p%data%n_real3x3))
       allocate(this%real3x3_ndims(p%data%n_real3x3))

       this%real3x3_var    = -1
       this%real3x3_ndims  = -1

       do i = 1, p%data%n_real3x3
          if (iand(p%data%tag_real3x3(i), F_TO_TRAJ) /= 0) then
             if (iand(p%data%tag_real3x3(i), F_CONSTANT) /= 0) then
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3x3(i), NF90_FLOAT, (/ this%spatial_dim, this%spatial_dim, this%atom_dim /), this%real3x3_var(i) ), ierror )
                this%real3x3_ndims(i)  = 2
             else
                CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3x3(i), NF90_FLOAT, (/ this%spatial_dim, this%spatial_dim, this%atom_dim, this%frame_dim /), this%real3x3_var(i) ), ierror )
                this%real3x3_ndims(i)  = 3
             endif

             CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real3x3_var(i), NC_UNITS_STR, p%data%unit_real3x3(i)), ierror )
             CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real3x3_var(i), NC_SCALE_FACTOR_STR, p%data%conv_real3x3(i)), ierror )
          endif
       enddo
    endif

    !
    ! Attributes
    !

    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, NF90_GLOBAL, "Conventions", "AMBER"), ierror )
    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, NF90_GLOBAL, "ConventionVersion", "1.0"), ierror )
    versionstr = &
         "Atomistica revision: " // trim(atomistica_revision) // &
         ", build date: " // trim(builddate) // &
         ", build host: " // trim(buildhost)
    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, NF90_GLOBAL, "program", "MDCore"), ierror )
    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, NF90_GLOBAL, "programVersion", trim(versionstr)), ierror )

    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%time_var, NC_UNITS_STR, "picosecond"), ierror )
    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%cell_lengths_var, NC_UNITS_STR, "angstrom"), ierror )
    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%cell_angles_var, NC_UNITS_STR, "degree"), ierror )

!    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%shear_dx_var, NC_UNITS_STR, "angstrom"), ierror )

    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%time_var, NC_SCALE_FACTOR_STR, time_to_fs/1000), ierror )
    CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%cell_lengths_var, NC_SCALE_FACTOR_STR, length_to_A), ierror )

    !
    ! Finished with definition
    !

    CHECK_NETCDF_ERROR( nf90_enddef(this%ncid), ierror )

#ifdef _MP
    CHECK_NETCDF_ERROR( nf90mpi_begin_indep_data(this%ncid), ierror )
    if (mpi_id() == ROOT) then
#endif

    !
    ! Write label variables
    !

    CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%spatial_var, "xyz"), ierror )
    CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%cell_spatial_var, "abc"), ierror )
    CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%cell_angular_var, (/ "alpha", "beta ", "gamma" /) ), ierror )

#ifdef _MP
    endif
    CHECK_NETCDF_ERROR( nf90mpi_end_indep_data(this%ncid), ierror )
#endif

    !
    ! Write initial configuration
    !

    this%frame_no = 1

    !
    ! Allocate buffers
    !

#ifndef _MP
    allocate(this%tmp_real(p%nat))
    allocate(this%tmp_integer(p%nat))
    allocate(this%tmp_real3(3, p%nat))
    allocate(this%tmp_real3x3(3, 3, p%nat))
#endif

    !
    ! Write constant information
    !

    call nc_write_constants(this, p)

  endsubroutine nc_create


  !>
  !! Open a NetCDF file
  !!
  !! Open the NetCDF file \param fn. If \mode is F_WRITE the file will be opened for
  !! write operations. Additionally specifying \param add_missing to true modifies the
  !! data structure to match the one given in \param p.
  !<
  subroutine nc_open(this, p, fn, mode, add_missing, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    character*(*), intent(in)         :: fn
    integer, intent(in), optional     :: mode
    logical, intent(in), optional     :: add_missing
    integer, intent(inout), optional  :: ierror

    ! ---

    integer                       :: xtype, i
    integer(kind=MPI_OFFSET_KIND) :: nat, ndims, nframes
    integer                       :: dimids(NF90_MAX_VAR_DIMS)

    logical                       :: in_define_mode

    ! ---

    if (.not. initialized(p)) then
       ! Default initialization
       call init(p)
    endif

    this%mode  = F_READ
    if (present(mode)) then
       this%mode  = mode
    endif

    if (this%mode == F_READ) then
#ifdef _MP
       CHECK_NETCDF_ERROR( nf90mpi_open(mod_communicator%mpi%communicator, fn, NF90_NOWRITE, MPI_INFO_NULL, this%ncid), ierror )
#else
       CHECK_NETCDF_ERROR( nf90_open(fn, NF90_NOWRITE, this%ncid), ierror )
#endif

       if (present(add_missing) .and. add_missing) then
          RAISE_ERROR("Missing fields can only be added if in write mode.", ierror)
       endif
    else
#ifdef _MP
       CHECK_NETCDF_ERROR( nf90mpi_open(mod_communicator%mpi%communicator, fn, NF90_WRITE, MPI_INFO_NULL, this%ncid), ierror )
#else
       CHECK_NETCDF_ERROR( nf90_open(fn, NF90_WRITE, this%ncid), ierror )
#endif
    endif

    CHECK_NETCDF_ERROR( nf90_inq_dimid(this%ncid, NC_FRAME_STR, this%frame_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_inq_dimid(this%ncid, NC_SPATIAL_STR, this%spatial_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_inq_dimid(this%ncid, NC_ATOM_STR, this%atom_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_inq_dimid(this%ncid, NC_CELL_SPATIAL_STR, this%cell_spatial_dim), ierror )
    CHECK_NETCDF_ERROR( nf90_inq_dimid(this%ncid, NC_CELL_ANGULAR_STR, this%cell_angular_dim), ierror )
    if (nf90_inq_dimid(this%ncid, NC_LABEL_STR, this%label_dim) /= NF90_NOERR) then
       this%label_dim = -1
    endif

    CHECK_NETCDF_ERROR( nf90_inquire_dimension(this%ncid, this%frame_dim, len=nframes), ierror )
    this%nframes = nframes
    CHECK_NETCDF_ERROR( nf90_inquire_dimension(this%ncid, this%atom_dim, len=nat), ierror )

    CHECK_NETCDF_ERROR( nf90_inquire_dimension(this%ncid, this%spatial_dim, len=ndims), ierror )
    if (ndims /= 3) then
       RAISE_ERROR("Something wrong: Dimensions of " // NC_SPATIAL_STR // " /= 3.", ierror)
    endif

    CHECK_NETCDF_ERROR( nf90_inquire_dimension(this%ncid, this%cell_spatial_dim, len=ndims), ierror )
    if (ndims /= 3) then
       RAISE_ERROR("Something wrong: Dimensions of " // NC_CELL_SPATIAL_STR // " /= 3.", ierror)
    endif

    CHECK_NETCDF_ERROR( nf90_inquire_dimension(this%ncid, this%cell_angular_dim, len=ndims), ierror )
    if (ndims /= 3) then
       RAISE_ERROR("Something wrong: Dimensions of " // NC_CELL_ANGULAR_STR // " /= 3.", ierror)
    endif

    if (nf90_inq_varid(this%ncid, NC_TIME_STR, this%time_var) /= NF90_NOERR) then
       this%time_var = -1
    endif
    if (nf90_inq_varid(this%ncid, NC_CELL_ORIGIN_STR, this%cell_origin_var) /= NF90_NOERR) then
       this%cell_origin_var = -1
    endif
    CHECK_NETCDF_ERROR( nf90_inq_varid(this%ncid, NC_CELL_LENGTHS_STR, this%cell_lengths_var), ierror )
    CHECK_NETCDF_ERROR( nf90_inq_varid(this%ncid, NC_CELL_ANGLES_STR, this%cell_angles_var), ierror )
    if (nf90_inq_varid(this%ncid, NC_SHEAR_DX_STR, this%shear_dx_var) /= NF90_NOERR) then
       this%shear_dx_var = -1
    endif

    this%Z_var  = -1

    in_define_mode  = .false.

    if (p%data%n_real > 0) then
       allocate(this%real_var(p%data%n_real))
       allocate(this%real_ndims(p%data%n_real))

       this%real_var    = -1
       this%real_ndims  = -1

       do i = 1, p%data%n_real
          if (iand(p%data%tag_real(i), F_TO_TRAJ) /= 0) then
             if (nf90_inq_varid(this%ncid, trim(p%data%name_real(i)), this%real_var(i)) == NF90_NOERR) then
                CHECK_NETCDF_ERROR( nf90_inquire_variable(this%ncid, this%real_var(i), xtype=xtype, ndims=this%real_ndims(i), dimids=dimids), ierror )

                if (.not. (xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE)) then
                   RAISE_ERROR("Data type mismatch: Expected floating point type for field '" // trim(p%data%name_real(i)) // "'.", ierror)
                endif

                if (this%real_ndims(i) == 1) then
                   if (.not. (dimids(1) == this%atom_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_real(i)) // "'.", ierror)
                   endif
                else if (this%real_ndims(i) == 2) then
                   if (.not. (dimids(1) == this%atom_dim .and. dimids(2) == this%frame_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_real(i)) // "'.", ierror)
                   endif
                else
                   RAISE_ERROR("Data type mismatch: Wrong number of dimensions for field '" // trim(p%data%name_real(i)) // "'.", ierror)
                endif

             else
                if (present(add_missing) .and. add_missing) then

                   WARN("Field '" // trim(p%data%name_real(i)) // "' was not found in the NetCDF-file and will be added.")

                   if (.not. in_define_mode) then
                      CHECK_NETCDF_ERROR( nf90_redef( this%ncid ), ierror )
                      in_define_mode  = .true.
                   endif

                   if (iand(p%data%tag_real(i), F_CONSTANT) /= 0) then
                      CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real(i), NF90_FLOAT, (/ this%atom_dim /), this%real_var(i) ), ierror )
                   else
                      CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real(i), NF90_FLOAT, (/ this%atom_dim, this%frame_dim /), this%real_var(i) ), ierror )
                   endif

                   CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real_var(i), NC_UNITS_STR, p%data%unit_real(i)), ierror )
                   CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real_var(i), NC_SCALE_FACTOR_STR, p%data%conv_real(i)), ierror )

                else

                   WARN("Field '" // trim(p%data%name_real(i)) // "' not found in the NetCDF-file.")
                
                   this%real_var(i)  = -1

                endif
             endif
          endif
       enddo
    else

       this%real_var    => NULL()
       this%real_ndims  => NULL()

    endif

    if (p%data%n_integer > 0) then
       allocate(this%integer_var(p%data%n_integer))
       allocate(this%integer_ndims(p%data%n_integer))

       this%integer_var    = -1
       this%integer_ndims  = -1

       do i = 1, p%data%n_integer
          if (iand(p%data%tag_integer(i), F_TO_TRAJ) /= 0) then
             if (nf90_inq_varid(this%ncid, trim(p%data%name_integer(i)), this%integer_var(i)) /= NF90_NOERR) then
                if (trim(p%data%alias_integer(i)) /= "*") then
                   if (nf90_inq_varid(this%ncid, trim(p%data%alias_integer(i)), this%integer_var(i)) /= NF90_NOERR) then
                      this%integer_var(i) = -1
                   endif
                endif
             endif
             if (this%integer_var(i) >= 0) then
                CHECK_NETCDF_ERROR_WITH_INFO( nf90_inquire_variable(this%ncid, this%integer_var(i), xtype=xtype, ndims=this%integer_ndims(i), dimids=dimids),  "Variable " // p%data%name_integer(i) // ": ", ierror )

                if (.not. (xtype == NF90_INT)) then
                   RAISE_ERROR("Data type mismatch: Expected integer type for field '" // trim(p%data%name_integer(i)) // "'.", ierror)
                endif

                if (this%integer_ndims(i) == 1) then
                   if (.not. (dimids(1) == this%atom_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_integer(i)) // "'.", ierror)
                   endif
                else if (this%integer_ndims(i) == 2) then
                   if (.not. (dimids(1) == this%atom_dim .and. dimids(2) == this%frame_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_integer(i)) // "'.", ierror)
                   endif
                else
                   RAISE_ERROR("Data type mismatch: Wrong number of dimensions for field '" // trim(p%data%name_integer(i)) // "'.", ierror)
                endif

                if (trim(p%data%name_integer(i)) == trim(Z_STR)) then
                   this%Z_var  = this%integer_var(i)
                endif

             else
                if (present(add_missing) .and. add_missing) then

                   WARN("Field '" // trim(p%data%name_integer(i)) // "' was not found in the NetCDF-file and will be added.")

                   if (.not. in_define_mode) then
                      CHECK_NETCDF_ERROR( nf90_redef( this%ncid ), ierror )
                      in_define_mode  = .true.
                   endif

                   if (iand(p%data%tag_integer(i), F_CONSTANT) /= 0) then
                      CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_integer(i), NF90_INT, (/ this%atom_dim /), this%integer_var(i) ), ierror )

                      if (trim(p%data%name_integer(i)) == trim(Z_STR)) then
                         this%Z_var  = this%integer_var(i)
                      endif
                   else
                      CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_integer(i), NF90_INT, (/ this%atom_dim, this%frame_dim /), this%integer_var(i) ), ierror )
                   endif

                else

                   WARN("Field '" // trim(p%data%name_integer(i)) // "' not found in the NetCDF-file.")

                   this%integer_var(i)  = -1

                endif
             endif
          endif
       enddo
    else

       this%integer_var    => NULL()
       this%integer_ndims  => NULL()

    endif

    if (p%data%n_real3 > 0) then
       allocate(this%real3_var(p%data%n_real3))
       allocate(this%real3_ndims(p%data%n_real3))

       this%real3_var    = -1
       this%real3_ndims  = -1

       do i = 1, p%data%n_real3
          if (iand(p%data%tag_real3(i), F_TO_TRAJ) /= 0) then
             if (nf90_inq_varid(this%ncid, trim(p%data%name_real3(i)), this%real3_var(i)) == NF90_NOERR) then
                CHECK_NETCDF_ERROR( nf90_inquire_variable(this%ncid, this%real3_var(i), xtype=xtype, ndims=this%real3_ndims(i), dimids=dimids), ierror )

                if (.not. (xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE)) then
                   RAISE_ERROR("Data type mismatch: Expected floating point type for field '" // trim(p%data%name_real3(i)) // "'.", ierror)
                endif

                if (this%real3_ndims(i) == 2) then
                   if (.not. (dimids(1) == this%spatial_dim .and. dimids(2) == this%atom_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_real3(i)) // "'.", ierror)
                   endif
                else if (this%real3_ndims(i) == 3) then
                   if (.not. (dimids(1) == this%spatial_dim .and. dimids(2) == this%atom_dim .and. dimids(3) == this%frame_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_real3(i)) // "'.", ierror)
                   endif
                else
                   RAISE_ERROR("Data type mismatch: Wrong number of dimensions for field '" // trim(p%data%name_real3(i)) // "'.", ierror)
                endif

             else
                if (present(add_missing) .and. add_missing) then

                   WARN("Field '" // trim(p%data%name_real3(i)) // "' was not found in the NetCDF-file and will be added.")

                   if (.not. in_define_mode) then
                      CHECK_NETCDF_ERROR( nf90_redef( this%ncid ), ierror )
                      in_define_mode  = .true.
                   endif

                   if (iand(p%data%tag_real3(i), F_CONSTANT) /= 0) then
                      CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3(i), NF90_FLOAT, (/ this%spatial_dim, this%atom_dim /), this%real3_var(i) ), ierror )
                   else
                      CHECK_NETCDF_ERROR( nf90_def_var( this%ncid, p%data%name_real3(i), NF90_FLOAT, (/ this%spatial_dim, this%atom_dim, this%frame_dim /), this%real3_var(i) ), ierror )
                   endif

                   CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real3_var(i), NC_UNITS_STR, p%data%unit_real3(i)), ierror )
                   CHECK_NETCDF_ERROR( nf90_put_att(this%ncid, this%real3_var(i), NC_SCALE_FACTOR_STR, p%data%conv_real3(i)), ierror )

                else

                   WARN("Field '" // trim(p%data%name_real3(i)) // "' not found in the NetCDF-file.")

                   this%real3_var(i)  = -1

                endif
             endif
          endif
       enddo
    else

       this%real3_var    => NULL()
       this%real3_ndims  => NULL()

    endif

    if (p%data%n_real3x3 > 0) then
       allocate(this%real3x3_var(p%data%n_real3x3))
       allocate(this%real3x3_ndims(p%data%n_real3x3))

       this%real3x3_var    = -1
       this%real3x3_ndims  = -1

       do i = 1, p%data%n_real3x3
          if (iand(p%data%tag_real3x3(i), F_TO_TRAJ) /= 0) then
             if (nf90_inq_varid(this%ncid, trim(p%data%name_real3x3(i)), this%real3x3_var(i)) == NF90_NOERR) then
                CHECK_NETCDF_ERROR( nf90_inquire_variable(this%ncid, this%real3x3_var(i), xtype=xtype, ndims=this%real3x3_ndims(i), dimids=dimids), ierror )

                if (.not. (xtype == NF90_FLOAT .or. xtype == NF90_DOUBLE)) then
                   RAISE_ERROR("Data type mismatch: Expected floating point type for field '" // trim(p%data%name_real3x3(i)) // "'.", ierror)
                endif

                if (this%real3x3_ndims(i) == 3) then
                   if (.not. (dimids(1) == this%spatial_dim .and. dimids(2) == this%spatial_dim .and. dimids(3) == this%atom_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_real3x3(i)) // "'.", ierror)
                   endif
                else if (this%real3x3_ndims(i) == 4) then
                   if (.not. (dimids(1) == this%spatial_dim .and. dimids(2) == this%spatial_dim .and. dimids(3) == this%atom_dim .and. dimids(4) == this%frame_dim)) then
                      RAISE_ERROR("Data type mismatch: Wrong type of dimensions for field '" // trim(p%data%name_real3x3(i)) // "'.", ierror)
                   endif
                else
                   RAISE_ERROR("Data type mismatch: Wrong number of dimensions for field '" // trim(p%data%name_real3x3(i)) // "'.", ierror)
                endif

             else
                WARN("Field '" // trim(p%data%name_real3x3(i)) // "' not found in the NetCDF-file.")

                this%real3x3_var(i)  = -1
             endif
          endif
       enddo
    else

       this%real3x3_var    => NULL()
       this%real3x3_ndims  => NULL()

    endif

    if (in_define_mode) then
       CHECK_NETCDF_ERROR( nf90_enddef(this%ncid), ierror )
    endif

    if (.not. allocated(p)) then
       i = nat
       call allocate(p, i)
    else
       if (nat /= p%nat) then
          RAISE_ERROR("Particles object was allocated, however the number of particles does not match input file.", ierror)
       endif
    endif

    do i = 1, p%data%n_real
       if (this%real_ndims(i) == 1) then
          !CHECK_NETCDF_ERROR( nf90_get_var(this%ncid, this%real_var(i), p%data%data_real(:, i), start = (/ 1 /), count = (/ nat /)), ierror )
       endif
    enddo

    do i = 1, p%data%n_integer
       if (this%integer_ndims(i) == 1) then
          !CHECK_NETCDF_ERROR( nf90_get_var(this%ncid, this%integer_var(i), p%data%data_integer(:, i), start = (/ 1 /), count = (/ nat /)), ierror )
       endif
    enddo

    do i = 1, p%data%n_real3
       if (this%real3_ndims(i) == 2) then
          !CHECK_NETCDF_ERROR( nf90_get_var(this%ncid, this%real3_var(i), p%data%data_real3(:, :, i), start = (/ 1, 1 /), count = (/ 3, nat /)), ierror )
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (this%real3x3_ndims(i) == 3) then
          !CHECK_NETCDF_ERROR( nf90_get_var(this%ncid, this%real3x3_var(i), p%data%data_real3x3(:, :, :, i), start = (/ 1, 1, 1 /), count = (/ 3, 3, nat /)), ierror )
       endif
    enddo

    do i = 1, p%nat
       if (p%Z(i) > 0 .and. p%Z(i) <= MAX_Z) then
          p%sym(i)  = ElementName(p%Z(i))
          p%m(i)    = ElementMass(p%Z(i))
       else
          RAISE_ERROR("Unknown element encountered.", ierror)
       endif

       p%global2local(i) = i
       p%index(i)        = i
    enddo

    call update_elements(p)

    this%frame_no  = 1

    !
    ! Allocate buffers
    !

#ifndef _MP
    allocate(this%tmp_real(p%nat))
    allocate(this%tmp_integer(p%nat))
    allocate(this%tmp_real3(3, p%nat))
    allocate(this%tmp_real3x3(3, 3, p%nat))
#endif

  endsubroutine nc_open


  !>
  !! Close a NetCDF file
  !!
  !! Close a NetCDF file
  !<
  subroutine nc_close(this, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    integer, intent(inout), optional  :: ierror

    ! ---

    CHECK_NETCDF_ERROR( nf90_close(this%ncid), ierror )

    if (associated(this%real_var)) then
       deallocate(this%real_var)
       deallocate(this%real_ndims)
    endif

    if (associated(this%integer_var)) then
       deallocate(this%integer_var)
       deallocate(this%integer_ndims)
    endif

    if (associated(this%real3_var)) then
       deallocate(this%real3_var)
       deallocate(this%real3_ndims)
    endif

    if (associated(this%real3x3_var)) then
       deallocate(this%real3x3_var)
       deallocate(this%real3x3_ndims)
    endif

#ifndef _MP
    deallocate(this%tmp_real)
    deallocate(this%tmp_integer)
    deallocate(this%tmp_real3)
    deallocate(this%tmp_real3x3)
#endif

  endsubroutine nc_close


  !>
  !! Retrieve only the time information from a frame
  !!
  !! Retrieve only the time information from a frame
  !<
  real(DP) function nc_get_time(this, in_it, ierror)
    implicit none

    type(nc_t), intent(in)            :: this
    integer, intent(in)               :: in_it
    integer, intent(inout), optional  :: ierror

    ! ---

    integer(kind=MPI_OFFSET_KIND) :: it
    real(DP)                      :: ti

    ! ---

    if (in_it < 0) then
       it  = this%nframes + 1 + in_it
    else 
       it  = in_it
    endif

    if (this%time_var /= -1) then
       !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%time_var, ti, start = (/ it /)), "While reading frame " // it // ": ", ierror )
    else
       ti = it
    endif

    nc_get_time  = ti

  endfunction nc_get_time


  !>
  !! Find the frame that contains the time step \param ti
  !!
  !! Find the frame that contains the time step \param ti
  !<
  integer function nc_find_frame(this, ti, ierror)
    implicit none

    type(nc_t), intent(in)            :: this
    real(DP), intent(in)              :: ti
    integer, intent(inout), optional  :: ierror

    ! ---

    real(DP)  :: ti1, ti2, ti3
    integer   :: it1, it2, it3

    ! ---

    nc_find_frame  = -1

    it1  = 1
    it2  = this%nframes

    ti1  = get_time(this, it1, ierror)
    PASS_ERROR(ierror)
    ti2  = get_time(this, it2, ierror)
    PASS_ERROR(ierror)

    do while (it2-it1 > 1)
       it3  = (it1+it2)/2
       ti3  = get_time(this, it3, ierror)
       PASS_ERROR(ierror)

       if (ti3 > ti) then
          it2  = it3
       else
          it1  = it3
       endif
    enddo

    nc_find_frame  = it1

  endfunction nc_find_frame


  !>
  !! Read a frame
  !!
  !! Read a frame
  !<
  subroutine nc_read_frame(this, in_it, ti, p, ierror)
    implicit none

    type(nc_t), intent(in)            :: this
    integer, intent(in)               :: in_it
    real(DP), intent(out)             :: ti
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: ierror

    ! ---

    integer :: i, it

    real(DP) :: o(3), l(3), a(3), cell(3, 3), cx, cy, cz

    ! ---

    if (in_it < 0) then
       it  = this%nframes + 1 + in_it
    else 
       it  = in_it
    endif

    if (this%time_var /= -1) then
       !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%time_var, ti, start = (/ it /)), "While reading frame " // it // ": ", ierror )
    else
       ti = it
    endif

    do i = 1, p%data%n_real
       if (this%real_ndims(i) == 2) then
          !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%real_var(i), p%data%data_real(:, i), start = (/ 1, it /), count = (/ p%nat, 1 /)), "While reading frame " // it // ": ", ierror )
       endif
    enddo

    do i = 1, p%data%n_integer
       if (this%integer_ndims(i) == 2) then
          !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%integer_var(i), p%data%data_integer(:, i), start = (/ 1, it /), count = (/ p%nat, 1 /)), "While reading frame " // it // ": ", ierror )
       endif
    enddo

    do i = 1, p%data%n_real3
       if (this%real3_ndims(i) == 3) then
          !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%real3_var(i), p%data%data_real3(:, :, i), start = (/ 1, 1, it /), count = (/ 3, p%nat, 1 /)), "While reading frame " // it // ": ", ierror )
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (this%real3x3_ndims(i) == 4) then
          !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%real3x3_var(i), p%data%data_real3x3(:, :, :, i), start = (/ 1, 1, 1, it /), count = (/ 3, 3, p%nat, 1 /)), "While reading frame " // it // ": ", ierror )
       endif
    enddo

    o = [ 0.0_DP, 0.0_DP, 0.0_DP ]
    if (this%cell_origin_var > 0) then
       !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%cell_origin_var, o, start = (/ 1, it /), count = (/ 3, 1 /)), "While reading frame " // it // ": ", ierror )
       do i = 1, p%nat
#ifdef IMPLICIT_R
          PNC3(p, i) = PNC3(p, i) - o
#else
          POS3(p, i) = POS3(p, i) - o
#endif
       enddo
    endif

    !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%cell_lengths_var, l, start = (/ 1, it /), count = (/ 3, 1 /) ), "While reading frame " // it // ": ", ierror )
    !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%cell_angles_var, a, start = (/ 1, it /), count = (/ 3, 1 /) ), "While reading frame " // it // ": ", ierror )

    if (this%shear_dx_var > 0) then
       !CHECK_NETCDF_ERROR_WITH_INFO( nf90_get_var(this%ncid, this%shear_dx_var, p%shear_dx, start = (/ 1, it /), count = (/ 3, 1 /)), "While reading frame " // it // ": ", ierror )
    endif

    a = a*PI/180.0;
    cx = cos(a(2));
    cy = (cos(a(1)) - cos(a(2))*cos(a(3)))/sin(a(3));
    cz = sqrt(1.0_DP - cx*cx - cy*cy);
    cell(1:3, 1) = [ l(1), 0.0_DP, 0.0_DP ]
    cell(1:3, 2) = [ l(2)*cos(a(3)), l(2)*sin(a(3)), 0.0_DP ]
    cell(1:3, 3) = [ l(3)*cx, l(3)*cy, l(3)*cz ]

    call set_cell(p, cell, error=ierror)
    PASS_ERROR_WITH_INFO("While reading frame " // it // ".", ierror)

    do i = 1, p%nat
       PNC3(p, i) = POS3(p, i)
    enddo
    call inbox(p)
    call pnc2pos(p)

    call I_changed_positions(p)
    call I_changed_other(p)

  endsubroutine nc_read_frame


  !>
  !! Write constant information, i.e. groups, atomic numbers, etc. to the NetCDF file
  !!
  !! Write constant information, i.e. groups, atomic numbers, etc. to the NetCDF file
  !<
  subroutine nc_write_constants(this, p, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: ierror

    ! ---

    integer(kind=MPI_OFFSET_KIND), parameter :: one = 1, three = 3

    integer(kind=MPI_OFFSET_KIND) :: natloc
    integer                       :: i
#ifdef _MP
    integer(kind=MPI_OFFSET_KIND) :: start
#else
    integer                       :: j
#endif

    ! ---

    if (this%mode /= F_WRITE) then
       RAISE_ERROR("File has not been opened for write access.", ierror)
    endif

    natloc = p%natloc
#ifdef _MP
    start = cumsum(mod_communicator%mpi, p%natloc, error=ierror)-p%natloc+1
    PASS_ERROR(ierror)
#endif

    !
    ! Write global stuff
    !

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real_var(i), p%data%data_real(:, i), start = (/ start /), count = (/ natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_real(j) = p%data%data_real(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real_var(i), this%tmp_real, start = (/ one /), count = (/ natloc /) ), ierror )
#endif
          endif
       endif
    enddo

    do i = 1, p%data%n_integer
       if (iand(p%data%tag_integer(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_integer(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%integer_var(i), p%data%data_integer(:, i), start = (/ start /), count = (/ natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_integer(j)  = p%data%data_integer(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%integer_var(i), this%tmp_integer, start = (/ one /), count = (/ natloc /) ), ierror )
#endif
          endif
       endif
    enddo

    do i = 1, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real3(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3_var(i), p%data%data_real3(:, :, i), start = (/ one, start /), count = (/ three, natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_real3(:, j)  = p%data%data_real3(:, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3_var(i), this%tmp_real3, start = (/ one, one /), count = (/ three, natloc /) ), ierror )
#endif

          endif
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real3x3(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3x3_var(i), p%data%data_real3x3(:, :, :, i), start = (/ one, one, one /), count = (/ three, three, natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_real3x3(:, :, j)  = p%data%data_real3x3(:, :, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3x3_var(i), this%tmp_real3x3, start = (/ one, one, one /), count = (/ three, three, natloc /) ), ierror )
#endif

          endif
       endif
    enddo

  endsubroutine nc_write_constants


  !>
  !! Write a single constant field, i.e. groups, atomic numbers, etc. to the NetCDF file
  !!
  !! Write a single constant field, i.e. groups, atomic numbers, etc. to the NetCDF file
  !<
  subroutine nc_write_constant(this, p, field_name, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    character(*), intent(in)          :: field_name
    integer, intent(inout), optional  :: ierror

    ! ---

    integer(kind=MPI_OFFSET_KIND) :: one = 1, three = 3

    integer(kind=MPI_OFFSET_KIND) :: natloc
    integer                       :: i, n
#ifdef _MP
    integer(kind=MPI_OFFSET_KIND) :: start
#else
    integer                       :: j
#endif 

    ! ---

    if (this%mode /= F_WRITE) then
       RAISE_ERROR("File has not been opened for write access.", ierror)
    endif

    n = 0
    natloc = p%natloc
#ifdef _MP
    start = cumsum(mod_communicator%mpi, p%natloc, error=ierror)-p%natloc+1
    PASS_ERROR(ierror)
#endif

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_real(i)) == trim(field_name)) then
          if (iand(p%data%tag_real(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real_var(i), p%data%data_real(:, i), start = (/ start /), count = (/ natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_real(j)  = p%data%data_real(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real_var(i), this%tmp_real, start = (/ one /), count = (/ natloc /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    do i = 1, p%data%n_integer
       if (iand(p%data%tag_integer(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_integer(i)) == trim(field_name)) then
          if (iand(p%data%tag_integer(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%integer_var(i), p%data%data_integer(:, i), start = (/ start /), count = (/ natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_integer(j)  = p%data%data_integer(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%integer_var(i), this%tmp_integer, start = (/ one /), count = (/ natloc /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    do i = 1, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_real3(i)) == trim(field_name)) then
          if (iand(p%data%tag_real3(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3_var(i), p%data%data_real3(:, :, i), start = (/ one, start /), count = (/ three, natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_real3(:, j)  = p%data%data_real3(:, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3_var(i), this%tmp_real3, start = (/ one, one /), count = (/ three, natloc /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_real3x3(i)) == trim(field_name)) then
          if (iand(p%data%tag_real3x3(i), F_CONSTANT) /= 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3x3_var(i), p%data%data_real3x3(:, :, :, i), start = (/ one, one, start /), count = (/ three, three, natloc /) ), ierror )
#else
             do j = 1, natloc
                this%tmp_real3x3(:, :, j)  = p%data%data_real3x3(:, :, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3x3_var(i), this%tmp_real3x3, start = (/ one, one, one /), count = (/ three, three, natloc /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    !
    ! Sync
    !

    CHECK_NETCDF_ERROR( nf90_sync(this%ncid), ierror )

    if (n == 0) then
       RAISE_ERROR("Field '" // trim(field_name) // "' not found.", ierror)
    else if (n > 1) then
       RAISE_ERROR("Internal error: Field '" // trim(field_name) // "' seems to exist more than once.", ierror)
    endif

  endsubroutine nc_write_constant


  !>
  !! Write a complete frame
  !!
  !! Write a complete frame
  !<
  subroutine nc_write_frame(this, ti, p, frame_no, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    real(DP), intent(in)              :: ti
    type(particles_t), intent(inout)  :: p
    integer, intent(in), optional     :: frame_no
    integer, intent(inout), optional  :: ierror

    ! ---

    integer(kind=MPI_OFFSET_KIND) :: one = 1, three = 3

    integer(kind=MPI_OFFSET_KIND) :: natloc, fno
    integer                       :: i
#ifdef _MP
    integer(kind=MPI_OFFSET_KIND) :: start
#else
    integer                       :: j
#endif

    real(DP) :: time
    real(DP) :: cell(3,3)
    real(DP) :: cell_lengths(3)
    real(DP) :: cell_angles(3)

    ! ---

    if (this%mode /= F_WRITE) then
       RAISE_ERROR("File has not been opened for write access.", ierror)
    endif

    natloc = p%natloc
#ifdef _MP
    start = cumsum(mod_communicator%mpi, p%natloc, error=ierror)-p%natloc+1
    PASS_ERROR(ierror)
#endif

    fno = this%frame_no
    if (present(frame_no)) then
       fno = frame_no
    endif

    call get_true_cell(p, cell, error=ierror)
    PASS_ERROR(ierror)

    do i = 1, 3
       cell_lengths(i) = sqrt(dot_product(cell(:, i), cell(:, i)))
    enddo

    do i = 1, 3
       cell_angles(i) = acos( &
            dot_product(cell(:, mod(i, 3)+1), cell(:, mod(i+1, 3)+1))/(cell_lengths(mod(i, 3)+1)*cell_lengths(mod(i+1, 3)+1)) &
            ) * 180 / PI
    enddo

#ifdef _MP
    CHECK_NETCDF_ERROR( nf90mpi_begin_indep_data(this%ncid), ierror )
    if (mpi_id() == ROOT) then
#endif

    if (this%time_var /= -1) then
       time = ti
       CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%time_var, time, start = (/ fno /)), ierror )
    endif

    CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%cell_lengths_var, cell_lengths, start = (/ one, fno /), count = (/ three, one /)), ierror )
    CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%cell_angles_var, cell_angles, start = (/ one, fno /), count = (/ three, one /)), ierror )

!    CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%shear_dx_var, p%shear_dx, start = (/ 1, fno /), count = (/ three, one /)), ierror )

    !
    ! Write particle data
    !

    !
    ! Attributes
    !

    do i = 1, p%data%n_real_attr
       if (iand(p%data%tag_real_attr(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real_attr(i), F_CONSTANT) == 0) then

             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real_attr_var(i), p%data%data_real_attr(i:i+1), start = (/ fno /), count = (/ one /) ), ierror )

          endif
       endif
    enddo

    do i = 1, p%data%n_integer_attr
       if (iand(p%data%tag_integer_attr(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_integer_attr(i), F_CONSTANT) == 0) then

             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%integer_attr_var(i), p%data%data_integer_attr(i:i+1), start = (/ fno /), count = (/ one /) ), ierror )

          endif
       endif
    enddo

    do i = 1, p%data%n_real3_attr
       if (iand(p%data%tag_real3_attr(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real3_attr(i), F_CONSTANT) == 0) then

             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3_attr_var(i), p%data%data_real3_attr(:, i), start = (/ one, fno /), count = (/ three, one /) ), ierror )

          endif
       endif
    enddo

#ifdef _MP
    endif
    CHECK_NETCDF_ERROR( nf90mpi_end_indep_data(this%ncid), ierror )
#endif

    !
    ! Fields
    !

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real_var(i), p%data%data_real(:, i), start = (/ start, fno /), count = (/ natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_real(j)  = p%data%data_real(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real_var(i), this%tmp_real, start = (/ one, fno /), count = (/ natloc, one /) ), ierror )
#endif

          endif
       endif
    enddo

    do i = 1, p%data%n_integer
       if (iand(p%data%tag_integer(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_integer(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%integer_var(i), p%data%data_integer(:, i), start = (/ start, fno /), count = (/ natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_integer(j)  = p%data%data_integer(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%integer_var(i), this%tmp_integer, start = (/ one, fno /), count = (/ natloc, one /) ), ierror )
#endif

          endif
       endif
    enddo

    do i = 1, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real3(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3_var(i), p%data%data_real3(:, :, i), start = (/ one, start, fno /), count = (/ three, natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_real3(:, j)  = p%data%data_real3(:, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3_var(i), this%tmp_real3, start = (/ one, one, fno /), count = (/ three, natloc, one /) ), ierror )
#endif

          endif
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_TO_TRAJ) /= 0) then
          if (iand(p%data%tag_real3x3(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3x3_var(i), p%data%data_real3x3(:, :, :, i), start = (/ one, one, start, fno /), count = (/ three, three, natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_real3x3(:, :, j)  = p%data%data_real3x3(:, :, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3x3_var(i), this%tmp_real3x3, start = (/ one, one, one, fno /), count = (/ three, three, natloc, one /) ), ierror )
#endif

          endif
       endif
    enddo

    !
    ! Sync
    !

    CHECK_NETCDF_ERROR( nf90_sync(this%ncid), ierror )

    if (.not. present(frame_no)) then
       if (this%frame_no > this%nframes) then
          this%nframes  = this%frame_no
       endif

       this%frame_no = this%frame_no + 1
    endif

  endsubroutine nc_write_frame


  !>
  !! Write a single field
  !!
  !! Write a single field
  !<
  subroutine nc_write_field(this, p, fno, field_name, ierror)
    implicit none

    type(nc_t),                    intent(inout) :: this
    type(particles_t),             intent(inout) :: p
    integer(kind=MPI_OFFSET_KIND), intent(in)    :: fno
    character(*),                  intent(in)    :: field_name
    integer,             optional, intent(out)   :: ierror

    ! ---

    integer(kind=MPI_OFFSET_KIND), parameter :: one = 1, three = 3

    integer(kind=MPI_OFFSET_KIND) :: natloc
    integer                       :: i, n
#ifdef _MP
    integer(kind=MPI_OFFSET_KIND) :: start
#else
    integer                       :: j
#endif 

    ! ---

    INIT_ERROR(ierror)

    if (this%mode /= F_WRITE) then
       RAISE_ERROR("File has not been opened for write access.", ierror)
    endif

    n = 0
    natloc = p%natloc
#ifdef _MP
    start = cumsum(mod_communicator%mpi, p%natloc, error=ierror)-p%natloc+1
    PASS_ERROR(ierror)
#endif

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_real(i)) == trim(field_name)) then
          if (iand(p%data%tag_real(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real_var(i), p%data%data_real(:, i), start = (/ start, fno /), count = (/ natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_real(j)  = p%data%data_real(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real_var(i), this%tmp_real, start = (/ one, fno /), count = (/ natloc, one /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    do i = 1, p%data%n_integer
       if (iand(p%data%tag_integer(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_integer(i)) == trim(field_name)) then
          if (iand(p%data%tag_integer(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%integer_var(i), p%data%data_integer(:, i), start = (/ start, fno /), count = (/ natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_integer(j)  = p%data%data_integer(p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%integer_var(i), this%tmp_integer, start = (/ one, fno /), count = (/ natloc, one /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    do i = 1, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_real3(i)) == trim(field_name)) then
          if (iand(p%data%tag_real3(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real3_var(i), p%data%data_real3(:, :, i), start = (/ one, start, fno /), count = (/ three, natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_real3(:, j)  = p%data%data_real3(:, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3_var(i), this%tmp_real3, start = (/ one, one, fno /), count = (/ three, natloc, one /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_TO_TRAJ) /= 0 .and. trim(p%data%name_real3x3(i)) == trim(field_name)) then
          if (iand(p%data%tag_real3x3(i), F_CONSTANT) == 0) then

#ifdef _MP
             CHECK_NETCDF_ERROR( nf90mpi_put_var_all(this%ncid, this%real_var(i), p%data%data_real3x3(:, :, :, i), start = (/ one, one, start, fno /), count = (/ three, three, natloc, one /) ), ierror)
#else
             do j = 1, natloc
                this%tmp_real3x3(:, :, j)  = p%data%data_real3x3(:, :, p%global2local(j), i)
             enddo
             CHECK_NETCDF_ERROR( nf90_put_var(this%ncid, this%real3x3_var(i), this%tmp_real3x3, start = (/ one, one, one, fno /), count = (/ three, natloc, one /) ), ierror )
#endif

             n  = n + 1

          endif
       endif
    enddo

    !
    ! Sync
    !

    CHECK_NETCDF_ERROR( nf90_sync(this%ncid), ierror )

    if (n == 0) then
       RAISE_ERROR("Field '" // trim(field_name) // "' not found in internal data structure (does the field have the F_TO_TRAJ flag?).", ierror)
    else if (n > 1) then
       RAISE_ERROR("Internal error: Field '" // trim(field_name) // "' seems to exist more than once.", ierror)
    endif

  endsubroutine nc_write_field

endmodule nc

#else

module nc
  use supplib

  use particles

  !use netcdf

  use versioninfo

  implicit none

  character(*), parameter, private  :: MODULE_STR = "NC_DUMMY"

  character(*), parameter, private  :: NC_FRAME_STR         = "frame"
  character(*), parameter, private  :: NC_SPATIAL_STR       = "spatial"
  character(*), parameter, private  :: NC_ATOM_STR          = "atom"
  character(*), parameter, private  :: NC_CELL_SPATIAL_STR  = "cell_spatial"
  character(*), parameter, private  :: NC_CELL_ANGULAR_STR  = "cell_angular"
  character(*), parameter, private  :: NC_LABEL_STR         = "label"

  character(*), parameter, private  :: NC_TIME_STR          = "time"
  character(*), parameter, private  :: NC_CELL_LENGTHS_STR  = "cell_lengths"
  character(*), parameter, private  :: NC_CELL_ANGLES_STR   = "cell_angles"

  character(*), parameter, private  :: NC_SHEAR_DX_STR      = "shear_dx"

  character(*), parameter, private  :: NC_UNITS_STR         = "units"
  character(*), parameter, private  :: NC_SCALE_FACTOR_STR  = "scale_factor"

  type nc_t

     !
     ! Mode (read/write) and NetCDF file handle
     !

     integer               :: mode
     integer               :: ncid

     !
     ! Total number of frames in file and current
     ! frame for consecutive writes
     !

     integer               :: nframes
     integer               :: frame_no

     !
     ! Amber convention
     !

     integer               :: frame_dim
     integer               :: spatial_dim
     integer               :: atom_dim
     integer               :: cell_spatial_dim
     integer               :: cell_angular_dim
     integer               :: label_dim

     integer               :: spatial_var
     integer               :: cell_spatial_var
     integer               :: cell_angular_var

     integer               :: time_var
     integer               :: cell_lengths_var
     integer               :: cell_angles_var

     !
     ! MDCore convention
     !

     integer               :: Z_var
     integer               :: shear_dx_var

     !
     ! Dynamic fields
     !

     integer, pointer      :: real_attr_var(:)
     integer, pointer      :: real_attr_ndims(:)
     integer, pointer      :: integer_attr_var(:)
     integer, pointer      :: integer_attr_ndims(:)
     integer, pointer      :: real3_attr_var(:)
     integer, pointer      :: real3_attr_ndims(:)
     integer, pointer      :: real3x3_attr_var(:)
     integer, pointer      :: real3x3_attr_ndims(:)

     integer, pointer      :: real_var(:)
     integer, pointer      :: real_ndims(:)
     integer, pointer      :: integer_var(:)
     integer, pointer      :: integer_ndims(:)
     integer, pointer      :: real3_var(:)
     integer, pointer      :: real3_ndims(:)
     integer, pointer      :: real3x3_var(:)
     integer, pointer      :: real3x3_ndims(:)

     !
     ! Temporary buffers
     !

     real(DP), pointer      :: tmp_real(:)
     integer, pointer       :: tmp_integer(:)
     real(DP), pointer      :: tmp_real3(:, :)

  endtype nc_t

  interface create
     module procedure nc_create
  endinterface

  interface open
     module procedure nc_open
  endinterface

  interface close
     module procedure nc_close
  endinterface

  interface get_time
     module procedure nc_get_time
  endinterface

  interface find_frame
     module procedure nc_find_frame
  endinterface

  interface read_frame
     module procedure nc_read_frame
  endinterface

  interface write_constant
     module procedure nc_write_constant
  endinterface

  interface write_frame
     module procedure nc_write_frame
  endinterface

  interface write_field
     module procedure nc_write_field
  endinterface

contains

  subroutine write_prmtop(p, fn, ierror)
    implicit none

    type(particles_t), intent(in)  :: p
    character(*), intent(in)       :: fn
    integer, intent(inout), optional :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine write_prmtop


  subroutine nc_create(this, p, fn, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    character(*), intent(in)          :: fn
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_create

  subroutine nc_open(this, p, fn, mode, add_missing, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    character*(*), intent(in)         :: fn
    integer, intent(in), optional     :: mode
    logical, intent(in), optional     :: add_missing
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_open

  subroutine nc_close(this, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_close

  real(DP) function nc_get_time(this, in_it, ierror)
    implicit none

    type(nc_t), intent(in)            :: this
    integer, intent(in)               :: in_it
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)
    nc_get_time = 0.0_DP

  endfunction nc_get_time

  integer function nc_find_frame(this, ti, ierror)
    implicit none

    type(nc_t), intent(in)            :: this
    real(DP), intent(in)              :: ti
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)
    nc_find_frame = 0

  endfunction nc_find_frame

  subroutine nc_read_frame(this, in_it, ti, p, ierror)
    implicit none

    type(nc_t), intent(in)            :: this
    integer, intent(in)               :: in_it
    real(DP), intent(out)             :: ti
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)
    ti = 0.0_DP

  endsubroutine nc_read_frame

  subroutine nc_write_constants(this, p, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_write_constants

  subroutine nc_write_constant(this, p, field_name, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(inout)  :: p
    character(*), intent(in)          :: field_name
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_write_constant

  subroutine nc_write_frame(this, ti, p, frame_no, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    real(DP), intent(in)              :: ti
    type(particles_t), intent(in)     :: p
    integer, intent(in), optional     :: frame_no
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_write_frame

  subroutine nc_write_field(this, p, fno, field_name, ierror)
    implicit none

    type(nc_t), intent(inout)         :: this
    type(particles_t), intent(in)     :: p
    integer, intent(in)               :: fno
    character(*), intent(in)          :: field_name
    integer, intent(inout), optional  :: ierror

    ! ---

    RAISE_ERROR("Recompile with NetCDF support if you want to use NetCDF.", ierror)

  endsubroutine nc_write_field

endmodule nc
#endif
