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
! @meta
!   shared
!   classtype:slicing_t classname:Slicing interface:callables
! @endmeta

!>
!! Compute averaged quantities over slices
!!
!! Compute averaged quantities over slices
!<

#include "macros.inc"

module slicing
  use supplib

  use particles
  use neighbors
  use dynamics

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  integer, parameter  :: n_dims       = 3
  integer, parameter  :: len_dim_str  = 15
  integer, parameter  :: ALL_DIMS     = 0

  ! This is need for xlf
  character(len_dim_str), parameter  :: STR_x            = CSTR("x")
  character(len_dim_str), parameter  :: STR_y            = CSTR("y")
  character(len_dim_str), parameter  :: STR_z            = CSTR("z")
  character(len_dim_str), parameter  :: dim_strs(n_dims) = &
       (/ STR_x, STR_y, STR_z /)

  integer, parameter  :: SL_AVERAGE    = 1
  integer, parameter  :: SL_HISTOGRAM  = 2

  public :: slicing_t
  type slicing_t

     logical            :: initialized = .false.

     !
     ! Flags
     !

     logical(BOOL)      :: compute_histograms = .false.

     !
     ! Geometry information
     !

     integer            :: n_bins

     integer            :: d
     integer            :: d2
     integer            :: d3
!     real(DP)           :: x1
!     real(DP)           :: x2

     real(DP)           :: smearing_length

     integer            :: n

     !
     ! Output frequency
     !

     real(DP)           :: freq
     real(DP)           :: ti
!     real(DP)           :: dx

     !
     ! Averages
     !

     type(Histogram1D)           :: spatial_avg_rho
     type(Histogram1D), pointer  :: spatial_avg_el(:)
     type(Histogram1D)           :: spatial_avg_T(3)

     real(DP), pointer           :: time_avg_rho(:)
     real(DP), pointer           :: time_avg_el(:, :)
     real(DP), pointer           :: time_avg_T(:, :)

     !
     ! Averages from generic particle data
     !

     integer                     :: n_out
     character(1000)             :: hdr_str(1000)

     type(Histogram1D), pointer  :: spatial_avg_real(:)
     type(Histogram1D), pointer  :: spatial_avg_real3(:, :)
     type(Histogram1D), pointer  :: spatial_avg_real3x3(:, :)

     real(DP), pointer           :: time_avg_real(:, :)
     real(DP), pointer           :: time_avg_real3(:, :, :)
     real(DP), pointer           :: time_avg_real3x3(:, :, :)

     real(DP), pointer           :: time_var_real(:, :)
     real(DP), pointer           :: time_var_real3(:, :, :)
     real(DP), pointer           :: time_var_real3x3(:, :, :)

     real(DP), pointer           :: time_avg_real_var(:, :)
     real(DP), pointer           :: time_avg_real3_var(:, :, :)
     real(DP), pointer           :: time_avg_real3x3_var(:, :, :)

     !
     ! Temperature histogram
     !

     integer            :: n_T_bins
     real(DP)           :: min_T
     real(DP)           :: max_T

     type(Histogram1D), pointer  :: T_hist(:, :)

     !
     ! Bond angle histogram
     !

     integer            :: n_angle_bins
     real(DP)           :: cutoff

     type(Histogram1D), pointer  :: angle_hist(:, :)

     !
     ! General histograms
     !

     logical, pointer   :: with_real(:)
     logical, pointer   :: with_real3(:)
     logical, pointer   :: with_real3x3(:)

     type(Histogram1D), pointer  :: hist_real(:, :)
     type(Histogram1D), pointer  :: hist_real3(:, :, :)
     type(Histogram1D), pointer  :: hist_real3x3(:, :, :)

     !
     ! Current velocities
     !

     real(DP), pointer  :: v(:, :)

  endtype slicing_t

  public :: init
  interface init
     module procedure slicing_init
  endinterface

  public :: del
  interface del
     module procedure slicing_del
  endinterface

  public :: invoke
  interface invoke
     module procedure slicing_invoke
  endinterface

  public :: register
  interface register
    module procedure slicing_register
  endinterface

contains

  !**********************************************************************
  ! Initialize a slicing object
  !**********************************************************************
  subroutine slicing_init(this)
    implicit none

    type(slicing_t), intent(inout)    :: this

    ! ---

    this%initialized = .false.

  endsubroutine slicing_init


  !**********************************************************************
  ! Initialize a slicing object
  !**********************************************************************
  subroutine slicing_internal_init(this, p, error)
    implicit none

    type(slicing_t), intent(inout)  :: this
    type(particles_t), intent(in)   :: p
    integer, optional, intent(out)  :: error

    ! ---

    integer                  :: i, j, un
    character(1000)          :: hlp_str

    integer, allocatable     :: n_bins_real(:)
    real(DP), allocatable    :: min_real(:)
    real(DP), allocatable    :: max_real(:)

    integer, allocatable     :: n_bins_real3(:)
    real(DP), allocatable    :: min_real3(:)
    real(DP), allocatable    :: max_real3(:)

    integer, allocatable     :: n_bins_real3x3(:)
    real(DP), allocatable    :: min_real3x3(:)
    real(DP), allocatable    :: max_real3x3(:)

    character(MAX_NAME_STR)  :: name

    ! ---

    INIT_ERROR(error)

    call prlog("- slicing_internal_init -")

    call log_memory_start("slicing_internal_init")

    call ptr_by_name(p%data, V_STR, this%v)

    this%d  = modulo(this%d-1, 3)+1
    this%d2 = modulo(this%d, 3)+1
    this%d3 = modulo(this%d+1, 3)+1

!    if (this%x1 < 0)   this%x1 = 0.0_DP
!    if (this%x2 < 0)   this%x2 = p%Abox(this%d, this%d)

!    this%dx = ( this%x2 - this%x1 ) / this%n_bins

!    write (ilog, '(5X,A,F20.10)')  "dx  = ", this%dx

    !
    ! Count number of (additional) output columns
    !

    this%n_out       = 1
    this%hdr_str(1)  = "i"
    this%hdr_str(2)  = "x"
    this%hdr_str(3)  = "density"
    this%hdr_str(4)  = "T(x)"
    this%hdr_str(5)  = "T(y)"
    this%hdr_str(6)  = "T(z)"
    this%hdr_str(7)  = "T"

    do i = 1, p%nel
       if (p%el2Z(i) > 0 .and. p%el2Z(i) <= MAX_Z) then
          write (hlp_str, '(1X,A)') trim(ElementName(p%el2Z(i)))
          this%hdr_str(this%n_out + 7) = hlp_str

          this%n_out = this%n_out + 1
       else
          RAISE_ERROR("Unknown element number encountered.", error)
       endif
    enddo

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_VERBOSE_ONLY) == 0) then
          write (hlp_str, '(A)') trim(p%data%name_real(i))
          this%hdr_str(this%n_out + 7)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real(i)), " - standard deviation"
          this%hdr_str(this%n_out + 8)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real(i)), " - standard deviation of the (spatial) average"
          this%hdr_str(this%n_out + 9)  = hlp_str

          this%n_out  = this%n_out + 3
       endif
    enddo

    ! Start from 2 to omit particle coordinates
    do i = 2, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_VERBOSE_ONLY) == 0) then
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(x)"
          this%hdr_str(this%n_out + 7)   = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(y)"
          this%hdr_str(this%n_out + 8)   = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(z)"
          this%hdr_str(this%n_out + 9)   = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(x) - average spatial variance"
          this%hdr_str(this%n_out + 10)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(y) - average spatial variance"
          this%hdr_str(this%n_out + 11)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(z) - average spatial variance"
          this%hdr_str(this%n_out + 12)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(x) - variance of the average"
          this%hdr_str(this%n_out + 13)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(y) - variance of the average"
          this%hdr_str(this%n_out + 14)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3(i)), "(z) - variance of the average"
          this%hdr_str(this%n_out + 15)  = hlp_str

          this%n_out  = this%n_out + 9
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_VERBOSE_ONLY) == 0) then
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(xx)"
          this%hdr_str(this%n_out + 7)   = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(yy)"
          this%hdr_str(this%n_out + 8)   = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(zz)"
          this%hdr_str(this%n_out + 9)   = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(xy)"
          this%hdr_str(this%n_out + 10)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(yz)"
          this%hdr_str(this%n_out + 11)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(zx)"
          this%hdr_str(this%n_out + 12)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(xx) - average spatial variance"
          this%hdr_str(this%n_out + 13)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(yy) - average spatial variance"
          this%hdr_str(this%n_out + 14)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(zz) - average spatial variance"
          this%hdr_str(this%n_out + 15)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(xy) - average spatial variance"
          this%hdr_str(this%n_out + 16)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(yz) - average spatial variance"
          this%hdr_str(this%n_out + 17)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(zx) - average spatial variance"
          this%hdr_str(this%n_out + 18)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(xx) - variance of the average"
          this%hdr_str(this%n_out + 19)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(yy) - variance of the average"
          this%hdr_str(this%n_out + 20)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(zz) - variance of the average"
          this%hdr_str(this%n_out + 21)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(xy) - variance of the average"
          this%hdr_str(this%n_out + 22)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(yz) - variance of the average"
          this%hdr_str(this%n_out + 23)  = hlp_str
          write (hlp_str, '(A,A)') trim(p%data%name_real3x3(i)), "(zx) - variance of the average"
          this%hdr_str(this%n_out + 24)  = hlp_str

          this%n_out  = this%n_out + 18
       endif
    enddo

    allocate(this%spatial_avg_el(p%nel))

    this%spatial_avg_real  => NULL()
    if (p%data%n_real > 0) then
       allocate(this%spatial_avg_real(p%data%n_real))
       call initialise(this%spatial_avg_real(:),       this%n_bins, 0.0_DP, p%Abox(this%d, this%d), this%smearing_length)
!       call log_memory_estimate(this%spatial_avg_real)

       allocate(this%time_avg_real(this%n_bins, p%data%n_real))
       allocate(this%time_var_real(this%n_bins, p%data%n_real))
       allocate(this%time_avg_real_var(this%n_bins, p%data%n_real))

!       call log_memory_estimate(this%time_avg_real)
!       call log_memory_estimate(this%time_var_real)
!       call log_memory_estimate(this%time_avg_real_var)

       this%time_avg_real(:, :)            = 0.0_DP
       this%time_var_real(:, :)            = 0.0_DP
       this%time_avg_real_var(:, :)        = 0.0_DP
    endif

    this%spatial_avg_real3  => NULL()
    if (p%data%n_real3 > 0) then
       allocate(this%spatial_avg_real3(3, p%data%n_real3))
       call initialise(this%spatial_avg_real3(:, :),   this%n_bins, 0.0_DP, p%Abox(this%d, this%d), this%smearing_length)
!       call log_memory_estimate(this%spatial_avg_real3)

       allocate(this%time_avg_real3(this%n_bins, 3, p%data%n_real3))
       allocate(this%time_var_real3(this%n_bins, 3, p%data%n_real3))
       allocate(this%time_avg_real3_var(this%n_bins, 3, p%data%n_real3))

!       call log_memory_estimate(this%time_avg_real3)
!       call log_memory_estimate(this%time_var_real3)
!       call log_memory_estimate(this%time_avg_real3_var)

       this%time_avg_real3(:, :, :)        = 0.0_DP
       this%time_var_real3(:, :, :)        = 0.0_DP
       this%time_avg_real3_var(:, :, :)    = 0.0_DP
    endif

    this%spatial_avg_real3x3  => NULL()
    if (p%data%n_real3x3 > 0) then
       allocate(this%spatial_avg_real3x3(6, p%data%n_real3x3))
       call initialise(this%spatial_avg_real3x3(:, :), this%n_bins, 0.0_DP, p%Abox(this%d, this%d), this%smearing_length)
!       call log_memory_estimate(this%spatial_avg_real3x3)

       allocate(this%time_avg_real3x3(this%n_bins, 6, p%data%n_real3x3))
       allocate(this%time_var_real3x3(this%n_bins, 6, p%data%n_real3x3))
       allocate(this%time_avg_real3x3_var(this%n_bins, 6, p%data%n_real3x3))

!       call log_memory_estimate(this%time_avg_real3x3)
!       call log_memory_estimate(this%time_var_real3x3)
!       call log_memory_estimate(this%time_avg_real3x3_var)

       this%time_avg_real3x3(:, :, :)      = 0.0_DP
       this%time_var_real3x3(:, :, :)      = 0.0_DP
       this%time_avg_real3x3_var(:, :, :)  = 0.0_DP
    endif

    call initialise(this%spatial_avg_rho,           this%n_bins, 0.0_DP, p%Abox(this%d, this%d), this%smearing_length)
    call initialise(this%spatial_avg_el(:),         this%n_bins, 0.0_DP, p%Abox(this%d, this%d), this%smearing_length)
    call initialise(this%spatial_avg_T(:),          this%n_bins, 0.0_DP, p%Abox(this%d, this%d), this%smearing_length)

!    call log_memory_estimate(this%spatial_avg_rho)
!    call log_memory_estimate(this%spatial_avg_el)
!    call log_memory_estimate(this%spatial_avg_T)

    allocate(this%time_avg_rho(this%n_bins))
    allocate(this%time_avg_el(this%n_bins, p%nel))
    allocate(this%time_avg_T(this%n_bins, 3))

!    call log_memory_estimate(this%time_avg_rho)
!    call log_memory_estimate(this%time_avg_el)
!    call log_memory_estimate(this%time_avg_T)

    this%time_avg_rho(:)                = 0.0_DP
    this%time_avg_el(:, :)              = 0.0_DP
    this%time_avg_T(:, :)               = 0.0_DP

    compute_histograms_1: if (this%compute_histograms) then
       this%with_real  => NULL()
       if (p%data%n_real > 0) then
          allocate(this%with_real(p%data%n_real))
          allocate(n_bins_real(p%data%n_real))
          allocate(min_real(p%data%n_real))
          allocate(max_real(p%data%n_real))

          this%with_real(:)     = .false.
          n_bins_real(:)        = 10
          min_real(:)           = -1.0_DP
          max_real(:)           = 1.0_DP
       endif

       this%with_real3  => NULL()
       if (p%data%n_real3 > 0) then
          allocate(this%with_real3(p%data%n_real3))
          allocate(n_bins_real3(p%data%n_real3))
          allocate(min_real3(p%data%n_real3))
          allocate(max_real3(p%data%n_real3))

          this%with_real3(:)    = .false.
          n_bins_real3(:)       = 10
          min_real3(:)          = -1.0_DP
          max_real3(:)          = 1.0_DP
       endif

       this%with_real3x3  => NULL()
       if (p%data%n_real3x3 > 0) then
          allocate(this%with_real3x3(p%data%n_real3x3))
          allocate(n_bins_real3x3(p%data%n_real3x3))
          allocate(min_real3x3(p%data%n_real3x3))
          allocate(max_real3x3(p%data%n_real3x3))

          this%with_real3x3(:)  = .false.
          n_bins_real3x3(:)     = 10
          min_real3x3(:)        = -1.0_DP
          max_real3x3(:)        = 1.0_DP
       endif

       un = fopen("histograms.dat", F_READ)
       if (un < 0) then
          RAISE_ERROR("Error opening 'histograms.dat'. Please provide that file.", error)
       endif
       read (un, *, iostat=i)  name
       do while (i == 0)
          j = 0
          if (p%data%n_real > 0) then
             j = index_by_name(p%data%n_real, p%data%name_real(:), name)
          endif
          if (j > 0) then
             this%with_real(j) = .true.
             read (un, *)  n_bins_real(j)
             read (un, *)  min_real(j)
             read (un, *)  max_real(j)

             write (hlp_str, '(A,A)') trim(p%data%name_real(j)), " - entropy"
             this%hdr_str(this%n_out + 7)  = hlp_str

             this%n_out  = this%n_out + 1
          else
             j = 0
             if (p%data%n_real3 > 0) then
                j = index_by_name(p%data%n_real3, p%data%name_real3(:), name)
             endif
             if (j > 0) then
                this%with_real3(j) = .true.
                read (un, *)  n_bins_real3(j)
                read (un, *)  min_real3(j)
                read (un, *)  max_real3(j)

                write (hlp_str, '(A,A)') trim(p%data%name_real3(j)), "(x) - entropy"
                this%hdr_str(this%n_out + 7)  = hlp_str
                write (hlp_str, '(A,A)') trim(p%data%name_real3(j)), "(y) - entropy"
                this%hdr_str(this%n_out + 8)  = hlp_str
                write (hlp_str, '(A,A)') trim(p%data%name_real3(j)), "(z) - entropy"
                this%hdr_str(this%n_out + 9)  = hlp_str

                this%n_out  = this%n_out + 3
             else
                j = 0
                if (p%data%n_real3x3 > 0) then
                   j = index_by_name(p%data%n_real3x3, p%data%name_real3x3(:), name)
                endif
                if (j > 0) then
                   this%with_real3x3(j) = .true.
                   read (un, *)  n_bins_real3x3(j)
                   read (un, *)  min_real3x3(j)
                   read (un, *)  max_real3x3(j)

                   write (hlp_str, '(A,A)') trim(p%data%name_real3x3(j)), "(xx) - entropy"
                   this%hdr_str(this%n_out + 7)   = hlp_str
                   write (hlp_str, '(A,A)') trim(p%data%name_real3x3(j)), "(yy) - entropy"
                   this%hdr_str(this%n_out + 8)   = hlp_str
                   write (hlp_str, '(A,A)') trim(p%data%name_real3x3(j)), "(zz) - entropy"
                   this%hdr_str(this%n_out + 9)   = hlp_str
                   write (hlp_str, '(A,A)') trim(p%data%name_real3x3(j)), "(xy) - entropy"
                   this%hdr_str(this%n_out + 10)  = hlp_str
                   write (hlp_str, '(A,A)') trim(p%data%name_real3x3(j)), "(yz) - entropy"
                   this%hdr_str(this%n_out + 11)  = hlp_str
                   write (hlp_str, '(A,A)') trim(p%data%name_real3x3(j)), "(zx) - entropy"
                   this%hdr_str(this%n_out + 12)  = hlp_str

                   this%n_out  = this%n_out + 6
                else
                   RAISE_ERROR("Unknown field: '" // trim(name) // "'.", error)
                endif
             endif
          endif
          read (un, *, iostat=i)  name
       enddo
       call fclose(un)

       if (this%n_T_bins > 0) then
          allocate(this%T_hist(3, this%n_bins))

          call initialise(this%T_hist(:, :), this%n_T_bins, this%min_T, this%max_T)

!          call log_memory_estimate(this%T_hist)
       endif

       if (this%n_angle_bins > 0) then
          allocate(this%angle_hist(this%n_bins, 3))

          call initialise(this%angle_hist(:, :), this%n_angle_bins, -90.0_DP, 90.0_DP, -1.0_DP, .true.)

!          call log_memory_estimate(this%angle_hist)

          write (hlp_str, '(A)') "bond_angles(xy) - entropy"
          this%hdr_str(this%n_out + 7)  = hlp_str
          write (hlp_str, '(A)') "bond_angles(yz) - entropy"
          this%hdr_str(this%n_out + 8)  = hlp_str
          write (hlp_str, '(A)') "bond_angles(zx) - entropy"
          this%hdr_str(this%n_out + 9)  = hlp_str

          this%n_out  = this%n_out + 3
       endif

       if (p%data%n_real > 0) then
          allocate(this%hist_real(this%n_bins, p%data%n_real))
       endif
       if (p%data%n_real3 > 0) then
          allocate(this%hist_real3(this%n_bins, 3, p%data%n_real3))
       endif
       if (p%data%n_real3x3 > 0) then
          allocate(this%hist_real3x3(this%n_bins, 6, p%data%n_real3x3))
       endif

       do j = 1, p%data%n_real
          if (this%with_real(j)) then
             write (ilog, '(5X,A40,A3,I10,2F20.10)')  p%data%name_real(j), " = ", n_bins_real(j), min_real(j), max_real(j)
             call initialise(this%hist_real(:, j), n_bins_real(j), min_real(j), max_real(j))
          endif
       enddo

       do j = 2, p%data%n_real3
          if (this%with_real3(j)) then
             write (ilog, '(5X,A40,A3,I10,2F20.10)')  p%data%name_real3(j), " = ", n_bins_real3(j), min_real3(j), max_real3(j)
             call initialise(this%hist_real3(:, :, j), n_bins_real3(j), min_real3(j), max_real3(j))
          endif
       enddo

       do j = 1, p%data%n_real3x3
          if (this%with_real3x3(j)) then
             write (ilog, '(5X,A40,A3,I10,2F20.10)')  p%data%name_real3x3(j), " = ", n_bins_real3x3(j), min_real3x3(j), max_real3x3(j)
             call initialise(this%hist_real3x3(:, :, j), n_bins_real3x3(j), min_real3x3(j), max_real3x3(j))
          endif
       enddo

       if (p%data%n_real > 0) then
          deallocate(n_bins_real)
          deallocate(min_real)
          deallocate(max_real)
       endif

       if (p%data%n_real3 > 0) then
          deallocate(n_bins_real3)
          deallocate(min_real3)
          deallocate(max_real3)
       endif

       if (p%data%n_real3x3 > 0) then
          deallocate(n_bins_real3x3)
          deallocate(min_real3x3)
          deallocate(max_real3x3)
       endif
    endif compute_histograms_1

    this%n_out  = this%n_out - 1

    call prlog("     n_out  = " // this%n_out)

    this%n   = 0
    this%ti  = 0.0_DP

    call log_memory_stop("slicing_internal_init")

    call prlog

  endsubroutine slicing_internal_init


  !**********************************************************************
  ! Delete a slicing object
  !**********************************************************************
  subroutine slicing_del(this)
    implicit none

    type(slicing_t), intent(inout)  :: this

    ! ---

    integer   :: j

    ! ---

    call finalise(this%spatial_avg_rho)
    call finalise(this%spatial_avg_el(:))
    call finalise(this%spatial_avg_T(:))

    deallocate(this%spatial_avg_el)

    if (associated(this%spatial_avg_real)) then
       call finalise(this%spatial_avg_real(:))

       deallocate(this%spatial_avg_real)
       deallocate(this%time_avg_real)
       deallocate(this%time_var_real)
       deallocate(this%time_avg_real_var)
    endif

    if (associated(this%spatial_avg_real3)) then
       call finalise(this%spatial_avg_real3(:, :))
       
       deallocate(this%spatial_avg_real3)
       deallocate(this%time_avg_real3)
       deallocate(this%time_var_real3)
       deallocate(this%time_avg_real3_var)
    endif

    if (associated(this%spatial_avg_real3x3)) then
       call finalise(this%spatial_avg_real3x3(:, :))

       deallocate(this%time_avg_real3x3)
       deallocate(this%time_var_real3x3)
       deallocate(this%time_avg_real3x3_var)
       deallocate(this%spatial_avg_real3x3)
    endif

    deallocate(this%time_avg_rho)
    deallocate(this%time_avg_el)
    deallocate(this%time_avg_T)

    compute_histograms_2: if (this%compute_histograms) then
       if (this%n_T_bins > 0) then
          call finalise(this%T_hist(:, :))
          deallocate(this%T_hist)
       endif

       if (this%n_angle_bins > 0) then
          call finalise(this%angle_hist(:, :))
          deallocate(this%angle_hist)
       endif

       do j = 1, size(this%hist_real, 2)
          if (this%with_real(j)) then
             call finalise(this%hist_real(:, j))
          endif
       enddo

       do j = 2, size(this%hist_real3, 3)
          if (this%with_real3(j)) then
             call finalise(this%hist_real3(:, :, j))
          endif
       enddo

       do j = 1, size(this%hist_real3x3, 3)
          if (this%with_real3x3(j)) then
             call finalise(this%hist_real3x3(:, :, j))
          endif
       enddo

       if (associated(this%with_real)) then
          deallocate(this%with_real)
          deallocate(this%hist_real)
       endif

       if (associated(this%with_real3)) then
          deallocate(this%with_real3)
          deallocate(this%hist_real3)
       endif

       if (associated(this%with_real3x3)) then
          deallocate(this%with_real3x3)
          deallocate(this%hist_real3x3)
       endif

    endif compute_histograms_2

  endsubroutine slicing_del


  !**********************************************************************
  ! Perform the measurement
  !**********************************************************************
  subroutine slicing_invoke(this, dyn, nl, error)
    implicit none

    type(slicing_t), intent(inout)  :: this
    type(dynamics_t), intent(in)    :: dyn
    type(neighbors_t), intent(in)   :: nl
    integer, optional, intent(out)  :: error

    ! ---

    integer          :: i, ni, j, k
    integer          :: bin1(dyn%p%natloc), bin2(dyn%p%natloc)

    real(DP)         :: x1(dyn%p%natloc), x2(dyn%p%natloc)
    real(DP)         :: r1(dyn%p%natloc), r2(dyn%p%natloc), r3(dyn%p%natloc)
    real(DP)         :: r_vol, dx
    real(DP)         :: theta, dr(3), abs_dr

    real(DP)         :: help(dyn%p%natloc)

    real(DP), allocatable  :: data(:)

    integer          :: un
    character(6)     :: fn
    character(1000)  :: fmt_str

    ! ---

    INIT_ERROR(error)

    if (.not. this%initialized) then
       call slicing_internal_init(this, dyn%p)
       this%initialized = .true.
    endif

    call timer_start("slicing_invoke")

    allocate(data(this%n_out))

    this%ti  = this%ti + dyn%dt

    dx       = dyn%p%Abox(this%d, this%d) / this%n_bins

    !$omp  parallel default(none) &
    !$omp& shared(dyn, this) &
    !$omp& private(help, i)

    !$omp sections
    !$omp section
    call set_bounds(this%spatial_avg_rho, 0.0_DP, dyn%p%Abox(this%d, this%d))
    call clear(this%spatial_avg_rho)
    call add(this%spatial_avg_rho, POS(dyn%p, 1:dyn%p%natloc, this%d), dyn%p%m(1:dyn%p%natloc))
    !$omp section
    help(1:dyn%p%natloc)  = dyn%p%m(1:dyn%p%natloc)*(VEC(this%v, 1:dyn%p%natloc, 1)*VEC(this%v, 1:dyn%p%natloc, 1))
    call set_bounds(this%spatial_avg_T(1), 0.0_DP, dyn%p%Abox(this%d, this%d))
    call clear(this%spatial_avg_T(1))
    call add(this%spatial_avg_T(1), POS(dyn%p, 1:dyn%p%natloc, this%d), help(1:dyn%p%natloc))
    !$omp section
    help(1:dyn%p%natloc)  = dyn%p%m(1:dyn%p%natloc)*(VEC(this%v, 1:dyn%p%natloc, 2)*VEC(this%v, 1:dyn%p%natloc, 2))
    call set_bounds(this%spatial_avg_T(2), 0.0_DP, dyn%p%Abox(this%d, this%d))
    call clear(this%spatial_avg_T(2))
    call add(this%spatial_avg_T(2), POS(dyn%p, 1:dyn%p%natloc, this%d), help(1:dyn%p%natloc))
    !$omp section
    help(1:dyn%p%natloc)  = dyn%p%m(1:dyn%p%natloc)*(VEC(this%v, 1:dyn%p%natloc, 3)*VEC(this%v, 1:dyn%p%natloc, 3))
    call set_bounds(this%spatial_avg_T(3), 0.0_DP, dyn%p%Abox(this%d, this%d))
    call clear(this%spatial_avg_T(3))
    call add(this%spatial_avg_T(3), POS(dyn%p, 1:dyn%p%natloc, this%d), help(1:dyn%p%natloc))
    !$omp endsections

    !$omp do
    do i = 1, dyn%p%nel
       call set_bounds(this%spatial_avg_el(i), 0.0_DP, dyn%p%Abox(this%d, this%d))
       call clear(this%spatial_avg_el(i))
       call add(this%spatial_avg_el(i), POS(dyn%p, 1:dyn%p%natloc, this%d), dyn%p%el(:) == i)
    enddo

    !$omp end parallel

!!    !$omp do
    do i = 1, dyn%p%data%n_real
       if (iand(dyn%p%data%tag_real(i), F_VERBOSE_ONLY) == 0) then
          call clear(this%spatial_avg_real(i))
          call set_bounds(this%spatial_avg_real(i), 0.0_DP, dyn%p%Abox(this%d, this%d))

          call add(this%spatial_avg_real(i), POS(dyn%p, 1:dyn%p%natloc, this%d), dyn%p%data%data_real(1:dyn%p%natloc, i))

#ifdef _MP
          call average(this%spatial_avg_real(i), mod_communicator%mpi)

          if (mod_communicator%mpi%my_proc == 0) then
#else
          call average(this%spatial_avg_real(i))
#endif

          this%time_avg_real(:, i)      = this%time_avg_real(:, i)     + this%spatial_avg_real(i)%h(:)*dyn%dt
          this%time_var_real(:, i)      = this%time_var_real(:, i)     + this%spatial_avg_real(i)%h(:)**2*dyn%dt
          this%time_avg_real_var(:, i)  = this%time_avg_real_var(:, i) + this%spatial_avg_real(i)%h_sq(:)*dyn%dt

#ifdef _MP
          endif
#endif
       endif
    enddo

!!    !$omp do
    do i = 2, dyn%p%data%n_real3
       if (iand(dyn%p%data%tag_real3(i), F_VERBOSE_ONLY) == 0) then
          call clear(this%spatial_avg_real3(:, i))
          call set_bounds(this%spatial_avg_real3(:, i), 0.0_DP, dyn%p%Abox(this%d, this%d))

          call add(this%spatial_avg_real3(1, i), POS(dyn%p, 1:dyn%p%natloc, this%d), dyn%p%data%data_real3(1, 1:dyn%p%natloc, i))
          call add(this%spatial_avg_real3(2, i), POS(dyn%p, 1:dyn%p%natloc, this%d), dyn%p%data%data_real3(2, 1:dyn%p%natloc, i))
          call add(this%spatial_avg_real3(3, i), POS(dyn%p, 1:dyn%p%natloc, this%d), dyn%p%data%data_real3(3, 1:dyn%p%natloc, i))

#ifdef _MP
          call average(this%spatial_avg_real3(1, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3(2, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3(3, i), mod_communicator%mpi)

          if (mod_communicator%mpi%my_proc == ROOT) then
#else
          call average(this%spatial_avg_real3(1, i))
          call average(this%spatial_avg_real3(2, i))
          call average(this%spatial_avg_real3(3, i))
#endif

          this%time_avg_real3(:, 1, i)      = this%time_avg_real3(:, 1, i)     + this%spatial_avg_real3(1, i)%h(:)*dyn%dt
          this%time_avg_real3(:, 2, i)      = this%time_avg_real3(:, 2, i)     + this%spatial_avg_real3(2, i)%h(:)*dyn%dt
          this%time_avg_real3(:, 3, i)      = this%time_avg_real3(:, 3, i)     + this%spatial_avg_real3(3, i)%h(:)*dyn%dt

          this%time_var_real3(:, 1, i)      = this%time_var_real3(:, 1, i)     + this%spatial_avg_real3(1, i)%h(:)**2*dyn%dt
          this%time_var_real3(:, 2, i)      = this%time_var_real3(:, 2, i)     + this%spatial_avg_real3(2, i)%h(:)**2*dyn%dt
          this%time_var_real3(:, 3, i)      = this%time_var_real3(:, 3, i)     + this%spatial_avg_real3(3, i)%h(:)**2*dyn%dt

          this%time_avg_real3_var(:, 1, i)  = this%time_avg_real3_var(:, 1, i) + this%spatial_avg_real3(1, i)%h_sq(:)*dyn%dt
          this%time_avg_real3_var(:, 2, i)  = this%time_avg_real3_var(:, 2, i) + this%spatial_avg_real3(2, i)%h_sq(:)*dyn%dt
          this%time_avg_real3_var(:, 3, i)  = this%time_avg_real3_var(:, 3, i) + this%spatial_avg_real3(3, i)%h_sq(:)*dyn%dt

#ifdef _MP
          endif
#endif
       endif
    enddo

!!    !$omp do
    do i = 1, dyn%p%data%n_real3x3
       if (iand(dyn%p%data%tag_real3x3(i), F_VERBOSE_ONLY) == 0) then
          call clear(this%spatial_avg_real3x3(:, i))
          call set_bounds(this%spatial_avg_real3x3(:, i), 0.0_DP, dyn%p%Abox(this%d, this%d))

          call add(this%spatial_avg_real3x3(1, i), POS(dyn%p, 1:dyn%p%natloc, this%d), &
               dyn%p%data%data_real3x3(1, 1, 1:dyn%p%natloc, i))
          call add(this%spatial_avg_real3x3(2, i), POS(dyn%p, 1:dyn%p%natloc, this%d), &
               dyn%p%data%data_real3x3(2, 2, 1:dyn%p%natloc, i))
          call add(this%spatial_avg_real3x3(3, i), POS(dyn%p, 1:dyn%p%natloc, this%d), &
               dyn%p%data%data_real3x3(3, 3, 1:dyn%p%natloc, i))
          call add(this%spatial_avg_real3x3(4, i), POS(dyn%p, 1:dyn%p%natloc, this%d), &
               ( dyn%p%data%data_real3x3(1, 2, 1:dyn%p%natloc, i) + dyn%p%data%data_real3x3(2, 1, 1:dyn%p%natloc, i) )/2)
          call add(this%spatial_avg_real3x3(5, i), POS(dyn%p, 1:dyn%p%natloc, this%d), &
               ( dyn%p%data%data_real3x3(2, 3, 1:dyn%p%natloc, i) + dyn%p%data%data_real3x3(3, 2, 1:dyn%p%natloc, i) )/2)
          call add(this%spatial_avg_real3x3(6, i), POS(dyn%p, 1:dyn%p%natloc, this%d), &
               ( dyn%p%data%data_real3x3(3, 1, 1:dyn%p%natloc, i) + dyn%p%data%data_real3x3(1, 3, 1:dyn%p%natloc, i) )/2)

#ifdef _MP
          call average(this%spatial_avg_real3x3(1, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3x3(2, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3x3(3, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3x3(4, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3x3(5, i), mod_communicator%mpi)
          call average(this%spatial_avg_real3x3(6, i), mod_communicator%mpi)

          if (mod_communicator%mpi%my_proc == 0) then
#else
          call average(this%spatial_avg_real3x3(1, i))
          call average(this%spatial_avg_real3x3(2, i))
          call average(this%spatial_avg_real3x3(3, i))
          call average(this%spatial_avg_real3x3(4, i))
          call average(this%spatial_avg_real3x3(5, i))
          call average(this%spatial_avg_real3x3(6, i))
#endif

          this%time_avg_real3x3(:, 1, i)      = this%time_avg_real3x3(:, 1, i)     + this%spatial_avg_real3x3(1, i)%h(:)*dyn%dt
          this%time_avg_real3x3(:, 2, i)      = this%time_avg_real3x3(:, 2, i)     + this%spatial_avg_real3x3(2, i)%h(:)*dyn%dt
          this%time_avg_real3x3(:, 3, i)      = this%time_avg_real3x3(:, 3, i)     + this%spatial_avg_real3x3(3, i)%h(:)*dyn%dt
          this%time_avg_real3x3(:, 4, i)      = this%time_avg_real3x3(:, 4, i)     + this%spatial_avg_real3x3(4, i)%h(:)*dyn%dt
          this%time_avg_real3x3(:, 5, i)      = this%time_avg_real3x3(:, 5, i)     + this%spatial_avg_real3x3(5, i)%h(:)*dyn%dt
          this%time_avg_real3x3(:, 6, i)      = this%time_avg_real3x3(:, 6, i)     + this%spatial_avg_real3x3(6, i)%h(:)*dyn%dt

          this%time_var_real3x3(:, 1, i)      = this%time_var_real3x3(:, 1, i)     + this%spatial_avg_real3x3(1, i)%h(:)**2*dyn%dt
          this%time_var_real3x3(:, 2, i)      = this%time_var_real3x3(:, 2, i)     + this%spatial_avg_real3x3(2, i)%h(:)**2*dyn%dt
          this%time_var_real3x3(:, 3, i)      = this%time_var_real3x3(:, 3, i)     + this%spatial_avg_real3x3(3, i)%h(:)**2*dyn%dt
          this%time_var_real3x3(:, 4, i)      = this%time_var_real3x3(:, 4, i)     + this%spatial_avg_real3x3(4, i)%h(:)**2*dyn%dt
          this%time_var_real3x3(:, 5, i)      = this%time_var_real3x3(:, 5, i)     + this%spatial_avg_real3x3(5, i)%h(:)**2*dyn%dt
          this%time_var_real3x3(:, 6, i)      = this%time_var_real3x3(:, 6, i)     + this%spatial_avg_real3x3(6, i)%h(:)**2*dyn%dt

          this%time_avg_real3x3_var(:, 1, i)  = this%time_avg_real3x3_var(:, 1, i) + this%spatial_avg_real3x3(1, i)%h_sq(:)*dyn%dt
          this%time_avg_real3x3_var(:, 2, i)  = this%time_avg_real3x3_var(:, 2, i) + this%spatial_avg_real3x3(2, i)%h_sq(:)*dyn%dt
          this%time_avg_real3x3_var(:, 3, i)  = this%time_avg_real3x3_var(:, 3, i) + this%spatial_avg_real3x3(3, i)%h_sq(:)*dyn%dt
          this%time_avg_real3x3_var(:, 4, i)  = this%time_avg_real3x3_var(:, 4, i) + this%spatial_avg_real3x3(4, i)%h_sq(:)*dyn%dt
          this%time_avg_real3x3_var(:, 5, i)  = this%time_avg_real3x3_var(:, 5, i) + this%spatial_avg_real3x3(5, i)%h_sq(:)*dyn%dt
          this%time_avg_real3x3_var(:, 6, i)  = this%time_avg_real3x3_var(:, 6, i) + this%spatial_avg_real3x3(6, i)%h_sq(:)*dyn%dt

#ifdef _MP
          endif
#endif

       endif
    enddo

!!    !$omp end parallel

    compute_histograms_3: if (this%compute_histograms) then

       x2(:)    = POS(dyn%p, :, this%d) / dx + 0.5_DP

       bin1(:)  = int(floor(x2(:)))
       bin2(:)  = bin1(:)+1
       x2(:)    = x2(:)-bin1(:)
       x1(:)    = 1.0_DP-x2(:)

       where (bin1 < 1)
          bin1 = bin1+this%n_bins
       endwhere
       where (bin2 > this%n_bins)
          bin2 = bin2-this%n_bins
       endwhere

       x1(:)    = x1(:)*dyn%dt
       x2(:)    = x2(:)*dyn%dt

       !$omp  parallel default(none) &
       !$omp& shared(bin1, bin2, dyn, nl, this, x1, x2) &
       !$omp& private(abs_dr, dr, i, ni, j, r1, r2, r3, theta)

       if (this%n_angle_bins > 0) then
          !$omp do
          do i = 1, dyn%p%natloc
             do ni = nl%seed(i), nl%last(i)
!                if (nl%abs_dr(ni) < this%cutoff) then
                DIST_SQ(dyn%p, nl, i, ni, dr, abs_dr)
                if (abs_dr < this%cutoff**2) then
!                   dr(:)  = VEC3(nl%dr, ni)
                   ! projected to the x-y plane
                   if (dr(1) /= 0.0_DP .or. dr(2) /= 0.0_DP) then
                      theta  = atan2(dr(2), dr(1)) * 180 / PI
                      call add(this%angle_hist(bin1(i), 1), theta, x1(i))
                      call add(this%angle_hist(bin2(i), 1), theta, x2(i))
                   endif
                   ! projected to the y-z plane
                   if (dr(2) /= 0.0_DP .or. dr(3) /= 0.0_DP) then
                      theta  = atan2(dr(3), dr(2)) * 180 / PI
                      call add(this%angle_hist(bin1(i), 2), theta, x1(i))
                      call add(this%angle_hist(bin2(i), 2), theta, x2(i))
                   endif
                   ! projected to the z-x plane
                   if (dr(1) /= 0.0_DP .or. dr(3) /= 0.0_DP) then
                      theta  = atan2(dr(1), dr(3)) * 180 / PI
                      call add(this%angle_hist(bin1(i), 3), theta, x1(i))
                      call add(this%angle_hist(bin2(i), 3), theta, x2(i))
                   endif
                endif
             enddo
          enddo
       endif

       !$omp do
       do j = 1, dyn%p%data%n_real
          if (this%with_real(j)) then
             do i = 1, dyn%p%natloc
                call add(this%hist_real(bin1(i), j), dyn%p%data%data_real(i, j), x1(i))
                call add(this%hist_real(bin2(i), j), dyn%p%data%data_real(i, j), x2(i))
             enddo
          endif
       enddo

       ! Omit positions
       !$omp do
       do j = 2, dyn%p%data%n_real3
          if (this%with_real3(j)) then
             do i = 1, dyn%p%natloc
                call add(this%hist_real3(bin1(i), 1, j), dyn%p%data%data_real3(1, i, j), x1(i))
                call add(this%hist_real3(bin2(i), 1, j), dyn%p%data%data_real3(1, i, j), x2(i))
                call add(this%hist_real3(bin1(i), 2, j), dyn%p%data%data_real3(2, i, j), x1(i))
                call add(this%hist_real3(bin2(i), 2, j), dyn%p%data%data_real3(2, i, j), x2(i))
                call add(this%hist_real3(bin1(i), 3, j), dyn%p%data%data_real3(3, i, j), x1(i))
                call add(this%hist_real3(bin2(i), 3, j), dyn%p%data%data_real3(3, i, j), x2(i))
             enddo
          endif
       enddo

       !$omp do
       do j = 1, dyn%p%data%n_real3x3
          if (this%with_real3x3(j)) then
             r1(1:dyn%p%natloc) = ( dyn%p%data%data_real3x3(1, 2, 1:dyn%p%natloc, j) + dyn%p%data%data_real3x3(2, 1, 1:dyn%p%natloc, j) )/2
             r2(1:dyn%p%natloc) = ( dyn%p%data%data_real3x3(2, 3, 1:dyn%p%natloc, j) + dyn%p%data%data_real3x3(3, 2, 1:dyn%p%natloc, j) )/2
             r3(1:dyn%p%natloc) = ( dyn%p%data%data_real3x3(3, 1, 1:dyn%p%natloc, j) + dyn%p%data%data_real3x3(1, 3, 1:dyn%p%natloc, j) )/2

             do i = 1, dyn%p%natloc
                call add(this%hist_real3x3(bin1(i), 1, j), dyn%p%data%data_real3x3(1, 1, i, j), x1(i))
                call add(this%hist_real3x3(bin2(i), 1, j), dyn%p%data%data_real3x3(1, 1, i, j), x2(i))
                call add(this%hist_real3x3(bin1(i), 2, j), dyn%p%data%data_real3x3(2, 2, i, j), x1(i))
                call add(this%hist_real3x3(bin2(i), 2, j), dyn%p%data%data_real3x3(2, 2, i, j), x2(i))
                call add(this%hist_real3x3(bin1(i), 3, j), dyn%p%data%data_real3x3(3, 3, i, j), x1(i))
                call add(this%hist_real3x3(bin2(i), 3, j), dyn%p%data%data_real3x3(3, 3, i, j), x2(i))

                call add(this%hist_real3x3(bin1(i), 4, j), r1(i), x1(i))
                call add(this%hist_real3x3(bin2(i), 4, j), r1(i), x2(i))
             
                call add(this%hist_real3x3(bin1(i), 5, j), r2(i), x1(i))
                call add(this%hist_real3x3(bin2(i), 5, j), r2(i), x2(i))

                call add(this%hist_real3x3(bin1(i), 6, j), r3(i), x1(i))
                call add(this%hist_real3x3(bin2(i), 6, j), r3(i), x2(i))
             enddo
          endif
       enddo

       !$omp end parallel

    endif compute_histograms_3

#ifdef _MP
    call average(this%spatial_avg_T(1), mod_communicator%mpi)
    call average(this%spatial_avg_T(2), mod_communicator%mpi)
    call average(this%spatial_avg_T(3), mod_communicator%mpi)

    if (mod_communicator%mpi%my_proc == 0) then
#else
    call average(this%spatial_avg_T(1))
    call average(this%spatial_avg_T(2))
    call average(this%spatial_avg_T(3))
#endif

    call mul(this%spatial_avg_T(1), 1.0_DP/K_to_energy)
    call mul(this%spatial_avg_T(2), 1.0_DP/K_to_energy)
    call mul(this%spatial_avg_T(3), 1.0_DP/K_to_energy)

#ifdef _MP
    endif
#endif

    compute_histograms_4: if (this%compute_histograms) then

       if (this%n_T_bins > 0) then
          do i = 1, this%n_bins
             do j = 1, 3
                call add(this%T_hist(j, i), this%spatial_avg_T(j)%h(i), dyn%dt)
             enddo
          enddo
       endif

    endif compute_histograms_4

    this%time_avg_rho(:)   = this%time_avg_rho(:)       + this%spatial_avg_rho%h(:)*dyn%dt
    do i = 1, dyn%p%nel
       this%time_avg_el(:, i)  = this%time_avg_el(:, i) + this%spatial_avg_el(i)%h(:)*dyn%dt
    enddo
    this%time_avg_T(:, 1)  = this%time_avg_T(:, 1)      + this%spatial_avg_T(1)%h(:)*dyn%dt
    this%time_avg_T(:, 2)  = this%time_avg_T(:, 2)      + this%spatial_avg_T(2)%h(:)*dyn%dt
    this%time_avg_T(:, 3)  = this%time_avg_T(:, 3)      + this%spatial_avg_T(3)%h(:)*dyn%dt

    if (this%ti >= this%freq) then

       this%n  = this%n + 1

       write (fn, '(I6.6)')  this%n

       !
       ! Write histograms to separate files
       !

       compute_histograms_5: if (this%compute_histograms) then

          if (this%n_T_bins > 0) then
             call write(reshape( this%T_hist(:, :), (/ 3*this%n_bins /) ), "T" // fn // ".out")
          endif

          if (this%n_angle_bins > 0) then
             call write( &
                  this%angle_hist(:, 1), &
                  "bond_angles_xy_" // fn // ".out" &
                  )
             call write( &
                  this%angle_hist(:, 2), &
                  "bond_angles_yz_" // fn // ".out" &
                  )
             call write( &
                  this%angle_hist(:, 3), &
                  "bond_angles_zx_" // fn // ".out" &
                  )
          endif

          do j = 1, dyn%p%data%n_real
             if (this%with_real(j)) then
                call write( &
                     this%hist_real(:, j), &
                     trim(dyn%p%data%name_real(j)) // "_" // fn // ".out" &
                     )
             endif
          enddo

          do j = 2, dyn%p%data%n_real3
             if (this%with_real3(j)) then
                call write( &
                     this%hist_real3(:, 1, j), &
                     trim(dyn%p%data%name_real3(j)) // "_x_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3(:, 2, j), &
                     trim(dyn%p%data%name_real3(j)) // "_y_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3(:, 3, j), &
                     trim(dyn%p%data%name_real3(j)) // "_z_" // fn // ".out" &
                     )
             endif
          enddo

          do j = 1, dyn%p%data%n_real3x3
             if (this%with_real3x3(j)) then
                call write( &
                     this%hist_real3x3(:, 1, j), &
                     trim(dyn%p%data%name_real3x3(j)) // "_xx_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3x3(:, 2, j), &
                     trim(dyn%p%data%name_real3x3(j)) // "_yy_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3x3(:, 3, j), &
                     trim(dyn%p%data%name_real3x3(j)) // "_zz_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3x3(:, 4, j), &
                     trim(dyn%p%data%name_real3x3(j)) // "_xy_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3x3(:, 5, j), &
                     trim(dyn%p%data%name_real3x3(j)) // "_yz_" // fn // ".out" &
                     )
                call write( &
                     this%hist_real3x3(:, 6, j), &
                     trim(dyn%p%data%name_real3x3(j)) // "_zx_" // fn // ".out" &
                     )
             endif
          enddo

       endif compute_histograms_5

       !
       ! Compute time averages
       !

       r_vol                   = 1.0_DP/(dyn%p%Abox(this%d2, this%d2)*dyn%p%Abox(this%d3, this%d3))
       this%time_avg_rho(:)    = r_vol * this%time_avg_rho(:) / this%ti
       this%time_avg_el(:, :)  = r_vol * this%time_avg_el(:, :) / this%ti

       this%time_avg_T(:, :)   = this%time_avg_T(:, :) / this%ti

       do i = 1, dyn%p%data%n_real
          if (iand(dyn%p%data%tag_real(i), F_VERBOSE_ONLY) == 0) then
             this%time_avg_real(:, i)      = this%time_avg_real(:, i)/this%ti
             this%time_var_real(:, i)      = this%time_var_real(:, i)/this%ti
             this%time_avg_real_var(:, i)  = this%time_avg_real_var(:, i)/this%ti
          endif
       enddo

       do i = 2, dyn%p%data%n_real3
          if (iand(dyn%p%data%tag_real3(i), F_VERBOSE_ONLY) == 0) then
             this%time_avg_real3(:, :, i)      = this%time_avg_real3(:, :, i)/this%ti
             this%time_var_real3(:, :, i)      = this%time_var_real3(:, :, i)/this%ti
             this%time_avg_real3_var(:, :, i)  = this%time_avg_real3_var(:, :, i)/this%ti
          endif
       enddo

       do i = 1, dyn%p%data%n_real3x3
          if (iand(dyn%p%data%tag_real3x3(i), F_VERBOSE_ONLY) == 0) then
             this%time_avg_real3x3(:, :, i)      = this%time_avg_real3x3(:, :, i)/this%ti
             this%time_var_real3x3(:, :, i)      = this%time_var_real3x3(:, :, i)/this%ti
             this%time_avg_real3x3_var(:, :, i)  = this%time_avg_real3x3_var(:, :, i)/this%ti
          endif
       enddo

       !
       ! Gather results from different processors
       !

#ifdef _MP
       call sum_in_place(mod_communicator%mpi, this%time_avg_rho)
       do i = 1, dyn%p%nel
          call sum_in_place(mod_communicator%mpi, this%time_avg_el(:, i))
       enddo
!!$       do i = 1, dyn%p%data%n_real
!!$          if (iand(dyn%p%data%tag_real(i), F_VERBOSE_ONLY) == 0) then
!!$             call dmp_sum_realarr(this%n_bins, this%time_avg_real(:, i), mod_communicator%mpi)
!!$             call dmp_sum_realarr(this%n_bins, this%time_var_real(:, i), mod_communicator%mpi)
!!$             call dmp_sum_realarr(this%n_bins, this%time_avg_real_var(:, i), mod_communicator%mpi)
!!$          endif
!!$       enddo
!!$       do i = 1, dyn%p%data%n_real3
!!$          if (iand(dyn%p%data%tag_real3(i), F_VERBOSE_ONLY) == 0) then
!!$             call dmp_sum_vecarr(this%n_bins, this%time_avg_real3(:, :, i), mod_communicator%mpi)
!!$             call dmp_sum_vecarr(this%n_bins, this%time_var_real3(:, :, i), mod_communicator%mpi)
!!$             call dmp_sum_vecarr(this%n_bins, this%time_avg_real3_var(:, :, i), mod_communicator%mpi)
!!$          endif
!!$       enddo
!!$       do i = 1, dyn%p%data%n_real3x3
!!$          EXIT_ON_ERROR("Implement for real3x3 and the variance.", i)
!!$       enddo

       if (mod_communicator%mpi%my_proc == ROOT) then
#endif

       !
       ! Normalize histograms for entropy calculation ONLY.
       ! This allows to later recombine differen histograms.
       !

       compute_histograms_6: if (this%compute_histograms) then

          if (this%n_T_bins > 0) then
             call normalize(this%T_hist(:, :))
          endif

          if (this%n_angle_bins > 0) then
             call normalize(this%angle_hist(:, :))
          endif

          do j = 1, dyn%p%data%n_real
             if (this%with_real(j)) then
                call normalize(this%hist_real(:, j))
             endif
          enddo

          do j = 2, dyn%p%data%n_real3
             if (this%with_real3(j)) then
                call normalize(this%hist_real3(:, :, j))
             endif
          enddo

          do j = 1, dyn%p%data%n_real3x3
             if (this%with_real3x3(j)) then
                call normalize(this%hist_real3x3(:, :, j))
             endif
          enddo

       endif compute_histograms_6

       !
       ! Write "slice_*.out"
       !

       write (fmt_str, '(A,I3.3,A)')  "(I10,", this%n_out+6, "ES20.10)"

       un = fopen("slice_" // fn // ".out", F_WRITE)
       write (un, '(A,ES20.10,A,ES20.10)')  "# t = ", dyn%ti-this%ti/2, ", dt = ", this%ti
       do i = 1, this%n_out + 7
          write (un, '(A,I3.3,A,A)')  "#   ", i, "   ", trim(this%hdr_str(i))
       enddo

       do i = 1, this%n_bins

          k = 1

          do j = 1, dyn%p%nel
             data(k)  = this%time_avg_el(i, j)
             k = k + 1
          enddo

          do j = 1, dyn%p%data%n_real
             if (iand(dyn%p%data%tag_real(j), F_VERBOSE_ONLY) == 0) then
                data(k)    = this%time_avg_real(i, j)
                data(k+1)  = sqrt( this%time_avg_real_var(i, j) - this%time_avg_real(i, j)**2 )
                data(k+2)  = sqrt( this%time_var_real(i, j) - this%time_avg_real(i, j)**2 )
                k = k + 3
             endif
          enddo

          do j = 2, dyn%p%data%n_real3
             if (iand(dyn%p%data%tag_real3(j), F_VERBOSE_ONLY) == 0) then
                data(k)    = this%time_avg_real3(i, 1, j)
                data(k+1)  = this%time_avg_real3(i, 2, j)
                data(k+2)  = this%time_avg_real3(i, 3, j)
                data(k+3)  = sqrt( this%time_avg_real3_var(i, 1, j) - this%time_avg_real3(i, 1, j)**2 )
                data(k+4)  = sqrt( this%time_avg_real3_var(i, 2, j) - this%time_avg_real3(i, 2, j)**2 )
                data(k+5)  = sqrt( this%time_avg_real3_var(i, 3, j) - this%time_avg_real3(i, 3, j)**2 )
                data(k+6)  = sqrt( this%time_var_real3(i, 1, j) - this%time_avg_real3(i, 1, j)**2 )
                data(k+7)  = sqrt( this%time_var_real3(i, 2, j) - this%time_avg_real3(i, 2, j)**2 )
                data(k+8)  = sqrt( this%time_var_real3(i, 3, j) - this%time_avg_real3(i, 3, j)**2 )
                k = k + 9
             endif
          enddo

          do j = 1, dyn%p%data%n_real3x3
             if (iand(dyn%p%data%tag_real3x3(j), F_VERBOSE_ONLY) == 0) then
                data(k)     = this%time_avg_real3x3(i, 1, j)
                data(k+1)   = this%time_avg_real3x3(i, 2, j)
                data(k+2)   = this%time_avg_real3x3(i, 3, j)
                data(k+3)   = this%time_avg_real3x3(i, 4, j)
                data(k+4)   = this%time_avg_real3x3(i, 5, j)
                data(k+5)   = this%time_avg_real3x3(i, 6, j)
                data(k+6)   = sqrt( this%time_avg_real3x3_var(i, 1, j) - this%time_avg_real3x3(i, 1, j)**2 )
                data(k+7)   = sqrt( this%time_avg_real3x3_var(i, 2, j) - this%time_avg_real3x3(i, 2, j)**2 )
                data(k+8)   = sqrt( this%time_avg_real3x3_var(i, 3, j) - this%time_avg_real3x3(i, 3, j)**2 )
                data(k+9)   = sqrt( this%time_avg_real3x3_var(i, 4, j) - this%time_avg_real3x3(i, 4, j)**2 )
                data(k+10)  = sqrt( this%time_avg_real3x3_var(i, 5, j) - this%time_avg_real3x3(i, 5, j)**2 )
                data(k+11)  = sqrt( this%time_avg_real3x3_var(i, 6, j) - this%time_avg_real3x3(i, 6, j)**2 )
                data(k+12)  = sqrt( this%time_var_real3x3(i, 1, j) - this%time_avg_real3x3(i, 1, j)**2 )
                data(k+13)  = sqrt( this%time_var_real3x3(i, 2, j) - this%time_avg_real3x3(i, 2, j)**2 )
                data(k+14)  = sqrt( this%time_var_real3x3(i, 3, j) - this%time_avg_real3x3(i, 3, j)**2 )
                data(k+15)  = sqrt( this%time_var_real3x3(i, 4, j) - this%time_avg_real3x3(i, 4, j)**2 )
                data(k+16)  = sqrt( this%time_var_real3x3(i, 5, j) - this%time_avg_real3x3(i, 5, j)**2 )
                data(k+17)  = sqrt( this%time_var_real3x3(i, 6, j) - this%time_avg_real3x3(i, 6, j)**2 )
                k = k + 18
             endif
          enddo

          compute_histograms_7: if (this%compute_histograms) then

             !
             ! Also output entropy of the histograms
             !

             do j = 1, dyn%p%data%n_real
                if (this%with_real(j)) then
                   data(k)  = entropy(this%hist_real(i, j))
                   k = k + 1
                endif
             enddo

             do j = 2, dyn%p%data%n_real3
                if (this%with_real3(j)) then
                   data(k)    = entropy(this%hist_real3(i, 1, j))
                   data(k+1)  = entropy(this%hist_real3(i, 2, j))
                   data(k+2)  = entropy(this%hist_real3(i, 3, j))
                   k = k + 3
                endif
             enddo
             
             do j = 1, dyn%p%data%n_real3x3
                if (this%with_real3x3(j)) then
                   data(k)    = entropy(this%hist_real3x3(i, 1, j))
                   data(k+1)  = entropy(this%hist_real3x3(i, 2, j))
                   data(k+2)  = entropy(this%hist_real3x3(i, 3, j))
                   data(k+3)  = entropy(this%hist_real3x3(i, 4, j))
                   data(k+4)  = entropy(this%hist_real3x3(i, 5, j))
                   data(k+5)  = entropy(this%hist_real3x3(i, 6, j))
                   k = k + 6
                endif
             enddo

             if (this%n_angle_bins > 0) then
                data(k)    = entropy(this%angle_hist(i, 1))
                data(k+1)  = entropy(this%angle_hist(i, 2))
                data(k+2)  = entropy(this%angle_hist(i, 3))
                k = k + 3
             endif
             
          endif compute_histograms_7

          write (un, fmt_str) &
               i, (i-0.5_DP)*dx, &
               this%time_avg_rho(i), &
               this%time_avg_T(i, :), sum(this%time_avg_T(i, :))/3, &
               data(:)
       enddo

       call fclose(un)

       !
       ! Clear histograms
       !

       compute_histograms_8: if (this%compute_histograms) then

          if (this%n_T_bins > 0) then
             call clear(this%T_hist(:, :))
          endif
          if (this%n_angle_bins > 0) then
             call clear(this%angle_hist(:, :))
          endif

          call clear(this%hist_real(:, :))
          call clear(this%hist_real3(:, :, :))
          call clear(this%hist_real3x3(:, :, :))

       endif compute_histograms_8

#ifdef _MP
       endif
#endif

       this%time_avg_T(:, :)              = 0.0_DP
       this%time_avg_rho(:)               = 0.0_DP
       this%time_avg_el(:, :)             = 0.0_DP
       if (dyn%p%data%n_real > 0) then
          this%time_avg_real(:, :)        = 0.0_DP
       endif
       if (dyn%p%data%n_real3 > 0) then
          this%time_avg_real3(:, :, :)    = 0.0_DP
       endif
       if (dyn%p%data%n_real3x3 > 0) then
          this%time_avg_real3x3(:, :, :)  = 0.0_DP
       endif

       this%ti                            = 0.0_DP

    endif

    deallocate(data)

    call timer_stop("slicing_invoke")

  endsubroutine slicing_invoke


  !****************************************************************
  ! Initialize the property list
  !****************************************************************
  subroutine slicing_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(slicing_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)         :: cfg
    type(c_ptr), intent(out)        :: m

    ! ---

    this%compute_histograms  = .false.
    this%n_bins              = 100
    this%freq                = -1.0_DP
    this%d                   = 3
!    this%x1      = -1
!    this%x2      = -1

    this%n_T_bins            = -1
    this%min_T               = 0.0_DP
    this%max_T               = 1000.0_DP

    this%n_angle_bins        = -1
    this%cutoff              = 1.85_DP

    this%smearing_length     = -1.0_DP

    m = ptrdict_register_section(cfg, CSTR("Slicing"), &
         CSTR("Determine space averaged quantities."))

    call ptrdict_register_boolean_property(m, c_loc(this%compute_histograms), CSTR("compute_histograms"), &
         CSTR("Compute histograms for temperature and velocity."))

    call ptrdict_register_enum_property(m, c_loc(this%d), &
         n_dims, len_dim_str, dim_strs(:), &
         CSTR("d"), &
         CSTR("Direction in which to slice: 'x', 'y', 'z'"))

!    call ptrdict_register_real_property(m, this%x1, CSTR("x1"), &
!         CSTR("Lower bound."))
!    call ptrdict_register_real_property(m, this%x2, CSTR("x2"), &
!         CSTR("Upper bound"))

    call ptrdict_register_integer_property(m, c_loc(this%n_bins), CSTR("n_bins"), &
         CSTR("Number of bins."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Time average (in real time)."))

    call ptrdict_register_real_property(m, c_loc(this%smearing_length), CSTR("sigma"), &
         CSTR("Interpolation length: If > 0 a Gaussian with this width is used for interpolation."))

    call ptrdict_register_integer_property(m, c_loc(this%n_T_bins), CSTR("n_T_bins"), &
         CSTR("Number of bins for T-histogram (do not compute if <= 0)."))
    call ptrdict_register_real_property(m, c_loc(this%min_T), CSTR("min_T"), &
         CSTR("Lower bound for T-histogram."))
    call ptrdict_register_real_property(m, c_loc(this%max_T), CSTR("max_T"), &
         CSTR("Upper bound for T-histogram."))

    call ptrdict_register_integer_property(m, c_loc(this%n_angle_bins), CSTR("n_angle_bins"), &
         CSTR("Number of bins for bond angle histogram (do not compute if <= 0)."))
    call ptrdict_register_real_property(m, c_loc(this%cutoff), CSTR("cutoff"), &
         CSTR("Cut-off for determination of bonds."))

  endsubroutine slicing_register

endmodule slicing
