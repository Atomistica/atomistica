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
!   classtype:output_time_t classname:OutputTime interface:callables
! @endmeta

!**********************************************************************
! Output the date, time and duration of time steps
!**********************************************************************

#include "macros.inc"

module output_time
  use libAtoms_module

  use io
  use logging

  use particles
  use dynamics
  use neighbors

  ! use main           ! just for "max_time", to calculate end date & time

  implicit none

  private

  public :: output_time_t
  type output_time_t

     real(DP)  :: freq  = -1.0_DP

     integer   :: un

     !
     ! Averaging
     !

     real(DP)  :: t

     real(DP)  :: begin             !< Starting time of the simulation
     real(DP)  :: last              !< The date and time of the last written time step

  endtype output_time_t

  ! ---

  public :: init
  interface init
     module procedure output_time_init
  endinterface

  public :: del
  interface del
     module procedure output_time_del
  endinterface

  public :: invoke
  interface invoke
     module procedure output_time_invoke
  endinterface

  public :: register
  interface register
    module procedure output_time_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine output_time_init(this)
    implicit none

    type(output_time_t), intent(inout)  :: this

    ! ---

    integer                             :: now(8)

    ! ---

#ifdef _MP
    if (mpi_id() == 0) then
#endif

    write (ilog, '(A)')            "- output_time_init -"
    write (ilog, '(5X,A,F20.10)')  "freq = ", this%freq

    this%un = fopen("time.out", F_WRITE)

    write (this%un, '(A1,1X,A7,1X,1A16,1A21,A16,A12,A12,A12,A16,A16)')  "#", "1:it", "2:time", "3:now", &
     "4:whole dur", "5:li dur", "6:av li dur", "7:av it dur", "8:1ns takes", "9:ends"

    this%t      = 0.0_DP

    ! At best, call here a function that returns the amount of seconds since 1970!!
    !call cpu_time(this%begin)
    !call PXFLOCALTIME(this%begin)
    !call localtime(this%begin)
    call date_and_time(values=now)
    this%begin = date_into_seconds(now)

    this%last = this%begin

    write (ilog, *)

#ifdef _MP
    endif
#endif

  endsubroutine output_time_init


  !>
  !! Destructor
  !!
  !! Delete a output_time object
  !<
  subroutine output_time_del(this)
    implicit none

    type(output_time_t), intent(inout)  :: this

    ! ---

#ifdef _MP
    if (mpi_id() == 0) then
#endif

    call fclose(this%un)

#ifdef _MP
    endif
#endif

  endsubroutine output_time_del


  !>
  !! Output the time data
  !!
  !! Output the time data
  !<
  subroutine output_time_invoke(this, dyn, nl, ierror)
    implicit none

    type(output_time_t), intent(inout)  :: this
    type(dynamics_t), intent(in)        :: dyn
    type(neighbors_t), intent(in)       :: nl
    integer, intent(inout), optional    :: ierror

    ! ---

    integer                             :: now(8)
    real(DP)                            :: nowcpu
    CHARACTER(LEN=8)                    :: mydate
    CHARACTER(LEN=10)                   :: mytime
    real(DP)                            :: realfreq

    ! ---

    if ((this%freq < 0) .or. (this%t >= this%freq)) then

      call date_and_time(date=mydate,time=mytime,values=now)

      ! "realfreq" is used for the calculation of the average durations of lines and steps
      if (this%freq .le. 0.0) then
       realfreq = dyn%dt
      else
       realfreq = this%freq
      end if

      ! At best, call here a function that gives seconds since 1970!!
      !call cpu_time(nowcpu)
      !call PXFLOCALTIME(nowcpu)
      !call localtime(nowcpu)
      nowcpu = date_into_seconds(now)

      !write (this%un, '(I9,X,ES20.10,4X,A8,X,A10,X)')  dyn%it, dyn%ti, mydate, mytime
      write (this%un, '(I9,X,ES16.9,2X,I2,A,I2,A,I4,X,I2,A,I2,A,I2,A16,A12,A12,A12,A16,A16)')  dyn%it, dyn%ti, &
       ! the actual time and date, when the line is written into the file
       now(3), '.', now(2), '.', now(1), now(5), ':', now(6), ':', now(7), &
       ! whole and line duration
       time_duration(this%begin, .true., nowcpu), time_duration(this%last, .false., nowcpu), &
       ! average line and step duration
       time_duration((nowcpu-this%begin)/MAX(1,dyn%it+1)/dyn%dt*realfreq), time_duration((nowcpu-this%begin)/MAX(1,dyn%it+1)), &
       ! how long takes 1ns to calculate
       time_duration(((1000000.0/10.18)/dyn%dt)*(nowcpu-this%begin)/MAX(1,dyn%it+1), .true.), &
       ! estimates the approximate time until the end of the simulation is reached
       !time_duration(((dyn%maxtime-(dyn%it+1)*dyn%dt)/dyn%dt)*(nowcpu-this%begin)/MAX(1,dyn%it+1), .true.) ! this starts at %it=0, but %ti>0 may be initially already
       time_duration((dyn%maxtime-dyn%ti)*(nowcpu-this%begin)/MAX(1,dyn%it+1)/dyn%dt, .true.)

      this%t    = 0.0_DP   ! Reset "%t", so that is has to reach "%freq" again before writing another line

      this%last = nowcpu   ! Save the date and time of this step: To calculate the duration until the next step

    endif

    this%t = this%t + dyn%dt

  endsubroutine output_time_invoke

  !>
  !! Convert seconds into time duration
  !!
  !! Convert seconds into time duration
  !<
  function time_duration(t1, days, t2)
    implicit none

    CHARACTER(LEN=16)                   :: time_duration
    real(DP), intent(in)                :: t1
    real(DP), intent(in), optional      :: t2
    logical, intent(in), optional       :: days

    ! ---

    integer, parameter                  :: ONESECOND = 1
    integer, parameter                  :: ONEMINUTE = 60 * ONESECOND
    integer, parameter                  :: ONEHOUR   = 60 * ONEMINUTE
    integer, parameter                  :: ONEDAY    = 24 * ONEHOUR

    ! ---

    real(DP)                            :: duration
    integer                             :: intduration ! type conversion for using "MOD" and "/"-operator
    logical                             :: usedays

    ! ---

    if (present(t2)) then
     duration = t2 - t1
    else
     duration = t1
    end if
    if (present(days)) then 
     usedays = days
    else
     usedays = .false.
    end if

    time_duration = ''

    duration = MAX(duration, 0.0)
    intduration = AnInt(duration)

    if (usedays) then
     write(time_duration, '(2X,I4,A,X,I2,A,I2,A,I2)') (intduration / ONEDAY), 'd', &
     (MOD(intduration, ONEDAY) / ONEHOUR), ':', &
     (MOD(intduration, ONEHOUR) / ONEMINUTE), ':', &
     (MOD(intduration, ONEMINUTE) / ONESECOND)
    else
     if (duration .lt. 10.0) then ! for small times: show milliseconds
      write(time_duration, '(4X,ES8.2)') duration
     else
      write(time_duration, '(4X,I2,A,I2,A,I2)') (MOD(intduration, ONEDAY) / ONEHOUR), ':', &
       (MOD(intduration, ONEHOUR) / ONEMINUTE), ':', &
       (MOD(intduration, ONEMINUTE) / ONESECOND)
     end if
    end if

    ! time_duration = trim(time_duration)

  end function time_duration

  !>
  !! Converts date and time into seconds
  !!
  !! Converts date and time into seconds (Be aware: It is not exact! There are hops when changing month or year, but fore approximations it is good)
  !<
  function date_into_seconds(dat)
   implicit none

   real(DP)                          :: date_into_seconds
   integer, intent(in)               :: dat(8)

   ! ---

   real(DP)                          :: help ! for non-integer calculations (here: for the milliseconds)

   ! the year "dat(1)" is skipped to make numbers smaller
   help = dat(8)
   help = help / 1000.0
   date_into_seconds = Real(dat(7)) + 60.0 * (Real(dat(6)) + 60.0 * (Real(dat(5)) + 24.0 * (Real(dat(3)) + 31.0 * Real(dat(2)))))
   date_into_seconds = date_into_seconds + help

  end function date_into_seconds


  subroutine output_time_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(output_time_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("OutputTime"), &
         CSTR("Output date, time and duration of time steps."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Output frequency (-1 means output every time step)."))

  endsubroutine output_time_register

endmodule output_time
