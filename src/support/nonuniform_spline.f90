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

!**********************************************************************
! Spline helper functions
!**********************************************************************

#include "macros.inc"

module nonuniform_spline
  use error_module
  use system_module
  use io

  implicit none

  ! ---

  private

  public :: spline_t
  type spline_t

     integer            :: n = -1     ! number of entries
     integer            :: ncol       ! number of columns
     real(DP)           :: cut        ! end point

     logical            :: associated = .false.

     real(DP), pointer  :: x(:)       ! x-positions
     real(DP), pointer  :: y(:, :)    ! value tables
     real(DP), pointer  :: d2y(:, :)  ! second derivatives table

  endtype spline_t


  public :: init, nonuniform_spline_init
  interface init
     module procedure nonuniform_spline_init
  endinterface

  public :: del
  interface del
     module procedure nonuniform_spline_del
  endinterface

  public :: associate
  interface associate
     module procedure nonuniform_spline_associate
  endinterface

  public :: read
  interface read
     module procedure nonuniform_spline_read, nonuniform_spline_read_fn
  endinterface

  public :: read2
  interface read2
     module procedure nonuniform_spline_read2, nonuniform_spline_read2_fn
  endinterface

  public :: write
  interface write
     module procedure nonuniform_spline_write
  endinterface write

  public :: interval
  interface interval
     module procedure nonuniform_spline_interval
  endinterface

  public :: f
  interface f
     module procedure nonuniform_spline_f, nonuniform_spline_f_unknown_interval
  endinterface

  public :: df
  interface df
     module procedure nonuniform_spline_df
  endinterface

  public :: f_and_df
  interface f_and_df
     module procedure nonuniform_spline_f_and_df
  endinterface

contains

  !**********************************************************************
  ! Initialize the spline table
  !**********************************************************************
  subroutine nonuniform_spline_init(this, nmax, n, x, ncol, y)
    implicit none

    type(spline_t), intent(out)  :: this
    integer, intent(in)          :: nmax
    integer, intent(in)          :: n
    real(DP), intent(in)         :: x(nmax)
    integer, intent(in)          :: ncol
    real(DP), intent(in)         :: y(nmax, ncol)

    ! ---

    integer  :: i, k
    real(DP) :: p, qn(ncol), un(ncol), u(n, ncol), sig

    ! ---

    this%associated = .false.

    this%n     = n
    this%ncol  = ncol

    allocate(this%x(n))
    allocate(this%y(n, ncol))
    allocate(this%d2y(n, ncol))

    this%x(1:n)    = x(1:n)
    this%y(1:n, :) = y(1:n, :)

    this%d2y(1, :) = 0.d0 ! natural spline, i.e. second derivatives vanishes
    u(1, :)     = 0.d0

    do i = 2, n-1
       sig  = (this%x(i)-this%x(i-1))/(this%x(i+1)-this%x(i-1))
       do k = 1, ncol
          p           = sig*this%d2y(i-1, k) + 2
          this%d2y(i, k) = (sig-1)/p   !-p/2
          !       u(i)       = 3*( tab(i+1)+tab(i-1)-2*tab(i) )/dx**2 - u(i-1)/(2*p)
          u(i, k)     = (6d0*((this%y(i+1, k)-this%y(i, k))/(this%x(i+1)-this%x(i)) &
                -(this%y(i, k)-this%y(i-1, k))/(this%x(i)-this%x(i-1)))/(this%x(i+1)-this%x(i-1))-sig*u(i-1, k))/p
       enddo
    enddo

!    qn = 0d0 !natural spline
!    un = 0d0
    qn(:)  = 0.5d0
    un(:)  = (3./(this%x(n)-this%x(n-1)))*(0.0-(this%y(n, :)-this%y(n-1, :))/(this%x(n)-this%x(n-1)))

    this%d2y(n, :) = (un(:)-qn(:)*u(n-1, :))/(qn(:)*this%d2y(n-1, :)+1d0)

    do i = 1, ncol
       do k = n-1, 1, -1
          this%d2y(k, i) = this%d2y(k, i)*this%d2y(k+1, i) + u(k, i)
       enddo
    enddo

    this%cut = this%x(n)

  endsubroutine nonuniform_spline_init


  !**********************************************************************
  ! Initialize the spline table
  !**********************************************************************
  subroutine nonuniform_spline_associate(this, w)
    type(spline_t), intent(inout) :: this
    type(spline_t), intent(in)    :: w


    ! --

    this%associated = .true.

    this%n = w%n
    this%ncol = w%ncol
    this%cut = w%cut

    this%x => w%x
    this%y => w%y
    this%d2y => w%d2y

  endsubroutine nonuniform_spline_associate


  !**********************************************************************
  ! Delete a spline table
  !**********************************************************************
  subroutine nonuniform_spline_del(this)
    implicit none

    type(spline_t), intent(inout)  :: this

    ! ---

    if (.not. this%associated) then
       deallocate(this%x)
       deallocate(this%y)
       deallocate(this%d2y)
    endif

  endsubroutine nonuniform_spline_del


  !**********************************************************************
  ! Read a spline table from a file
  ! The format has to be the following
  !    1:   N                           Number of grid points
  !   2-:   x  s1  s2  s3  ...          Values
  ! *this* is the spline object,
  ! *un* the file, *ncol* the total number of columns in the file.
  ! If *ndata* is present the number of grid points is not read
  ! from the file.
  ! *xconv* and *yconv* can be used to convert the values.
  !**********************************************************************
  subroutine nonuniform_spline_read(this, un, ncol, xconv, yconv, ndata, xdata, ierror)
    implicit none

    type(spline_t), intent(out)       :: this
    integer, intent(in)               :: un
    integer, intent(in)               :: ncol
    real(DP), intent(in), optional    :: xconv
    real(DP), intent(in), optional    :: yconv(ncol-1)
    integer, intent(in), optional     :: ndata
    real(DP), intent(in), optional    :: xdata(*)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer                :: i, n, info
    real(DP), allocatable  :: x(:)
    real(DP), allocatable  :: data(:, :)

    ! ---

    if (ncol < 2) then
       RAISE_ERROR("Number of columns in input file is smaller than two.", ierror)
    endif

    if (present(ndata)) then
       n = ndata
    else
       read (un, *, iostat=info)  n
       if (info /= 0) then
          RAISE_ERROR("Error reading number of lines to follow.", ierror)
       endif
    endif

    allocate(x(n))
    allocate(data(n, ncol-1))

    x     = 0.0_DP
    data  = 0.0_DP

    do i = 1, n
       if (present(xdata)) then
          x(i) = xdata(i)
          read (un, *, iostat=info)  data(i, :)
       else
          read (un, *, iostat=info)  x(i), data(i, :)
       endif
       if (info /= 0) then
          if (i == 1) then
             RAISE_ERROR("Error reading " // i // "th line (of " // n // " lines) of data.", ierror)
          else
             RAISE_ERROR("Error reading " // i // "th line (of " // n // " lines) of data. Last data: " // data(i-1, :), ierror)
          endif
       endif

       if (present(xconv))  x(i)        = x(i)*xconv
       if (present(yconv))  data(i, :)  = data(i, :)*yconv
    enddo

    call nonuniform_spline_init(this, n, n, x, ncol-1, data)

    deallocate(x)
    deallocate(data)

  endsubroutine nonuniform_spline_read


  !**********************************************************************
  ! Read with filename given
  !**********************************************************************
  subroutine nonuniform_spline_read_fn(this, fn, ncol, xconv, yconv, ierror)
    implicit none

    type(spline_t), intent(out)       :: this
    character(*), intent(in)          :: fn
    integer, intent(in)               :: ncol
    real(DP), intent(in), optional    :: xconv
    real(DP), intent(in), optional    :: yconv(ncol-1)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: un

    ! ---

    un = fopen(fn, F_READ, ierror)
    PASS_ERROR(ierror)
    call nonuniform_spline_read(this, un, ncol, xconv, yconv, ierror)
    PASS_ERROR(ierror)
    call fclose(un, ierror)
    PASS_ERROR(ierror)

  endsubroutine nonuniform_spline_read_fn


  !**********************************************************************
  ! Read a spline table from a file
  !   1-:   x  s1  s2  s3  ...          Values
  !    X:  empty line
  ! *this* is the spline object,
  ! *un* the file, *ncol* the total number of columns in the file.
  ! *xconv* and *yconv* can be used to convert the values.
  !**********************************************************************
  subroutine nonuniform_spline_read2(this, un, ncol, xconv, yconv, ierror)
    implicit none

    type(spline_t), intent(out)       :: this
    integer, intent(in)               :: un
    integer, intent(in)               :: ncol
    real(DP), intent(in), optional    :: xconv
    real(DP), intent(in), optional    :: yconv(ncol-1)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer, parameter     :: buffer_size = 1000

    character(1000)        :: line

    integer                :: n, info
    real(DP), allocatable  :: x(:)
    real(DP), allocatable  :: data(:, :)

    ! ---

    if (ncol < 2) then
       RAISE_ERROR("Number of columns in input file is smaller than two.", ierror)
    endif

    allocate(x(buffer_size))
    allocate(data(buffer_size, ncol-1))

    x     = 0.0_DP
    data  = 0.0_DP

    n = 0
    read (un, '(A1000)', iostat=info)  line
    do while (info == 0 .and. len_trim(line) >= 3)
       n = n+1
       x(n)        = 0.0_DP
       data(n, :)  = 0.0_DP
       read (line, *, iostat=info)  x(n), data(n, :)

       if (info /= 0) then
          write (*, '(A,A)')  "line = ", trim(line)
          write (*, '(A,I5)') "ncol = ", ncol
          RAISE_ERROR("Number of columns given smaller than number of columns in input file?", ierror)
       endif

       if (present(xconv))  x(n)        = x(n)*xconv
       if (present(yconv))  data(n, :)  = data(n, :)*yconv

       read (un, '(A1000)', iostat=info)  line
    enddo

    if (n == 0) then
       RAISE_ERROR("No data found.", ierror)
    endif
    call nonuniform_spline_init(this, buffer_size, n, x(:), ncol-1, data(:, :))

    deallocate(x)
    deallocate(data)

  endsubroutine nonuniform_spline_read2


  !**********************************************************************
  ! Read with filename given
  !**********************************************************************
  subroutine nonuniform_spline_read2_fn(this, fn, ncol, xconv, yconv, ierror)
    implicit none

    type(spline_t), intent(out)       :: this
    character(*), intent(in)          :: fn
    integer, intent(in)               :: ncol
    real(DP), intent(in), optional    :: xconv
    real(DP), intent(in), optional    :: yconv(ncol-1)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: un

    ! ---

    un = fopen(fn, F_READ, ierror)
    PASS_ERROR(ierror)
    call nonuniform_spline_read2(this, un, ncol, xconv, yconv, ierror)
    PASS_ERROR(ierror)
    call fclose(un, ierror)
    PASS_ERROR(ierror)

  endsubroutine nonuniform_spline_read2_fn


  !**********************************************************************
  ! Output the spline table to a file
  !**********************************************************************
  subroutine nonuniform_spline_write(this, fn)
    implicit none

    type(spline_t), intent(in)  :: this
    character*(*), intent(in)   :: fn

    ! ---

    integer        :: un, i
    character(80)  :: fmt
    
    ! ---

    write (fmt, '(A,I2.2,A)')  "(", this%ncol+1, "ES20.10)"

    un = fopen(fn, F_WRITE)
    do i = 1, this%n
       write (un, trim(fmt))  this%x(i), this%y(i, :)
    enddo
    call fclose(un)

  endsubroutine nonuniform_spline_write


  !**********************************************************************
  ! Find a such that x(a) < x < x(a+1)
  !**********************************************************************
  integer function nonuniform_spline_interval(this, x, ierror)
    implicit none

    type(spline_t), intent(in)        :: this
    real(DP), intent(in)              :: x
    integer, intent(inout), optional  :: ierror

    ! ---

    integer        :: a, b, k

    ! ---

    !
    ! if x out of interpolation region
    !

    if (x < this%x(1) .or. x > this%x(this%n)) then
       a = 1
       b = 1
       nonuniform_spline_interval = 1
       RAISE_ERROR("x = "//x//" sits outside of the interval ["//this%x(1)//", "//this%x(this%n)//"] for which the spline is defined.", ierror)
    endif

    a = 1
    b = this%n

    do while (b-a > 1)
       k = (b+a)/2

       if (this%x(k) > x)then
          b = k
       else
          a = k
       endif
    enddo

    if (a == b) then
       RAISE_ERROR("a == b in interval(2)", ierror)
    endif

    nonuniform_spline_interval = a

  endfunction nonuniform_spline_interval


  !**********************************************************************
  ! Return the function value
  !**********************************************************************
  real(DP) function nonuniform_spline_f(this, col, x, i)
    implicit none

    type(spline_t), intent(in)  :: this
    integer, intent(in)         :: col
    real(DP), intent(in)        :: x
    integer, intent(in)         :: i         ! interval in which to find x

    ! ---

    real(DP)  :: dx, dx2, A, B

    ! ---


    dx        = this%x(i+1)-this%x(i)
    dx2       = dx**2

    A         = (this%x(i+1)-x)/dx
    B         = (x-this%x(i))/dx

    nonuniform_spline_f  =  A*this%y(i, col)+B*this%y(i+1, col)+((A**3-A)*this%d2y(i, col)+(B**3-B)*this%d2y(i+1, col))*dx2/6.

  endfunction nonuniform_spline_f


  !**********************************************************************
  ! Return the function value
  !**********************************************************************
  real(DP) function nonuniform_spline_f_unknown_interval(this, col, x)
    implicit none

    type(spline_t), intent(in)  :: this
    integer, intent(in)         :: col
    real(DP), intent(in)        :: x

    ! ---

    integer   :: i
    real(DP)  :: dx, dx2, A, B

    ! ---

    i         = interval(this, x)

    dx        = this%x(i+1)-this%x(i)
    dx2       = dx**2

    A         = (this%x(i+1)-x)/dx
    B         = (x-this%x(i))/dx

    nonuniform_spline_f_unknown_interval  =  A*this%y(i, col)+B*this%y(i+1, col)+((A**3-A)*this%d2y(i, col)+(B**3-B)*this%d2y(i+1, col))*dx2/6.

  endfunction nonuniform_spline_f_unknown_interval


  !**********************************************************************
  ! Return the derivative of the function
  !**********************************************************************
  real(DP) function nonuniform_spline_df(this, col, x, i)
    implicit none

    type(spline_t), intent(in)  :: this
    integer, intent(in)         :: col
    real(DP), intent(in)        :: x
    integer, intent(in)         :: i         ! interval in which to find x

    ! ---

    real(DP)  :: dx, dx2, A, B

    ! ---

    dx         = this%x(i+1)-this%x(i)
    dx2        = dx**2

    A          = (this%x(i+1)-x)/dx
    B          = (x-this%x(i))/dx

    nonuniform_spline_df  = (this%y(i+1, col)-this%y(i, col))/dx-(3*A**2-1)/6*dx*this%d2y(i, col)+(3*B**2-1)/6*dx*this%d2y(i+1, col)

  endfunction nonuniform_spline_df


  !**********************************************************************
  ! Return the function value and derivative
  !**********************************************************************
  subroutine nonuniform_spline_f_and_df(this, col, x, i, f, df)
    implicit none

    type(spline_t), intent(in)  :: this
    integer, intent(in)         :: col
    real(DP), intent(in)        :: x
    integer, intent(in)         :: i         ! interval in which to find x
    real(DP), intent(out)       :: f
    real(DP), intent(out)       :: df

    ! ---

    real(DP)  :: dx, dx2, A, B

    ! ---

    dx   = this%x(i+1)-this%x(i)
    dx2  = dx**2

    A    = (this%x(i+1)-x)/dx
    B    = (x-this%x(i))/dx

    f    =  A*this%y(i, col)+B*this%y(i+1, col)+((A**3-A)*this%d2y(i, col)+(B**3-B)*this%d2y(i+1, col))*dx2/6.
    df   = (this%y(i+1, col)-this%y(i, col))/dx-(3*A**2-1)/6*dx*this%d2y(i, col)+(3*B**2-1)/6*dx*this%d2y(i+1, col)

  endsubroutine nonuniform_spline_f_and_df

endmodule nonuniform_spline
