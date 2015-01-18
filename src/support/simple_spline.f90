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
!! Simple spline: Constant x-spacing, one column only
!<

#include "macros.inc"

module simple_spline
  use error_module
  use system_module
  use io

#ifndef SILENT
  use logging
#endif

  implicit none

  private

  public :: simple_spline_t
  type simple_spline_t

     integer            :: n = -1  ! number of entries
     real(DP)           :: x0
     real(DP)           :: cut     ! end point
     real(DP)           :: dx

     logical            :: associated = .false.

     real(DP), pointer  :: y(:)        => NULL()   ! value tables
     real(DP), pointer  :: d2y(:)      => NULL()   ! second derivatives table

     real(DP), pointer  :: coeff1(:)   => NULL()  ! spline coefficients
     real(DP), pointer  :: coeff2(:)   => NULL()  ! spline coefficients
     real(DP), pointer  :: coeff3(:)   => NULL()  ! spline coefficients
     real(DP), pointer  :: dcoeff1(:)  => NULL()  ! spline coefficients
     real(DP), pointer  :: dcoeff2(:)  => NULL()  ! spline coefficients
     real(DP), pointer  :: dcoeff3(:)  => NULL()  ! spline coefficients

  endtype simple_spline_t

  public :: init
  interface init
     module procedure simple_spline_init, simple_spline_init_from_func
  endinterface

  public :: del
  interface del
     module procedure simple_spline_del
  endinterface

  public :: associate
  interface associate
     module procedure simple_spline_associate
  endinterface

#ifndef SILENT
  public :: read
  interface read
     module procedure simple_spline_read
  endinterface

  public :: write
  interface write
     module procedure simple_spline_write
  endinterface
#endif

  public :: func
  interface func
     module procedure simple_spline_f, simple_spline_f_array
  endinterface

  public :: dfunc
  interface dfunc
     module procedure simple_spline_df
  endinterface

  public :: f_and_df
  interface f_and_df
     module procedure simple_spline_f_and_df
  endinterface

  public :: square_x_axis
  interface square_x_axis
     module procedure simple_spline_square_x_axis
  endinterface

  public :: scale_y_axis
  interface scale_y_axis
     module procedure simple_spline_scale_y_axis
  endinterface

  public :: log_memory_estimate
  interface log_memory_estimate
     module procedure simple_spline_log_memory_estimate
  endinterface

contains

  !>
  !! Initialize the Spline table
  !!
  !! Initialize the Spline table
  !<
  subroutine simple_spline_init(this, n, x0, dx, y)
    implicit none

    type(simple_spline_t), intent(out)  :: this
    integer, intent(in)                 :: n
    real(DP), intent(in)                :: x0
    real(DP), intent(in)                :: dx
    real(DP), intent(in)                :: y(n)

    ! ---

    real(DP), parameter  :: sig = 0.5

    integer  :: i, k
    real(DP) :: p, qn, un, u(n)

    ! ---

    this%associated  = .false.

    this%n     = n
    this%x0    = x0
    this%dx    = dx

    allocate(this%y(n))
    allocate(this%d2y(n))
    allocate(this%coeff1(n-1))
    allocate(this%coeff2(n-1))
    allocate(this%coeff3(n-1))
    allocate(this%dcoeff1(n-1))
    allocate(this%dcoeff2(n-1))
    allocate(this%dcoeff3(n-1))

    this%y(:)   = y(:)

    this%d2y(1) = 0.0_DP ! natural simple_spline, i.e. second derivatives vanishes
    u(1)        = 0.0_DP

    do i = 2, n-1
       p        = sig*this%d2y(i-1) + 2
       this%d2y(i) = (sig-1)/p   !-p/2
       u(i)     = (6.0_DP*((this%y(i+1)-this%y(i))/dx &
            -(this%y(i)-this%y(i-1))/dx)/(2*dx)-sig*u(i-1))/p
    enddo

    qn = 0.0_DP !natural simple_spline
    un = 0.0_DP
!    qn  = 0.5_DP
!    un  = 3.0_DP*(0.0_DP-(this%y(n)-this%y(n-1)))/(dx**2)

    this%d2y(n) = (un-qn*u(n-1))/(qn*this%d2y(n-1)+1.0_DP)

    do k = n-1, 1, -1
       this%d2y(k) = this%d2y(k)*this%d2y(k+1) + u(k)
    enddo

    do k = 1, n-1
       this%coeff1(k)   = this%y(k+1) - this%y(k) - ( 2*this%d2y(k) + this%d2y(k+1) )*this%dx**2/6
       this%coeff2(k)   = this%d2y(k) * this%dx**2/2
       this%coeff3(k)   = ( this%d2y(k+1) - this%d2y(k) )*this%dx**2/6
    enddo

    this%dcoeff1  = this%coeff1/this%dx
    this%dcoeff2  = 2*this%coeff2/this%dx
    this%dcoeff3  = 3*this%coeff3/this%dx

    this%cut = this%x0+dx*(n-1)

  endsubroutine simple_spline_init


  !>
  !! Initialize the from a function
  !!
  !! Initialize the from a function
  !<
  subroutine simple_spline_init_from_func(this, n, x0, cut, func, arg1, arg2, arg3, arg4, arg5, arg6)
    implicit none

    type(simple_spline_t), intent(out)  :: this
    integer, intent(in)                 :: n
    real(DP), intent(in)                :: x0
    real(DP), intent(in)                :: cut
    real(DP), external                  :: func
    real(DP), intent(in), optional      :: arg1
    real(DP), intent(in), optional      :: arg2
    real(DP), intent(in), optional      :: arg3
    real(DP), intent(in), optional      :: arg4
    real(DP), intent(in), optional      :: arg5
    real(DP), intent(in), optional      :: arg6

    ! ---

    real(DP)  :: y(n), dx
    integer   :: i

    ! ---

    dx = (cut-x0)/(n-1)

    ! FIXME! This is rather clumsy!

    narg1: if (present(arg1)) then

       narg2: if (present(arg2)) then

          narg3: if (present(arg3)) then

             narg4: if (present(arg4)) then

                narg5: if (present(arg5)) then

                   nargs6: if (present(arg6)) then

                      do i = 1, n
                         y(i) = func(x0+(i-1)*dx, arg1, arg2, arg3, arg4, arg5, arg6)
                      enddo

                   else

                      do i = 1, n
                         y(i) = func(x0+(i-1)*dx, arg1, arg2, arg3, arg4, arg5)
                      enddo

                   endif nargs6

                else

                   do i = 1, n
                      y(i) = func(x0+(i-1)*dx, arg1, arg2, arg3, arg4)
                   enddo

                endif narg5

             else

                do i = 1, n
                   y(i) = func(x0+(i-1)*dx, arg1, arg2, arg3)
                enddo

             endif narg4

          else
             
             do i = 1, n
                y(i) = func(x0+(i-1)*dx, arg1, arg2)
             enddo

          endif narg3

       else

          do i = 1, n
             y(i) = func(x0+(i-1)*dx, arg1)
          enddo
          
       endif narg2

    else

       do i = 1, n
          y(i) = func(x0+(i-1)*dx)
       enddo
       
    endif narg1

    call simple_spline_init(this, n, x0, dx, y)

  endsubroutine simple_spline_init_from_func


  !>
  !! Delete a Spline table
  !!
  !! Delete a Spline table
  !<
  elemental subroutine simple_spline_del(this)
    implicit none

    type(simple_spline_t), intent(inout)  :: this

    ! ---

    if (.not. this%associated) then
       if (associated(this%y))        deallocate(this%y)
       if (associated(this%d2y))      deallocate(this%d2y)
       if (associated(this%coeff1))   deallocate(this%coeff1)
       if (associated(this%coeff2))   deallocate(this%coeff2)
       if (associated(this%coeff3))   deallocate(this%coeff3)
       if (associated(this%dcoeff1))  deallocate(this%dcoeff1)
       if (associated(this%dcoeff2))  deallocate(this%dcoeff2)
       if (associated(this%dcoeff3))  deallocate(this%dcoeff3)

       this%y       => NULL()
       this%d2y     => NULL()
       this%coeff1  => NULL()
       this%coeff2  => NULL()
       this%coeff3  => NULL()
       this%dcoeff1 => NULL()
       this%dcoeff2 => NULL()
       this%dcoeff3 => NULL()
    endif

  endsubroutine simple_spline_del


  !>
  !! Associate this Spline to another one
  !!
  !! This will lead to a shallow copy of the spline, i.e. not the data
  !! is copied but simply pointers to the other spline's data.
  !<
  subroutine simple_spline_associate(this, w)
    implicit none

    type(simple_spline_t), intent(inout) :: this
    type(simple_spline_t), intent(in)    :: w

    ! --

    this%associated = .true.

    this%n        = w%n
    this%x0       = w%x0
    this%cut      = w%cut
    this%dx       = w%dx

    this%y        => w%y
    this%d2y      => w%d2y

    this%coeff1   => w%coeff1
    this%coeff2   => w%coeff2
    this%coeff3   => w%coeff3

    this%dcoeff1  => w%dcoeff1
    this%dcoeff2  => w%dcoeff2
    this%dcoeff3  => w%dcoeff3

  endsubroutine simple_spline_associate


  !>
  !! Return the function value
  !!
  !! Return the function value
  !<
  function simple_spline_f(this, x, ierror)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(DP), intent(in)               :: x
    real(DP)                           :: simple_spline_f
    integer, intent(inout), optional   :: ierror

    ! ---

#if 0
    real(DP)  :: dx, dx2, A, B, xf
    integer   :: i

    ! ---

    dx        = this%dx
    dx2       = dx**2

    if (x == this%cut) then
       xf        = this%n
       i         = floor(this%n-1)
    else
       xf        = (x-this%x0)/dx+1
       i         = floor(xf)
    endif

    if (i < 1 .or. i >= this%n) then
       simple_spline_f  = 1.0_DP
       RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    A         = (i+1)-xf
    B         = xf-i

    simple_spline_f  =  A*this%y(i)+B*this%y(i+1)+((A**3-A)*this%d2y(i)+(B**3-B)*this%d2y(i+1))*dx2/6.
#endif

    real(DP)  :: dx, B, xf
    integer   :: i

    ! ---

    dx  = this%dx

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = floor(xf)
    endif

    if (i < 1 .or. i >= this%n) then
       simple_spline_f  = 1.0_DP
       RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    B                = xf-i
    simple_spline_f  = this%y(i) + B*(this%coeff1(i) + B*(this%coeff2(i) + B*this%coeff3(i)))

  endfunction simple_spline_f


  !>
  !! Return the function value for an array of arguments
  !!
  !! Return the function value for an array of arguments
  !<
  function simple_spline_f_array(this, x)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(DP), intent(in)               :: x(:)
    real(DP)                           :: simple_spline_f_array(lbound(x, 1):ubound(x, 1))

    ! ---

    integer  :: i

    ! ---

    do i = lbound(x, 1), ubound(x, 1)
       simple_spline_f_array(i)  = simple_spline_f(this, x(i))
    enddo

  endfunction simple_spline_f_array


  !>
  !! Return the derivative of the function
  !!
  !! Return the derivative of the function
  !<
  function simple_spline_df(this, x, ierror)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(DP), intent(in)               :: x
    real(DP)                           :: simple_spline_df
    integer, intent(inout), optional   :: ierror

    ! ---

#if 0
    real(DP)  :: dx, dx2, A, B, xf
    integer   :: i

    ! ---

    dx        = this%dx
    dx2       = dx**2

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = xf
    endif

    if (i < 1 .or. i >= this%n) then
       simple_spline_df  = 1.0_DP
       RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    A         = (i+1)-xf
    B         = xf-i

    simple_spline_df  = (this%y(i+1)-this%y(i))/dx-(3*A**2-1)/6*dx*this%d2y(i)+(3*B**2-1)/6*dx*this%d2y(i+1)
#endif

    real(DP)  :: dx, B, xf
    integer   :: i

    ! ---

    dx  = this%dx

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = floor(xf)
    endif

    if (i < 1 .or. i >= this%n) then
       simple_spline_df  = 1.0_DP
       RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    B                 = xf-i
    simple_spline_df  = this%dcoeff1(i) + B*(this%dcoeff2(i) + B*this%dcoeff3(i))

  endfunction simple_spline_df


  !>
  !! Return the function value and derivative
  !!
  !! Return the function value and derivative
  !<
  subroutine simple_spline_f_and_df(this, x, f, df, ierror)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    real(DP), intent(in)               :: x
    real(DP), intent(out)              :: f
    real(DP), intent(out)              :: df
    integer, intent(inout), optional   :: ierror

    ! ---

#if 0
    real(DP)  :: dx, dx2, A, B, xf
    integer   :: i

    ! ---

    dx   = this%dx
    dx2  = dx**2

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = xf
    endif

    if (i < 1 .or. i >= this%n) then
       f   = 1.0_DP
       df  = 1.0_DP
       RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    A   = (i+1)-xf
    B   = xf-i

    f   =  A*this%y(i) + B*this%y(i+1) + ((A**3-A)*this%d2y(i) + (B**3-B)*this%d2y(i+1))*dx2/6.
    df  = (this%y(i+1)-this%y(i))/dx-(3*A**2-1)/6*dx*this%d2y(i)+(3*B**2-1)/6*dx*this%d2y(i+1)
#endif

    real(DP)  :: dx, B, xf
    integer   :: i

    ! ---

    dx  = this%dx

    if (x == this%cut) then
       xf        = this%n
       i         = this%n-1
    else
       xf        = (x-this%x0)/dx+1
       i         = floor(xf)
    endif

    if (i < 1 .or. i >= this%n) then
       f   = 1.0_DP
       df  = 1.0_DP
       RAISE_ERROR("x = " // x // " outside of the defined interval.", ierror)
    endif

    B   = xf-i
    f   = this%y(i) + B*(this%coeff1(i) + B*(this%coeff2(i) + B*this%coeff3(i)))
    df  = this%dcoeff1(i) + B*(this%dcoeff2(i) + B*this%dcoeff3(i))

  endsubroutine simple_spline_f_and_df


#ifndef SILENT
  !>
  !! Read the spline table from a file
  !!
  !! Read the spline table from a file
  !<
  subroutine simple_spline_read(this, f, n, x0, dx, pad)
    implicit none

    type(simple_spline_t), intent(inout)  :: this
    integer,               intent(in)     :: f
    integer,               intent(in)     :: n
    real(DP),              intent(in)     :: x0
    real(DP),              intent(in)     :: dx
    real(DP),    optional, intent(in)     :: pad(:)

    ! ---

    real(DP)  :: y(n)
    real(DP), allocatable  :: ytmp(:)

    ! ---

    call read_ascii(f, y)
    if (present(pad)) then
       allocate(ytmp(n+size(pad)))
       ytmp(1:n)   = y
       ytmp(n+1:)  = pad
       call init(this, n+size(pad), x0, dx, ytmp)
       deallocate(ytmp)
    else
       call init(this, n, x0, dx, y)
    endif

  endsubroutine simple_spline_read
  


  !>
  !! Output the spline table to a file
  !!
  !! Output the spline table to a file. The file has three colums containing
  !! the x-value, the function value and its derivative. The optional parameter *dx*
  !! specifies the spacing at which the data is writen to the file. If it
  !! is not specified, the natural spacing (i.e., the one that is used for internal
  !! data storage and that has been specified in the init method) is used.
  !<
  subroutine simple_spline_write(this, fn, dx)
    implicit none

    type(simple_spline_t), intent(in)  :: this
    character*(*), intent(in)          :: fn
    real(DP), intent(in), optional     :: dx

    ! ---

    integer       :: un
    character(80) :: fmt
    real(DP)      :: x, f, df
    
    ! ---

    write (fmt, '(A,I2.2,A)')  "(3ES20.10)"

    un = fopen(fn, mode=F_WRITE)

    x  = this%x0
    do while (x <= this%cut)
       call simple_spline_f_and_df(this, x, f, df)
       write (un, trim(fmt))  x, f, df

       if (present(dx)) then
          x  = x + dx
       else
          x  = x + this%dx
       endif
    enddo

    call fclose(un)

  endsubroutine simple_spline_write
#endif


  !>
  !! Compute the square of the x-axis
  !!
  !! Rescales the x-axis nonlinearly to contain its square. This
  !! is intended to reduce the computation of sqrts in distances.
  !<
  subroutine simple_spline_square_x_axis(this, n, ierror)
    implicit none

    type(simple_spline_t),           intent(inout)  :: this
    integer,               optional, intent(in)     :: n
    integer,               optional, intent(inout)  :: ierror

    ! ---

    integer   :: i, my_n
    real(DP)  :: x0, dx

    real(DP), allocatable  :: y(:)

    ! ---

    x0  = this%x0**2
    dx  = (this%x0+this%dx)**2 - x0

    if (present(n)) then
       my_n = n
    else
       my_n = int( ( this%cut**2 - x0 )/dx )+1
    endif

    dx  = (this%cut**2 - x0)/(my_n-1)

!    call print("simple_spline_square_x_axis : n = " // my_n // ", dx = " // dx, PRINT_VERBOSE)

    allocate(y(my_n))

    do i = 1, my_n
       y(i) = func(this, sqrt(x0 + (i-1)*dx), ierror=ierror)
       PASS_ERROR(ierror)
    enddo

    call del(this)
    call init(this, my_n, x0, dx, y)

    deallocate(y)

  endsubroutine simple_spline_square_x_axis


  !>
  !! Scale the y-axis
  !!
  !! Rescales the y-axis by a value
  !<
  subroutine simple_spline_scale_y_axis(this, fac)
    implicit none

    type(simple_spline_t), intent(inout)  :: this
    real(DP),              intent(in)     :: fac

    ! ---

    this%y        = fac*this%y
    this%d2y      = fac*this%d2y
    this%coeff1   = fac*this%coeff1
    this%coeff2   = fac*this%coeff2
    this%coeff3   = fac*this%coeff3
    this%dcoeff1  = fac*this%dcoeff1
    this%dcoeff2  = fac*this%dcoeff2
    this%dcoeff3  = fac*this%dcoeff3

  endsubroutine simple_spline_scale_y_axis


  !>
  !! Log the memory require
  !!
  !! Log the memory requirement of the Spline object
  !<
  subroutine simple_spline_log_memory_estimate(this)
    implicit none

    type(simple_spline_t), intent(in)  :: this

    ! ---

#ifndef SILENT
    call log_memory_estimate(this%y, "y")
    call log_memory_estimate(this%d2y, "d2y")
#endif

  endsubroutine simple_spline_log_memory_estimate

endmodule
