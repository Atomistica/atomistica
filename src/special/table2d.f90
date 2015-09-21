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
!! 2D cubic spline interpolation
!<

#include "macros.inc"

module table2d
  use libAtoms_module

  use logging, only: ilog

  implicit none

  private

  public :: table2d_t
  type table2d_t

     integer               :: nx = 1
     integer               :: ny = 1 

     integer               :: nboxs

     real(DP), allocatable :: coeff(:, :, :)

  endtype table2d_t

  integer, parameter, private  :: npara = 4*4   ! 4^dim
  integer, parameter, private  :: ncorn = 4

  public :: init
  interface init
     module procedure table2d_init
  endinterface

  public :: del
  interface del
     module procedure table2d_del
  endinterface

  public :: eval
  interface eval
     module procedure table2d_eval
  endinterface

!  interface print
!     module procedure table2d_print, table2d_print_un
!  endinterface

!  interface prlog
!     module procedure table2d_prlog
!  endinterface

  public :: table2d_prlog

contains

  !>
  !! generates the coefficients for bicubic interpolation of fch(ni,nj)
  !! copyright: Keith Beardmore 30/11/93.
  !!            Lars Pastewka 05/07
  !<
  subroutine table2d_init(t, nx, ny, values, dvdx, dvdy, error)
    implicit none

    type(table2d_t),    intent(inout) :: t
    integer,            intent(in)    :: nx
    integer,            intent(in)    :: ny
    real(DP),           intent(in)    :: values(0:, 0:)
    real(DP), optional, intent(in)    :: dvdx(0:nx, 0:ny), dvdy(0:nx, 0:ny)
    integer,  optional, intent(inout) :: error

    ! ---

    !
    ! calculate 2-d cubic parameters within each box.
    !
    ! normalised coordinates.
    !  4--<--3
    !  |     ^
    !  v     |
    !  1-->--2
    !

    integer, parameter     :: ix1(ncorn) = (/ 0,1,1,0 /)
    integer, parameter     :: ix2(ncorn) = (/ 0,0,1,1 /)

    real(DP)               :: A(npara, npara)
    real(DP), allocatable  :: B(:, :)
    integer                :: ipiv(npara)

    integer                :: icorn, irow, icol, ibox, nx1, nx2
    integer                :: npow1, npow2, npow1m, npow2m
    integer                :: i, j, nhbox, ncbox, info

    ! ---

    ! Bounds checking

    if (lbound(values, 1) /= 0 .or. ubound(values, 1) /= nx) then
       RAISE_ERROR("First index of *values* must run from 0 to " // nx // ", but does run from " // lbound(values, 1) // " to " // ubound(values, 1) // ".", error)
    endif
    if (lbound(values, 2) /= 0 .or. ubound(values, 2) /= ny) then
       RAISE_ERROR("Second index of *values* must run from 0 to " // ny // ", but does run from " // lbound(values, 2) // " to " // ubound(values, 2) // ".", error)
    endif

    t%nx     = nx
    t%ny     = ny
    t%nboxs  = nx*ny

    if (allocated(t%coeff))  deallocate(t%coeff)
    allocate(t%coeff(t%nboxs, 4, 4))
    allocate(B(npara, t%nboxs))

    !
    ! for each box, create and solve the matrix equatoion.
    !    / values of  \     /              \     / function and \
    !  a |  products  | * x | coefficients | = b |  derivative  |
    !    \within cubic/     \ of 2d cubic  /     \    values    /
    !

    !
    ! construct the matrix.
    ! this is the same for all boxes as coordinates are normalised.
    ! loop through corners.
    !

    do icorn = 1, ncorn
       irow  = icorn
       nx1   = ix1(icorn)
       nx2   = ix2(icorn)
       ! loop through powers of variables.
       do npow1 = 0, 3
          do npow2 = 0, 3

             npow1m = npow1-1
             if (npow1m < 0)  npow1m = 0
             npow2m = npow2-1
             if (npow2m < 0)  npow2m=0

             icol = 1+4*npow1+npow2

             !   values of products within cubic and derivatives.
             A(irow   ,icol) = 1.0_DP*(      nx1**npow1         *nx2**npow2 )
             A(irow+4 ,icol) = 1.0_DP*(npow1*nx1**npow1m        *nx2**npow2 )
             A(irow+8 ,icol) = 1.0_DP*(      nx1**npow1  * npow2*nx2**npow2m)
             A(irow+12,icol) = 1.0_DP*(npow1*nx1**npow1m * npow2*nx2**npow2m)

          enddo
       enddo
    enddo

    !
    ! construct the 16 r.h.s. vectors ( 1 for each box ).
    ! loop through boxes.
    !

    B = 0.0_DP
    do nhbox = 0, nx-1
       do ncbox = 0, ny-1

          icol = 1+ny*nhbox+ncbox

          do icorn = 1, ncorn

             irow = icorn
             nx1  = ix1(icorn)+nhbox
             nx2  = ix2(icorn)+ncbox
             !   values of function and derivatives at corner.
             B(irow   ,icol) = values(nx1, nx2)
             !   all derivatives are supposed to be zero
             if (present(dvdx)) then
                B(irow+ ncorn  ,icol) = dvdx(nx1, nx2)
             endif
             if (present(dvdy)) then
                B(irow+ 2*ncorn,icol) = dvdy(nx1, nx2)
             endif
          enddo
       enddo
    enddo

    !
    ! solve by gauss-jordan elimination with full pivoting.
    !

!    call gaussjhc(a,n,np,b,m,mp)
    call dgesv(npara, t%nboxs, A, npara, ipiv, B, npara, info)

    if (info /= 0) then
       RAISE_ERROR("dgesv failed.", error)
    endif

    !
    ! get the coefficient values.
    !

    do ibox = 1, t%nboxs
       icol = ibox
       do i = 1, 4
          do j = 1, 4
             irow = 4*(i-1)+j

             t%coeff(ibox, i, j) = B(irow, icol)
          enddo
       enddo
    enddo

    deallocate(B)

  endsubroutine table2d_init


  !>
  !! Free memory allocated for the spline coefficients
  !<
  elemental subroutine table2d_del(t)
    implicit none

    type(table2d_t), intent(inout)  :: t

    ! ---

    if (allocated(t%coeff)) then
       deallocate(t%coeff)
    endif

  endsubroutine table2d_del


  !>
  !! Compute function values and derivatives
  !!
  !! bicubic interpolation of hch.
  !! assumes 0.0 <= nhi,nci < 4.0
  !! copyright: Keith Beardmore 30/11/93.
  !!            Lars Pastewka 05/07
  !!
  !<
  subroutine table2d_eval(t, nhi, nci, hch, dhchdh, dhchdc)
    implicit none

    type(table2d_t), intent(in)  :: t
    real(DP),        intent(in)  :: nhi
    real(DP),        intent(in)  :: nci

    real(DP),        intent(out) :: hch
    real(DP),        intent(out) :: dhchdh
    real(DP),        intent(out) :: dhchdc

    ! ---

    real(DP)  :: x1, x2, coefij
    real(DP)  :: shch, shchdc

    integer   :: nhbox, ncbox, i, j, ibox

    ! ---

!    write (*, *)  nhi, nci

    nhbox = int( nhi )
    if (nhbox < 0)      nhbox = 0
    if (nhbox >= t%nx)  nhbox = t%nx-1
    ncbox = int( nci )
    if (ncbox < 0)      ncbox = 0
    if (ncbox >= t%ny)  ncbox = t%ny-1

    !
    !   find which box we're in and convert to normalised coordinates.
    !

    ibox = 1+t%ny*nhbox+ncbox
    x1   = nhi - nhbox
    x2   = nci - ncbox

!!$    if (x1 == 0.0 .and. x2 == 0.0) then
!!$
!!$       hch    = t%coeff(ibox, 1, 1)
!!$       dhchdh = 0.0_DP
!!$       dhchdc = 0.0_DP
!!$
!!$    else

       hch    = 0.0_DP
       dhchdh = 0.0_DP
       dhchdc = 0.0_DP
       do i = 4, 1, -1
          shch   = 0.0_DP
          shchdc = 0.0_DP
          do j = 4, 1, -1
                         coefij = t%coeff(ibox, i, j)
                         shch   =   shch*x2 +       coefij
             if (j > 1)  shchdc = shchdc*x2 + (j-1)*coefij
          enddo
                      hch    = hch   *x1 +       shch
          if (i > 1)  dhchdh = dhchdh*x1 + (i-1)*shch
                      dhchdc = dhchdc*x1 +       shchdc
       enddo

!!$    endif

  endsubroutine table2d_eval


  !>
  !! Print to screen
  !!
  !! Print to screen
  !<
  subroutine table2d_print(this, indent)
    implicit none

    type(table2d_t), intent(in)    :: this
    integer, intent(in), optional  :: indent

    ! ---

    call table2d_print_un(7, this)

  endsubroutine table2d_print


  !>
  !! Print to log file
  !!
  !! Print to log file
  !<
  subroutine table2d_prlog(this, indent)
    implicit none

    type(table2d_t), intent(in)    :: this
    integer, intent(in), optional  :: indent

    ! ---

    call table2d_print_un(ilog, this, indent)

  endsubroutine table2d_prlog


  !>
  !! Print to unit
  !!
  !! Print to unit
  !<
  subroutine table2d_print_un(un, this, indent)
    implicit none

    integer, intent(in)            :: un
    type(table2d_t), intent(in)    :: this
    integer, intent(in), optional  :: indent
    
    ! ---

    integer               :: i, j, k
    real(DP)              :: row(0:this%nx), val, dummy1, dummy2
    character(1000)       :: fmtstart, fmt

    ! ---

    if (present(indent)) then
       fmt = "(" // (indent+5) // "X," // (this%nx+1) // "I20)"
    else
       fmt = "(5X," // (this%nx+1) // "I20)"
    endif

    write (un, fmt)  (/ ( i, i=0, this%nx ) /)

    if (present(indent)) then
       fmtstart = "(" // indent // "X,I3,' -'"
    else
       fmtstart = "(4I,1X"
    endif

    do j = 0, this%ny
       fmt = fmtstart
       k = 0
       do i = 0, this%nx
          call eval(this, i*1.0_DP, j*1.0_DP, val, dummy1, dummy2)
          if (abs(val) > 1e-12) then
             fmt = trim(fmt) // ",ES20.10"
             row(k) = val
             k = k+1
          else
             fmt = trim(fmt) // ",'       ----------   '"
          endif
       enddo
       fmt = trim(fmt) // ")"

       write (un, fmt)  j, row(0:k-1)
    enddo

  endsubroutine table2d_print_un

endmodule table2d
