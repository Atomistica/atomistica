!! ======================================================================
!! MDCORE - Interatomic potential library
!! https://github.com/pastewka/mdcore
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
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

     integer            :: nx = 1
     integer            :: ny = 1

     integer            :: nboxs

     real(DP), pointer  :: coeff(:, :, :)  => NULL()

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
  subroutine table2d_init(t, nx, ny, values, ierror)
    implicit none

    type(table2d_t), intent(inout)    :: t
    integer, intent(in)               :: nx
    integer, intent(in)               :: ny
    real(DP), intent(in)              :: values(0:nx, 0:ny)
    integer, intent(inout), optional  :: ierror

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

    t%nx     = nx
    t%ny     = ny
    t%nboxs  = nx*ny

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

    do nhbox = 0, nx-1
       do ncbox = 0, ny-1

          icol = 1+ny*nhbox+ncbox

          do icorn = 1, ncorn

             irow = icorn
             nx1  = ix1(icorn)+nhbox
             nx2  = ix2(icorn)+ncbox
             !   values of function and derivatives at corner.
             B(irow   ,icol) = values(nx1, nx2)
             B(irow+4 ,icol) = 0.0
             B(irow+8 ,icol) = 0.0
             B(irow+12,icol) = 0.0

          enddo
       enddo
    enddo

    !
    ! solve by gauss-jordan elimination with full pivoting.
    !

!    call gaussjhc(a,n,np,b,m,mp)
    call dgesv(npara, t%nboxs, A, npara, ipiv, B, npara, info)

    if (info /= 0) then
       RAISE_ERROR("dgesv failed.", ierror)
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

    deallocate(t%coeff)

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
    real(DP), intent(in)         :: nhi
    real(DP), intent(in)         :: nci

    real(DP), intent(out)        :: hch
    real(DP), intent(out)        :: dhchdh
    real(DP), intent(out)        :: dhchdc

    ! ---

    real(DP)  :: x1, x2, coefij
    real(DP)  :: shch, shchdc

    integer   :: nhbox, ncbox, i, j, ibox

    ! ---

!    write (*, *)  nhi, nci

    nhbox = int( nhi )
    ncbox = int( nci )

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

       hch    = 0.0
       dhchdh = 0.0
       dhchdc = 0.0
       do i = 4, 1, -1
          shch   = 0.0
          shchdc = 0.0
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

    integer               :: i, j
    real(DP)              :: row(0:this%nx-1), dummy1, dummy2
    character(1000)       :: fmt

    ! ---

    if (present(indent)) then
       fmt = "(" // (indent+5) // "X," // (this%nx+1) // "I20)"
    else
       fmt = "(5X," // (this%nx+1) // "I20)"
    endif

    write (un, fmt)  (/ ( i, i=0, this%nx-1 ) /)

    if (present(indent)) then
       fmt = "(" // indent // "X,I3,' -'," // (this%nx+1) // "ES20.10)"
    else
       fmt = "(4I,1X," // (this%nx+1) // "ES20.10)"
    endif

    do i = 0, this%ny-1
       do j = 0, this%nx-1
          call eval(this, i*1.0_DP, j*1.0_DP, row(j), dummy1, dummy2)
       enddo

       write (un, fmt)  i, row
    enddo

  endsubroutine table2d_print_un

endmodule table2d
