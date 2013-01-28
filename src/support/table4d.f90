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
!! 4D cubic spline interpolation
!<

#include "macros.inc"

module table4d
  use libAtoms_module

  use logging, only: ilog

  private

  public :: table4d_t
  type table4d_t

     integer            :: nx = 1
     integer            :: ny = 1
     integer            :: nz = 1
     integer            :: nt = 1

     integer            :: nboxs

     real(DP), pointer  :: coeff(:, :, :, :, :)  => NULL()

  endtype table4d_t

  integer, parameter, private  :: npara = 4*4*4*4   ! 4^dim
  integer, parameter, private  :: ncorn = 16        ! 2^dim

  public :: init
  interface init
     module procedure table4d_init
  endinterface

  public :: del
  interface del
     module procedure table4d_del
  endinterface

  public :: f
  interface f
     module procedure table4d_f
  endinterface

contains

  !>
  !! generates the coefficients for bicubic interpolation of fch(ni,nj)
  !! copyright: Keith Beardmore 30/11/93.
  !!            Lars Pastewka 05/07
  !<
  subroutine table4d_init(t, nx, ny, nz, nt, values, ierror)
    implicit none

    type(table4d_t), intent(inout)    :: t
    integer, intent(in)               :: nx
    integer, intent(in)               :: ny
    integer, intent(in)               :: nz
    integer, intent(in)               :: nt
    real(DP), intent(in)              :: values(0:nx, 0:ny, 0:nz, 0:nt)
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

    integer, parameter     :: ix1(ncorn) = (/ 0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0 /)
    integer, parameter     :: ix2(ncorn) = (/ 0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1 /)
    integer, parameter     :: ix3(ncorn) = (/ 0,0,0,0,1,1,1,1,0,0,0,0,1,1,1,1 /)
    integer, parameter     :: ix4(ncorn) = (/ 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 /)

    real(DP)               :: A(npara, npara)
    real(DP), allocatable  :: B(:, :)
    integer                :: ipiv(npara)

    integer                :: icorn, irow, icol, ibox, nx1, nx2, nx3, nx4
    integer                :: npow1, npow2, npow3, npow4, npow1m, npow2m, npow3m, npow4m
    integer                :: i, j, k, l, nibox, njbox, ncbox, ndbox, info

    ! ---

    t%nx     = nx
    t%ny     = ny
    t%nz     = nz
    t%nt     = nt
    t%nboxs  = nx*ny*nz*nt

    allocate(t%coeff(t%nboxs, 4, 4, 4, 4))
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
       irow = icorn
       nx1  = ix1(icorn)
       nx2  = ix2(icorn)
       nx3  = ix3(icorn)
       nx4  = ix4(icorn)
       ! loop through powers of variables.
       do npow1 = 0, 3
          do npow2 = 0, 3
             do npow3 = 0, 3
                do npow4 = 0, 3
                   npow1m = npow1-1
                   if (npow1m < 0)  npow1m=0
                   npow2m = npow2-1
                   if (npow2m < 0)  npow2m=0
                   npow3m = npow3-1
                   if (npow3m < 0)  npow3m=0
                   npow4m = npow4-1
                   if (npow4m < 0)  npow4m=0
                   icol = 1+4*4*4*npow1+4*4*npow2+4*npow3+npow4
                   ! values of products within cubic and derivatives.
                   A(irow   ,icol)=1.0*         &
                          (       nx1**npow1    &
                                 *nx2**npow2    &
                                 *nx3**npow3    &
                                 *nx4**npow4    )
                   A(irow+ncorn ,icol)=1.0*     &
                          ( npow1*nx1**npow1m   &
                                 *nx2**npow2    &
                                 *nx3**npow3    &
                                 *nx4**npow4    )
                   A(irow+2*ncorn,icol)=1.0*    &
                          (       nx1**npow1    &
                           *npow2*nx2**npow2m   &
                                 *nx3**npow3    &
                                 *nx4**npow4    )
                   A(irow+3*ncorn,icol)=1.0*    &
                          (       nx1**npow1    &
                                 *nx2**npow2    &
                           *npow3*nx3**npow3m   &
                                 *nx4**npow4    )     
                   A(irow+4*ncorn,icol)=1.0*    &
                        (       nx1**npow1      &
                               *nx2**npow2      &
                               *nx3**npow3      &
                         *npow4*nx4**npow4m     )
                   A(irow+5*ncorn,icol)=1.0*    &
                        ( npow1*nx1**npow1m     &
                         *npow2*nx2**npow2m     &
                               *nx3**npow3      &
                               *nx4**npow4      )
                   A(irow+6*ncorn,icol)=1.0*    &
                        ( npow1*nx1**npow1m     &
                               *nx2**npow2      &
                         *npow3*nx3**npow3m     &
                               *nx4**npow4      )
                   A(irow+7*ncorn,icol)=1.0*    &
                        ( npow1*nx1**npow1m     &
                               *nx2**npow2      &
                               *nx3**npow3      &
                         *npow4*nx4**npow4m     )
                   A(irow+8*ncorn,icol)=1.0*    &
                        (       nx1**npow1      &
                         *npow2*nx2**npow2m     &
                         *npow3*nx3**npow3m     &
                               *nx4**npow4      )
                   A(irow+9*ncorn,icol)=1.0*    &
                        (       nx1**npow1      &
                         *npow2*nx2**npow2m     &
                               *nx3**npow3      &
                         *npow4*nx4**npow4m     )
                   A(irow+10*ncorn,icol)=1.0*   &
                        (       nx1**npow1      &
                               *nx2**npow2      &
                         *npow3*nx3**npow3m     &
                         *npow4*nx4**npow4m     )
                   A(irow+11*ncorn,icol)=1.0*   &
                        ( npow1*nx1**npow1m     &
                         *npow2*nx2**npow2m     &
                         *npow3*nx3**npow3m     &
                               *nx4**npow4      )
                   A(irow+12*ncorn,icol)=1.0*   &
                        ( npow1*nx1**npow1m     &
                         *npow2*nx2**npow2m     &
                               *nx3**npow3      &
                         *npow4*nx4**npow4m     )
                   A(irow+13*ncorn,icol)=1.0*   &
                        ( npow1*nx1**npow1m     &
                               *nx2**npow2      &
                         *npow3*nx3**npow3m     &
                         *npow4*nx4**npow4m     )
                   A(irow+14*ncorn,icol)=1.0*   &
                        (       nx1**npow1      &
                         *npow2*nx2**npow2m     &
                         *npow3*nx3**npow3m     &
                         *npow4*nx4**npow4m     )     
                   A(irow+15*ncorn,icol)=1.0*   &
                        ( npow1*nx1**npow1m     &
                         *npow2*nx2**npow2m     &
                         *npow3*nx3**npow3m     &
                         *npow4*nx4**npow4m     )     
                enddo
             enddo
          enddo
       enddo
    enddo

    !
    ! construct the 16 r.h.s. vectors ( 1 for each box ).
    ! loop through boxes.
    !

    do nibox = 0, nx-1
       do njbox = 0, ny-1
          do ncbox = 0, nz-1
             do ndbox = 0, nt-1
                icol = 1+t%nx*(t%ny*(t%nz*ndbox+ncbox)+njbox)+nibox
                do icorn = 1, ncorn
                   irow = icorn
                   nx1  = ix1(icorn)+nibox
                   nx2  = ix2(icorn)+njbox
                   nx3  = ix3(icorn)+ncbox
                   nx4  = ix4(icorn)+ndbox
                   ! values of function and derivatives at corner.
                   B(irow         ,icol) = values(nx1, nx2, nx3, nx4)
                   !   all derivatives are supposed to be zero
                   B(irow+ ncorn  ,icol) = 0.0
                   B(irow+ 2*ncorn,icol) = 0.0
                   B(irow+ 3*ncorn,icol) = 0.0
                   B(irow+ 4*ncorn,icol) = 0.0
                   B(irow+ 5*ncorn,icol) = 0.0
                   B(irow+ 6*ncorn,icol) = 0.0
                   B(irow+ 7*ncorn,icol) = 0.0
                   B(irow+ 8*ncorn,icol) = 0.0
                   B(irow+ 9*ncorn,icol) = 0.0
                   B(irow+10*ncorn,icol) = 0.0
                   B(irow+11*ncorn,icol) = 0.0
                   B(irow+12*ncorn,icol) = 0.0
                   B(irow+13*ncorn,icol) = 0.0
                   B(irow+14*ncorn,icol) = 0.0
                   B(irow+15*ncorn,icol) = 0.0
                enddo
             enddo
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
       do i= 1, 4
          do j = 1, 4
             do k = 1, 4
                do l = 1, 4
                   irow=4*4*4*(i-1)+4*4*(j-1)+4*(k-1)+l
                   t%coeff(ibox,i,j,k,l) = b(irow,icol)
                enddo
             enddo
          enddo
       enddo
    enddo

    deallocate(B)

  endsubroutine table4d_init


  !>
  !! Free memory allocated for the spline coefficients
  !<
  elemental subroutine table4d_del(t)
    implicit none

    type(table4d_t), intent(inout)  :: t

    ! ---

    deallocate(t%coeff)

  endsubroutine table4d_del


  !>
  !! Compute function values and derivatives
  !!
  !! bicubic interpolation of hch.
  !! assumes 0.0 <= nhi,nci < 4.0
  !! copyright: Keith Beardmore 30/11/93.
  !!            Lars Pastewka 05/07
  !!
  !<
  subroutine table4d_f(t, nti, ntj, nconji, nconjj, fcc, dfccdi, dfccdj, &
       dfccdc, dfccdd)
    implicit none

    type(table4d_t), intent(in)  :: t
    real(DP), intent(in)         :: nti
    real(DP), intent(in)         :: ntj
    real(DP), intent(in)         :: nconji
    real(DP), intent(in)         :: nconjj
    real(DP), intent(out)        :: fcc
    real(DP), intent(out)        :: dfccdi
    real(DP), intent(out)        :: dfccdj
    real(DP), intent(out)        :: dfccdc
    real(DP), intent(out)        :: dfccdd

    ! ---

    integer   :: nibox, njbox, ncbox, ndbox, ibox, i, j, k, l
    real(DP)  :: x1, x2, x3, x4
    real(DP)  :: sfcc, sfccdj, sfccdc, sfccdd
    real(DP)  :: tfcc, tfccdc, tfccdd
    real(DP)  :: ufcc, ufccdd
    real(DP)  :: coefij

    !
    !   find which box we're in and convert to normalised coordinates.
    !

    nibox = int( nti )
    x1    = nti - nibox
    njbox = int( ntj )
    x2    = ntj - njbox
    ncbox = int( nconji )
    x3    = nconji - ncbox
    ndbox = int( nconjj )
    x4    = nconjj - ndbox

    ibox = 1+t%nx*(t%ny*(t%nz*ndbox+ncbox)+njbox)+nibox

    fcc    = 0.0
    dfccdi = 0.0
    dfccdj = 0.0
    dfccdc = 0.0
    dfccdd = 0.0
    do i = 4, 1, -1
       sfcc   = 0.0
       sfccdj = 0.0
       sfccdc = 0.0
       sfccdd = 0.0
       do j = 4, 1, -1
          tfcc   = 0.0
          tfccdc = 0.0
          tfccdd = 0.0
          do k = 4, 1, -1
             ufcc   = 0.0
             ufccdd = 0.0
             do l=4,1,-1
                            coefij = t%coeff(ibox,i,j,k,l)
                            ufcc   =   ufcc*x4+       coefij
                if (l > 1)  ufccdd = ufccdd*x4+ (l-1)*coefij
             enddo
                         tfcc   =   tfcc*x3+       ufcc
             if (k > 1)  tfccdc = tfccdc*x3+ (k-1)*ufcc
                         tfccdd = tfccdd*x3+       ufccdd
          enddo
                       sfcc   = sfcc   *x2+       tfcc
          if  (j > 1)  sfccdj = sfccdj *x2+ (j-1)*tfcc
                       sfccdc = sfccdc *x2+       tfccdc
                       sfccdd = sfccdd *x2+       tfccdd
       enddo
                   fcc    = fcc    *x1+       sfcc
       if (i > 1)  dfccdi = dfccdi *x1+ (i-1)*sfcc
                   dfccdj = dfccdj *x1+       sfccdj
                   dfccdc = dfccdc *x1+       sfccdc
                   dfccdd = dfccdd *x1+       sfccdd
    enddo

  endsubroutine table4d_f

endmodule table4d
