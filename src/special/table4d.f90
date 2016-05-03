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
!! 4D cubic spline interpolation
!<

#include "macros.inc"

module table4d
  use libAtoms_module

  use logging, only: ilog

  private

  public :: table4d_t
  type table4d_t

     integer               :: nx = 1
     integer               :: ny = 1
     integer               :: nz = 1
     integer               :: nt = 1

     integer               :: nboxs

     real(DP), allocatable :: coeff(:, :, :, :, :)

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

  public :: eval
  interface eval
     module procedure table4d_eval
  endinterface

  public :: table4d_prlog

contains

  !>
  !! generates the coefficients for bicubic interpolation of fch(ni,nj)
  !! copyright: Keith Beardmore 30/11/93.
  !!            Lars Pastewka 05/07
  !<
  subroutine table4d_init(t, nx, ny, nz, nt, values, dvdx, dvdy, dvdz, dvdt, &
                          error)
    implicit none

    type(table4d_t), intent(inout)    :: t
    integer, intent(in)               :: nx
    integer, intent(in)               :: ny
    integer, intent(in)               :: nz
    integer, intent(in)               :: nt
    real(DP), intent(in)              :: values(0:, 0:, 0:, 0:)
    real(DP), optional, intent(in)    :: dvdx(0:nx, 0:ny, 0:nz, 0:nt)
    real(DP), optional, intent(in)    :: dvdy(0:nx, 0:ny, 0:nz, 0:nt)
    real(DP), optional, intent(in)    :: dvdz(0:nx, 0:ny, 0:nz, 0:nt)
    real(DP), optional, intent(in)    :: dvdt(0:nx, 0:ny, 0:nz, 0:nt)
    integer, intent(inout), optional  :: error

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

    ! Bounds checking

    if (lbound(values, 1) /= 0 .or. ubound(values, 1) /= nx) then
       RAISE_ERROR("First index of *values* must run from 0 to " // nx // ", but does run from " // lbound(values, 1) // " to " // ubound(values, 1) // ".", error)
    endif
    if (lbound(values, 2) /= 0 .or. ubound(values, 2) /= ny) then
       RAISE_ERROR("Second index of *values* must run from 0 to " // ny // ", but does run from " // lbound(values, 2) // " to " // ubound(values, 2) // ".", error)
    endif
    if (lbound(values, 3) /= 0 .or. ubound(values, 3) /= nz) then
       RAISE_ERROR("Third index of *values* must run from 0 to " // nz // ", but does run from " // lbound(values, 3) // " to " // ubound(values, 3) // ".", error)
    endif
    if (lbound(values, 4) /= 0 .or. ubound(values, 4) /= nt) then
       RAISE_ERROR("Fourth index of *values* must run from 0 to " // nt // ", but does run from " // lbound(values, 4) // " to " // ubound(values, 4) // ".", error)
    endif

    if (present(dvdx)) then
       if (lbound(dvdx, 1) /= 0 .or. ubound(dvdx, 1) /= nx) then
          RAISE_ERROR("First index of *dvdx* must run from 0 to " // nx // ", but does run from " // lbound(dvdx, 1) // " to " // ubound(dvdx, 1) // ".", error)
       endif
       if (lbound(dvdx, 2) /= 0 .or. ubound(dvdx, 2) /= ny) then
          RAISE_ERROR("Second index of *dvdx* must run from 0 to " // ny // ", but does run from " // lbound(dvdx, 2) // " to " // ubound(dvdx, 2) // ".", error)
       endif
       if (lbound(dvdx, 3) /= 0 .or. ubound(dvdx, 3) /= nz) then
          RAISE_ERROR("Third index of *dvdx* must run from 0 to " // nz // ", but does run from " // lbound(dvdx, 3) // " to " // ubound(dvdx, 3) // ".", error)
       endif
       if (lbound(dvdx, 4) /= 0 .or. ubound(dvdx, 4) /= nt) then
          RAISE_ERROR("Third index of *dvdx* must run from 0 to " // nt // ", but does run from " // lbound(dvdx, 3) // " to " // ubound(dvdx, 3) // ".", error)
       endif
    endif

    if (present(dvdy)) then
       if (lbound(dvdy, 1) /= 0 .or. ubound(dvdy, 1) /= nx) then
          RAISE_ERROR("First index of *dvdy* must run from 0 to " // nx // ", but does run from " // lbound(dvdy, 1) // " to " // ubound(dvdy, 1) // ".", error)
       endif
       if (lbound(dvdy, 2) /= 0 .or. ubound(dvdy, 2) /= ny) then
          RAISE_ERROR("Second index of *dvdy* must run from 0 to " // ny // ", but does run from " // lbound(dvdy, 2) // " to " // ubound(dvdy, 2) // ".", error)
       endif
       if (lbound(dvdy, 3) /= 0 .or. ubound(dvdy, 3) /= nz) then
          RAISE_ERROR("Third index of *dvdy* must run from 0 to " // nz // ", but does run from " // lbound(dvdy, 3) // " to " // ubound(dvdy, 3) // ".", error)
       endif
       if (lbound(dvdy, 4) /= 0 .or. ubound(dvdy, 4) /= nt) then
          RAISE_ERROR("Third index of *dvdy* must run from 0 to " // nt // ", but does run from " // lbound(dvdy, 3) // " to " // ubound(dvdy, 3) // ".", error)
       endif
    endif

    if (present(dvdz)) then
       if (lbound(dvdz, 1) /= 0 .or. ubound(dvdz, 1) /= nx) then
          RAISE_ERROR("First index of *dvdz* must run from 0 to " // nx // ", but does run from " // lbound(dvdz, 1) // " to " // ubound(dvdz, 1) // ".", error)
       endif
       if (lbound(dvdz, 2) /= 0 .or. ubound(dvdz, 2) /= ny) then
          RAISE_ERROR("Second index of *dvdz* must run from 0 to " // ny // ", but does run from " // lbound(dvdz, 2) // " to " // ubound(dvdz, 2) // ".", error)
       endif
       if (lbound(dvdz, 3) /= 0 .or. ubound(dvdz, 3) /= nz) then
          RAISE_ERROR("Third index of *dvdz* must run from 0 to " // nz // ", but does run from " // lbound(dvdz, 3) // " to " // ubound(dvdz, 3) // ".", error)
       endif
       if (lbound(dvdz, 4) /= 0 .or. ubound(dvdz, 4) /= nt) then
          RAISE_ERROR("Third index of *dvdz* must run from 0 to " // nt // ", but does run from " // lbound(dvdz, 3) // " to " // ubound(dvdz, 3) // ".", error)
       endif
    endif

    if (present(dvdt)) then
       if (lbound(dvdt, 1) /= 0 .or. ubound(dvdt, 1) /= nx) then
          RAISE_ERROR("First index of *dvdt* must run from 0 to " // nx // ", but does run from " // lbound(dvdt, 1) // " to " // ubound(dvdt, 1) // ".", error)
       endif
       if (lbound(dvdt, 2) /= 0 .or. ubound(dvdt, 2) /= ny) then
          RAISE_ERROR("Second index of *dvdt* must run from 0 to " // ny // ", but does run from " // lbound(dvdt, 2) // " to " // ubound(dvdt, 2) // ".", error)
       endif
       if (lbound(dvdt, 3) /= 0 .or. ubound(dvdt, 3) /= nz) then
          RAISE_ERROR("Third index of *dvdt* must run from 0 to " // nz // ", but does run from " // lbound(dvdt, 3) // " to " // ubound(dvdt, 3) // ".", error)
       endif
       if (lbound(dvdt, 4) /= 0 .or. ubound(dvdt, 4) /= nt) then
          RAISE_ERROR("Third index of *dvdt* must run from 0 to " // nt // ", but does run from " // lbound(dvdt, 3) // " to " // ubound(dvdt, 3) // ".", error)
       endif
    endif

    ! ---

    t%nx     = nx
    t%ny     = ny
    t%nz     = nz
    t%nt     = nt
    t%nboxs  = nx*ny*nz*nt

    if (allocated(t%coeff))  deallocate(t%coeff)
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

    B = 0.0_DP
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
                   if (present(dvdx)) then
                      B(irow+ ncorn  ,icol) = dvdx(nx1, nx2, nx3, nx4)
                   endif
                   if (present(dvdy)) then
                      B(irow+ 2*ncorn,icol) = dvdy(nx1, nx2, nx3, nx4)
                   endif
                   if (present(dvdz)) then
                      B(irow+ 3*ncorn,icol) = dvdz(nx1, nx2, nx3, nx4)
                   endif
                   if (present(dvdt)) then
                      B(irow+ 4*ncorn,icol) = dvdt(nx1, nx2, nx3, nx4)
                   endif
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
       RAISE_ERROR("dgesv failed.", error)
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

    if (allocated(t%coeff)) then
       deallocate(t%coeff)
    endif

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
  subroutine table4d_eval(t, nti, ntj, nconji, nconjj, fcc, dfccdi, dfccdj, &
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
    if (nibox < 0)      nibox = 0
    if (nibox >= t%nx)  nibox = t%nx-1
    njbox = int( ntj )
    if (njbox < 0)      njbox = 0
    if (njbox >= t%ny)  njbox = t%ny-1
    ncbox = int( nconji )
    if (ncbox < 0)      ncbox = 0
    if (ncbox >= t%nz)  ncbox = t%nz-1
    ndbox = int( nconjj )
    if (ndbox < 0)      ndbox = 0
    if (ndbox >= t%nt)  ndbox = t%nt-1

    ibox = 1+t%nx*(t%ny*(t%nz*ndbox+ncbox)+njbox)+nibox
    x1   = nti - nibox
    x2   = ntj - njbox
    x3   = nconji - ncbox
    x4   = nconjj - ndbox

    fcc    = 0.0_DP
    dfccdi = 0.0_DP
    dfccdj = 0.0_DP
    dfccdc = 0.0_DP
    dfccdd = 0.0_DP
    do i = 4, 1, -1
       sfcc   = 0.0_DP
       sfccdj = 0.0_DP
       sfccdc = 0.0_DP
       sfccdd = 0.0_DP
       do j = 4, 1, -1
          tfcc   = 0.0_DP
          tfccdc = 0.0_DP
          tfccdd = 0.0_DP
          do k = 4, 1, -1
             ufcc   = 0.0_DP
             ufccdd = 0.0_DP
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
    
  endsubroutine table4d_eval


  !>
  !! Print to screen
  !!
  !! Print to screen
  !<
  subroutine table4d_print(this, indent)
    implicit none

    type(table4d_t),   intent(in) :: this
    integer, optional, intent(in) :: indent

    ! ---

    call table4d_print_un(7, this)

  endsubroutine table4d_print


  !>
  !! Print to log file
  !!
  !! Print to log file
  !<
  subroutine table4d_prlog(this, indent)
    implicit none

    type(table4d_t),   intent(in) :: this
    integer, optional, intent(in) :: indent

    ! ---

    call table4d_print_un(ilog, this, indent)

  endsubroutine table4d_prlog


  !>
  !! Print to unit
  !!
  !! Print to unit
  !<
  subroutine table4d_print_un(un, this, indent)
    implicit none

    integer,           intent(in) :: un
    type(table4d_t),   intent(in) :: this
    integer, optional, intent(in) :: indent
    
    ! ---

    integer          :: i, j, k, l, m
    real(DP)         :: row(0:this%nx), val, dummy1, dummy2, dummy3, dummy4
    character(1000)  :: fmt, fmthdr, fmtstart

    ! ---

    if (present(indent)) then
       fmthdr = "(" // (indent) // "X,A15,I10," // this%nx // "I20)"
    else
       fmthdr = "(A15,I10," // this%nx // "I20)"
    endif

    if (present(indent)) then
       fmtstart = "(" // indent // "X,I3,' -'"
    else
       fmtstart = "(4I,1X"
    endif

    do l = 0, this%nt
       do k = 0, this%nz
          write (un, fmthdr)  "[:,:,"//k//","//l//"]", &
               (/ ( i, i=0, this%nx ) /)
          do j = 0, this%ny
             fmt = fmtstart
             m = 0
             do i = 0, this%nx
                call eval(this, i*1.0_DP, j*1.0_DP, k*1.0_DP, l*1.0_DP, &
                     val, dummy1, dummy2, dummy3, dummy4)
                if (abs(val) > 1e-12) then
                   fmt = trim(fmt) // ",ES20.10"
                   row(m) = val
                   m = m+1
                else
                   fmt = trim(fmt) // ",'       ----------   '"
                endif
             enddo
             fmt = trim(fmt) // ")"

             write (un, fmt)  j, row(0:m-1)
          enddo

          write (un, *)
       enddo

       write (un, *)
    enddo

  endsubroutine table4d_print_un

endmodule table4d
