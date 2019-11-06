!! ======================================================================
!! Atomistica - Interatomic potential library and molecular dynamics code
!! https://github.com/Atomistica/atomistica
!!
!! Copyright (2005-2020) Lars Pastewka <lars.pastewka@imtek.uni-freiburg.de>
!! and others. See the AUTHORS file in the top-level Atomistica directory.
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
!! Smooth Particle-Mesh-Ewald kernel
!!
!! Smooth Particle-Mesh-Ewald kernel
!!
!! Code adopted from ORAC under GPL license,
!! http://www.chim.unifi.it/orac/
!<

#include "macros.inc"

module pme_kernel
#ifdef _OPENMP
  use omp_lib
#endif

  use supplib

#ifdef HAVE_FFTW3
  use, intrinsic :: iso_c_binding
  include 'fftw3.f03'
#endif

  private

  public :: pme_grid_t
  type pme_grid_t

     !
     ! PME stuff
     !

     integer                   :: numatoms, order, nfft1, nfft2, nfft3
     integer                   :: nfftdim1, nfftdim2, nfftdim3, ntheta

     !
     ! PME permanent space
     !

     real(DP), allocatable     :: bsp_mod1(:), bsp_mod2(:), bsp_mod3(:)

     !
     ! PME scratch space
     !

     real(DP), allocatable     :: theta1(:), theta2(:), theta3(:)
     real(DP), allocatable     :: dtheta1(:), dtheta2(:), dtheta3(:)
     real(DP), allocatable     :: fr1(:), fr2(:), fr3(:)

#ifdef HAVE_FFTW3
     logical                   :: fftw_is_initialized = .false.
     type(C_PTR)               :: plan_forward, plan_backward
     type(C_PTR)               :: Q_ptr
     complex(C_DOUBLE_COMPLEX), pointer  :: Q(:, :, :)
#else
     integer                   :: nfftable, nffwork
     integer                   :: sizfftable, sizffwork
     real(DP), allocatable     :: fftable(:), ffwork(:)
     complex(DP), allocatable  :: Q(:, :, :)
#endif

  endtype pme_grid_t


  public :: init
  interface init
     module procedure pme_grid_init
  endinterface

  public :: del
  interface del
     module procedure pme_grid_del
  endinterface

  public :: potential_and_field
  interface potential_and_field
     module procedure pme_grid_potential_and_field
  endinterface

contains

  subroutine pme_grid_init(this, grid, numatoms, order, error)
    implicit none

    type(pme_grid_t),  intent(inout)  :: this
    integer,           intent(in)     :: grid(3)
    integer,           intent(in)     :: numatoms
    integer,           intent(in)     :: order
    integer, optional, intent(out)    :: error

    ! ---

#ifndef HAVE_FFTW3
    integer  :: sfft,sffw
    real(DP)  :: dummy(1)
#endif

    ! ---

    INIT_ERROR(error)

    call del(this)

! INPUT  
!      nfft1,nfft2,nfft3,numatoms,order
!      nfft1,nfft2,nfft3 are the dimensions of the charge grid array
!      numatoms is number of atoms
!      order is the order of B-spline interpolation

! OUTPUT
!      nfftable,nffwork,ntheta,nQ
!      nfftable is permanent 3d fft table storage
!      nffwork is temporary 3d fft work storage
!      ntheta is size of arrays theta1-3 dtheta1-3
!      nheap is total size of permanent storage
!      nstack is total size of temporary storage

! This routine computes the above output parameters needed for 
! heap or stack allocation.

    this%numatoms = numatoms
    this%order = order
    this%nfft1 = grid(1)
    this%nfft2 = grid(2)
    this%nfft3 = grid(3)

    this%ntheta = this%numatoms*this%order

#ifdef HAVE_FFTW3
#ifdef _OPENMP
    write (ilog, '(5X,A)')  "Using multi-threaded (OpenMP) FFTW3 with " // omp_get_max_threads() // " threads."
    if (fftw_init_threads() == 0) then
      RAISE_ERROR("Could not initialize multi-threaded FFTW3.", error)
    endif
    call fftw_plan_with_nthreads(omp_get_max_threads())
#endif

    this%nfftdim1 = this%nfft1
    this%nfftdim2 = this%nfft2
    this%nfftdim3 = this%nfft3
    this%Q_ptr = fftw_alloc_complex( &
         int(this%nfft1*this%nfft2*this%nfft3, C_SIZE_T))
    call c_f_pointer(this%Q_ptr, this%Q, [this%nfft1,this%nfft2,this%nfft3])
#else
    call get_fftdims(this%nfft1, this%nfft2, this%nfft3, &
         this%nfftdim1, this%nfftdim2, this%nfftdim3, &
         this%nfftable, this%nffwork, this%sizfftable, this%sizffwork)

    allocate(this%Q(this%nfftdim1, this%nfftdim2, this%nfftdim3))
    allocate(this%fftable(this%sizfftable), this%ffwork(this%sizffwork))
#endif

    allocate(this%bsp_mod1(this%nfft1), this%bsp_mod2(this%nfft2))
    allocate(this%bsp_mod3(this%nfft3))

    allocate(this%theta1(this%ntheta), this%theta2(this%ntheta))
    allocate(this%theta3(this%ntheta), this%dtheta1(this%ntheta))
    allocate(this%dtheta2(this%ntheta), this%dtheta3(this%ntheta))

    allocate(this%fr1(this%numatoms), this%fr2(this%numatoms))
    allocate(this%fr3(this%numatoms))

    call load_bsp_moduli(this%bsp_mod1,this%bsp_mod2,this%bsp_mod3, &
         this%nfft1, this%nfft2, this%nfft3, this%order)
#ifdef HAVE_FFTW3
    this%plan_forward = fftw_plan_dft_3d( &
         this%nfft3, this%nfft2, this%nfft1, &
         this%Q, this%Q, &
         FFTW_FORWARD, FFTW_ESTIMATE)
    this%plan_backward = fftw_plan_dft_3d( &
         this%nfft3, this%nfft2, this%nfft1, &
         this%Q, this%Q, &
         FFTW_BACKWARD, FFTW_ESTIMATE)
    this%fftw_is_initialized = .true.
#else
    call fft_setup(dummy, this%fftable, this%ffwork, &
         this%nfft1, this%nfft2, this%nfft3, &
         this%nfftdim1, this%nfftdim2, this%nfftdim3, &
         this%nfftable, this%nffwork)
#endif

  endsubroutine pme_grid_init



  subroutine pme_grid_del(this)
    implicit none

    type(pme_grid_t), intent(inout)  :: this
    
    ! ---

    if (allocated(this%bsp_mod1))  deallocate(this%bsp_mod1)
    if (allocated(this%bsp_mod2))  deallocate(this%bsp_mod2)
    if (allocated(this%bsp_mod3))  deallocate(this%bsp_mod3)

    if (allocated(this%theta1))  deallocate(this%theta1)
    if (allocated(this%theta2))  deallocate(this%theta2)
    if (allocated(this%theta3))  deallocate(this%theta3)
        
    if (allocated(this%dtheta1))  deallocate(this%dtheta1)
    if (allocated(this%dtheta2))  deallocate(this%dtheta2)
    if (allocated(this%dtheta3))  deallocate(this%dtheta3)

    if (allocated(this%fr1))  deallocate(this%fr1)
    if (allocated(this%fr2))  deallocate(this%fr2)
    if (allocated(this%fr3))  deallocate(this%fr3)

#ifdef HAVE_FFTW3
    if (this%fftw_is_initialized) then
       call fftw_destroy_plan(this%plan_forward)
       call fftw_destroy_plan(this%plan_backward)
       this%Q => NULL()
       call fftw_free(this%Q_ptr)
#ifdef _OPENMP
       call fftw_cleanup_threads()
#endif
       this%fftw_is_initialized = .false.
    endif
#else
    if (allocated(this%Q))  deallocate(this%Q)
    if (allocated(this%fftable)) deallocate(this%fftable)
    if (allocated(this%ffwork)) deallocate(this%ffwork)
#endif

  endsubroutine pme_grid_del



  subroutine pme_grid_potential_and_field(this, &
       x, y, z, charge, recip, volume, ewald_coeff, &
       eer, virial, phi, dx, dy, dz)
    implicit none

! INPUT 
!       numatoms:  number of atoms
!       x,y,z:   atomic coords
!       charge  atomic charges
!       recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
!       volume: the volume of the unit cell
!       ewald_coeff:   ewald convergence parameter
!       order: the order of Bspline interpolation. E.g. cubic is order 4
!          fifth degree is order 6 etc. The order must be an even number 
!          and at least 4.
!       nfft1,nfft2,nfft3: the dimensions of the charge grid array

    type(pme_grid_t),   intent(inout)  :: this
    real(DP),           intent(in)     :: x(this%numatoms)
    real(DP),           intent(in)     :: y(this%numatoms)
    real(DP),           intent(in)     :: z(this%numatoms)
    real(DP),           intent(in)     :: charge(this%numatoms)
    real(DP),           intent(in)     :: recip(3,3), volume, ewald_coeff

! OUTPUT
!       eer:  ewald reciprocal or k-space  energy
!       phi:  potential incremented by k-space sum
!       dx,dy,dz:  field incremented by k-space sum
!       virial:  virial due to k-space sum (valid for atomic scaling;
!                rigid molecule virial needs a correction term not
!                computed here
    real(DP),           intent(out)    :: eer
    real(DP),           intent(out)    :: virial(3,3)
    real(DP),           intent(inout)  :: phi(this%numatoms)
    real(DP), optional, intent(inout)  :: dx(this%numatoms)
    real(DP), optional, intent(inout)  :: dy(this%numatoms)
    real(DP), optional, intent(inout)  :: dz(this%numatoms)

    ! ---

    call get_scaled_fractionals( &
         this%numatoms,x,y,z,recip,this%nfft1,this%nfft2,this%nfft3, &
         this%fr1,this%fr2,this%fr3)
    call get_bspline_coeffs( &
         this%numatoms,this%fr1,this%fr2,this%fr3,this%order, &
         this%theta1,this%theta2,this%theta3, &
         this%dtheta1,this%dtheta2,this%dtheta3)
    call fill_charge_grid( &
         this%numatoms,charge,this%theta1,this%theta2,this%theta3, &
         this%fr1,this%fr2,this%fr3,this%order, &
         this%nfft1,this%nfft2,this%nfft3, &
         this%nfftdim1,this%nfftdim2,this%nfftdim3, &
         this%Q)
#ifdef HAVE_FFTW3
    call fftw_execute_dft(this%plan_backward, this%Q, this%Q)
#else
    call fft_back( &
         this%Q,this%fftable,this%ffwork,this%nfft1,this%nfft2,this%nfft3, &
         this%nfftdim1,this%nfftdim2,this%nfftdim3,this%nfftable,this%nffwork)
#endif
    call energy_and_virial_sum( &
         this%Q,ewald_coeff,volume,recip, &
         this%bsp_mod1,this%bsp_mod2,this%bsp_mod3, &
         this%nfft1,this%nfft2,this%nfft3, &
         this%nfftdim1,this%nfftdim2,this%nfftdim3, &
         eer,virial)
#ifdef HAVE_FFTW3
    call fftw_execute_dft(this%plan_forward, this%Q, this%Q)
#else
    call fft_forward( &
         this%Q,this%fftable,this%ffwork,this%nfft1,this%nfft2,this%nfft3, &
         this%nfftdim1,this%nfftdim2,this%nfftdim3,this%nfftable,this%nffwork)
#endif

    if (present(dx) .and. present(dy) .and. present(dz)) then

       call potential_and_field_sum( &
            this%numatoms,recip,volume, &
            this%theta1,this%theta2,this%theta3, &
            this%dtheta1,this%dtheta2,this%dtheta3, &
            phi,dx,dy,dz, &
            this%fr1,this%fr2,this%fr3, &
            this%order, &
            this%nfft1,this%nfft2,this%nfft3, &
            this%nfftdim1,this%nfftdim2,this%nfftdim3, &
            this%Q)

    else

       call potential_sum( &
            this%numatoms,volume, &
            this%theta1,this%theta2,this%theta3, &
            phi, &
            this%fr1,this%fr2,this%fr3, &
            this%order, &
            this%nfft1,this%nfft2,this%nfft3, &
            this%nfftdim1,this%nfftdim2,this%nfftdim3, &
            this%Q)

    endif

  endsubroutine pme_grid_potential_and_field



  subroutine get_scaled_fractionals(numatoms,x,y,z,recip,nfft1,nfft2,nfft3, &
       fr1,fr2,fr3)
    implicit none

! INPUT:
!      numatoms: number of atoms
!      x,y,z: arrays of cartesian coords
!      recip: the 3x3 array of reciprocal vectors stored as columns
! OUTPUT:
!     fr1,fr2,fr3 the scaled and shifted fractional coords

    integer   :: numatoms,nfft1,nfft2,nfft3
    real(DP)  :: x(numatoms),y(numatoms),z(numatoms),recip(3,3)
    real(DP)  :: fr1(numatoms),fr2(numatoms),fr3(numatoms)

    ! ---

    integer   :: n
    real(DP)  :: w

    ! ---

    !$omp  parallel do default(none) &
    !$omp& shared(x, y, z, fr1, fr2, fr3) &
    !$omp& firstprivate(recip, nfft1, nfft2, nfft3, numatoms) &
    !$omp& private(w)
    do n = 1,numatoms
       w = x(n)*recip(1,1)+y(n)*recip(2,1)+z(n)*recip(3,1)
       fr1(n) = nfft1*(w - anint(w) + 0.5d0)
       w = x(n)*recip(1,2)+y(n)*recip(2,2)+z(n)*recip(3,2)
       fr2(n) = nfft2*(w - anint(w) + 0.5d0)
       w = x(n)*recip(1,3)+y(n)*recip(2,3)+z(n)*recip(3,3)
       fr3(n) = nfft3*(w - anint(w) + 0.5d0)
    enddo
    
  endsubroutine get_scaled_fractionals



  subroutine load_bsp_moduli(bsp_mod1,bsp_mod2,bsp_mod3, &
       nfft1,nfft2,nfft3,order)
    implicit none

    integer   :: nfft1,nfft2,nfft3,order
    real(DP)  :: bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)

    ! ---

    integer, parameter  :: MAXORDER  = 25
    integer, parameter  :: MAXNFFT   = 1000

    real(DP)  :: array(MAXORDER),darray(MAXORDER),w
    real(DP)  :: bsp_arr(MAXNFFT)
    integer   :: i,maxn

! this routine loads the moduli of the inverse DFT of the B splines
! bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
! Order is the order of the B spline approx.

    if ( order .gt. MAXORDER )then
       write(6,*)'order too large! check on MAXORDER'
       stop
    endif
    maxn = max(nfft1,nfft2,nfft3)
    if ( maxn .gt. MAXNFFT )then 
       write(6,*)'nfft1-3 too large! check on MAXNFFT'
       stop
    endif
    w = 0.d0
    call fill_bspline(w,order,array,darray)
    do i = 1,maxn
       bsp_arr(i) = 0.d0
    enddo
    do i = 2,order+1
       bsp_arr(i) = array(i-1)
    enddo
    call dftmod(bsp_mod1,bsp_arr,nfft1)
    call dftmod(bsp_mod2,bsp_arr,nfft2)
    call dftmod(bsp_mod3,bsp_arr,nfft3)

  endsubroutine load_bsp_moduli


    
  subroutine dftmod(bsp_mod,bsp_arr,nfft)
    implicit none

    integer   :: nfft
    real(DP)  :: bsp_mod(nfft),bsp_arr(nfft)

! Computes the modulus of the discrete fourier transform of bsp_arr,
!  storing it into bsp_mod

    integer   :: j,k
    real(DP)  :: sum1,sum2,twopi,arg,tiny

      ! ---

!      twopi = 2.d0*3.14159265358979323846
    twopi = 2.0_DP * PI
    tiny = 1.d-7
    do k = 1,nfft
       sum1 = 0.d0
       sum2 = 0.d0
       do j = 1,nfft
          arg = twopi*(k-1)*(j-1)/nfft
          sum1 = sum1 + bsp_arr(j)*dcos(arg)
          sum2 = sum2 + bsp_arr(j)*dsin(arg)
       enddo
       bsp_mod(k) = sum1**2 + sum2**2
    enddo
    do k = 1,nfft
       if ( bsp_mod(k) .lt. tiny ) &
            bsp_mod(k) = 0.5d0*(bsp_mod(k-1) + bsp_mod(k+1))
    enddo

  endsubroutine dftmod



  subroutine fill_charge_grid( &
       numatoms,charge,theta1,theta2,theta3,fr1,fr2,fr3, &
       order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)

!---------------------------------------------------------------------
! INPUT:
!      numatoms:  number of atoms
!      charge: the array of atomic charges
!      theta1,theta2,theta3: the spline coeff arrays
!      fr1,fr2,fr3 the scaled and shifted fractional coords
!      nfft1,nfft2,nfft3: the charge grid dimensions
!      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
!      order: the order of spline interpolation
! OUTPUT:
!      Q the charge grid
!---------------------------------------------------------------------

    implicit none
    
    integer      :: numatoms,order,nfft1,nfft2,nfft3
    integer      :: nfftdim1,nfftdim2,nfftdim3
    real(DP)     :: fr1(numatoms),fr2(numatoms),fr3(numatoms)
    real(DP)     :: theta1(order,numatoms),theta2(order,numatoms)
    real(DP)     :: theta3(order,numatoms),charge(numatoms)
    complex(DP)  :: Q(nfftdim1,nfftdim2,nfftdim3)

    ! ---

    integer   :: n,ith1,ith2,ith3,i0,j0,k0,i,j,k
    real(DP)  :: prod

    ! ---

    Q = 0.0_DP
      
    do n = 1,numatoms
       k0 = int(fr3(n)) - order
       do ith3 = 1,order
          k0 = k0 + 1
          k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
          j0 = int(fr2(n)) - order
          do ith2 = 1,order
             j0 = j0 + 1
             j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
             prod = theta2(ith2,n)*theta3(ith3,n)*charge(n)
             i0 = int(fr1(n)) - order
             do ith1 = 1,order
                i0 = i0 + 1
                i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
                Q(i,j,k) = Q(i,j,k) + theta1(ith1,n) * prod
             enddo
          enddo
       enddo
    enddo
    
  endsubroutine fill_charge_grid



  subroutine get_bspline_coeffs( &
       numatoms,fr1,fr2,fr3,order, &
       theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
!---------------------------------------------------------------------
! INPUT:
!      numatoms: number of atoms
!      fr1,fr2,fr3 the scaled and shifted fractional coords
!      order: the order of spline interpolation
! OUTPUT
!      theta1,theta2,theta3: the spline coeff arrays
!      dtheta1,dtheta2,dtheta3: the 1st deriv of spline coeff arrays
!---------------------------------------------------------------------
    implicit none

    integer   :: numatoms,order
    real(DP)  :: fr1(numatoms),fr2(numatoms),fr3(numatoms)
    real(DP)  :: theta1(order,numatoms),theta2(order,numatoms)
    real(DP)  :: theta3(order,numatoms),dtheta1(order,numatoms)
    real(DP)  :: dtheta2(order,numatoms),dtheta3(order,numatoms)

    ! ---

    real(DP)  :: w
    integer   :: n

    !$omp  parallel do default(none) &
    !$omp& shared(fr1, fr2, fr3) &
    !$omp& shared(theta1, dtheta1, theta2, dtheta2, theta3, dtheta3) &
    !$omp& firstprivate(order, numatoms) &
    !$omp& private(w)
    do n = 1,numatoms
       w = fr1(n)-int(fr1(n))
       call fill_bspline(w,order,theta1(1,n),dtheta1(1,n))
       w = fr2(n)-int(fr2(n))
       call fill_bspline(w,order,theta2(1,n),dtheta2(1,n))
       w = fr3(n)-int(fr3(n))
       call fill_bspline(w,order,theta3(1,n),dtheta3(1,n))
    enddo

  endsubroutine get_bspline_coeffs



  subroutine fill_bspline(w,order,array,darray)
! use standard B-spline recursions: see doc file
    implicit none

    integer   :: order
    real(DP)  :: w,array(order),darray(order)
    
    ! ---

    integer  :: k

    ! ---

! do linear case
    call bsp_init(array,w,order)
! compute standard b-spline recursion
    do k = 3,order-1
       call bsp_one_pass(array,w,k)
    enddo
! perform standard b-spline differentiation
    call bsp_diff(array,darray,order)
! one more recursion
    call bsp_one_pass(array,w,order)

  endsubroutine fill_bspline



  subroutine bsp_init(c,x,order)
    implicit none

    integer   :: order
    real(DP)  :: c(order),x
    
    ! ---

    c(order) = 0.d0
    c(2) = x
    c(1) = 1.d0 - x

  endsubroutine bsp_init



  subroutine bsp_one_pass(c,x,k)
    implicit none

    real(DP)  :: c(*),x
    integer   :: k

    ! ---

    real(DP)  :: div
    integer   :: j

    ! ---

    div = 1.d0 / (k-1)
    c(k) = div*x*c(k-1)
    do j = 1,k-2
       c(k-j) = div*((x+j)*c(k-j-1) + (k-j-x)*c(k-j))
    enddo
    c(1) = div*(1-x)*c(1)
   
  endsubroutine bsp_one_pass



  subroutine bsp_diff(c,d,order)
    implicit none

    real(DP)  :: c(*),d(*)
    integer   :: order

    ! ---

    integer  :: j

    ! ---

    d(1) = -c(1)
    do j = 2,order
       d(j) = c(j-1) - c(j)
    enddo

  endsubroutine bsp_diff



  subroutine energy_and_virial_sum( &
       Q,ewaldcof,volume,recip,bsp_mod1,bsp_mod2,bsp_mod3, &
       nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,eer,vir)
    implicit none

    integer      :: nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3
    complex(DP)  :: Q(nfftdim1,nfftdim2,nfftdim3)
    real(DP)     :: bsp_mod1(nfft1),bsp_mod2(nfft2),bsp_mod3(nfft3)
    real(DP)     :: ewaldcof,volume
    real(DP)     :: eer,vir(3,3)
    real(DP)     :: recip(3,3)

    ! ---

    real(DP)  :: fac,denom,eterm,vterm,energy
    integer   :: k1,k2,k3,m1,m2,m3,nff,ind,jnd,indtop
    integer   :: nf1,nf2,nf3
    real(DP)  :: mhat1,mhat2,mhat3,msq,struc2

    ! ---

    indtop = nfft1*nfft2*nfft3
    fac = pi**2/ewaldcof**2
    nff = nfft1*nfft2
    nf1 = nfft1/2
    if ( 2*nf1 .lt. nfft1 )nf1 = nf1+1
    nf2 = nfft2/2
    if ( 2*nf2 .lt. nfft2 )nf2 = nf2+1
    nf3 = nfft3/2
    if ( 2*nf3 .lt. nfft3 )nf3 = nf3+1
    energy = 0.d0
    DO k1 = 1,3
       DO k2 = 1,3
          vir(k1,k2) = 0.0D0
       END DO
    END DO

    do ind = 1,indtop-1

! get k1,k2,k3 from the relationship
!           ind = (k1-1) + (k2-1)*nfft1 + (k3-1)*nfft2*nfft1

       k3 = ind/nff + 1
       jnd = ind - (k3-1)*nff
       k2 = jnd/nfft1 + 1
       k1 = jnd - (k2-1)*nfft1 +1
       m1 = k1 - 1
       if ( k1 .gt. nf1 )m1 = k1 - 1 - nfft1
       m2 = k2 - 1
       if ( k2 .gt. nf2 )m2 = k2 - 1 - nfft2
       m3 = k3 - 1
       if ( k3 .gt. nf3 )m3 = k3 - 1 - nfft3
       mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
       mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
       mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
       msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
       denom = bsp_mod1(k1)*bsp_mod2(k2)*bsp_mod3(k3)*msq
       eterm = dexp(-fac*msq)/denom
       vterm = 2.d0*(fac*msq + 1.d0)/msq
       struc2 = Q(k1,k2,k3)*conjg(Q(k1,k2,k3))
       energy = energy + eterm * struc2
       vir(1,1) = vir(1,1) + eterm * struc2 * (vterm*mhat1*mhat1 - 1.d0)
       vir(1,2) = vir(1,2) + eterm * struc2 * (vterm*mhat1*mhat2)
       vir(1,3) = vir(1,3) + eterm * struc2 * (vterm*mhat1*mhat3)
       vir(2,2) = vir(2,2) + eterm * struc2 * (vterm*mhat2*mhat2 - 1.d0)
       vir(2,3) = vir(2,3) + eterm * struc2 * (vterm*mhat2*mhat3)
       vir(3,3) = vir(3,3) + eterm * struc2 * (vterm*mhat3*mhat3 - 1.d0)
       Q(k1,k2,k3) = eterm * Q(k1,k2,k3)

    enddo

    fac = 1.0_DP/(2*pi*volume)
    eer = fac*energy
    vir(2,1)=vir(1,2)
    vir(3,1)=vir(1,3)
    vir(3,2)=vir(2,3)
   
    DO k1 = 1,3
       DO k2 = 1,3
          vir(k1,k2) = fac*vir(k1,k2)
       END DO
    END DO

  endsubroutine energy_and_virial_sum



  subroutine potential_and_field_sum( &
       numatoms,recip,volume,theta1,theta2,theta3, &
       dtheta1,dtheta2,dtheta3,phi,Ex,Ey,Ez,fr1,fr2,fr3, &
       order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
    implicit none

    integer      :: numatoms,order,nfft1,nfft2,nfft3
    integer      :: nfftdim1,nfftdim2,nfftdim3
    real(DP)     :: recip(3,3), volume
    real(DP)     :: fr1(numatoms),fr2(numatoms),fr3(numatoms)
    real(DP)     :: phi(numatoms),Ex(numatoms),Ey(numatoms),Ez(numatoms)
    real(DP)     :: theta1(order,numatoms),theta2(order,numatoms)
    real(DP)     :: theta3(order,numatoms)
    real(DP)     :: dtheta1(order,numatoms),dtheta2(order,numatoms)
    real(DP)     :: dtheta3(order,numatoms)
    complex(DP)  :: Q(nfftdim1,nfftdim2,nfftdim3)

    ! ---

    integer   :: n,ith1,ith2,ith3,i0,j0,k0,i,j,k
    real(DP)  :: f0,f1,f2,f3,term,fac
!      real(DP)  :: pi

!$DOACROSS LOCAL(f1,f2,f3,k0,k,j0,j,i0,i,term,n,ith1,ith2,ith3),
!$&  SHARE(numatoms,fr1,fr2,fr3,charge,Q,Ex,Ey,Ez,recip,order,
!$&   nfft1,nfft2,nfft3,theta1,theta2,theta3,dtheta1,dtheta2,dtheta3)
!      pi = 3.14159265358979323846d0

    fac = 1.0_DP/(pi*volume)
    !$omp  parallel do default(none) &
    !$omp& shared(fr1, fr2, fr3, phi, Ex, Ey, Ez, Q) &
    !$omp& shared(theta1, dtheta1, theta2, dtheta2, theta3, dtheta3) &
    !$omp& firstprivate(fac, nfft1, nfft2, nfft3, order, recip, numatoms) &
    !$omp& private(n, ith1, ith2, ith3, i0, j0, k0, i, j, k) &
    !$omp& private(f0, f1, f2, f3, term)
    do n = 1,numatoms
       f0 = 0.d0
       f1 = 0.d0
       f2 = 0.d0
       f3 = 0.d0
       k0 = int(fr3(n)) - order
       do ith3 = 1,order
          k0 = k0 + 1
          k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
          j0 = int(fr2(n)) - order
          do ith2 = 1,order
             j0 = j0 + 1
             j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
             i0 = int(fr1(n)) - order
             do ith1 = 1,order
                i0 = i0 + 1
                i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
! --- pas
                  !           term = charge(n)*Q(1,i,j,k)
                term = real(Q(i,j,k), DP)
! force is negative of grad
                f0 = f0 + term * theta1(ith1,n) * &
                     theta2(ith2,n) * theta3(ith3,n)
                f1 = f1 - nfft1 * term * dtheta1(ith1,n) * &
                     theta2(ith2,n) * theta3(ith3,n)
                f2 = f2 - nfft2 * term * theta1(ith1,n) * &
                     dtheta2(ith2,n) * theta3(ith3,n)
                f3 = f3 - nfft3 * term * theta1(ith1,n) * &
                     theta2(ith2,n) * dtheta3(ith3,n)
             enddo
          enddo
       enddo
       phi(n) = phi(n) + fac*f0
       Ex(n) = Ex(n) + fac*(recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3)
       Ey(n) = Ey(n) + fac*(recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3)
       Ez(n) = Ez(n) + fac*(recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3)
    enddo

  endsubroutine potential_and_field_sum



  subroutine potential_sum( &
       numatoms,volume,theta1,theta2,theta3, &
       phi,fr1,fr2,fr3, &
       order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,Q)
    implicit none

    integer      :: numatoms,order,nfft1,nfft2,nfft3
    integer      :: nfftdim1,nfftdim2,nfftdim3
    real(DP)     :: volume
    real(DP)     :: fr1(numatoms),fr2(numatoms),fr3(numatoms)
    real(DP)     :: phi(numatoms)
    real(DP)     :: theta1(order,numatoms),theta2(order,numatoms)
    real(DP)     :: theta3(order,numatoms)
    complex(DP)  :: Q(nfftdim1,nfftdim2,nfftdim3)

    ! ---

    integer   :: n,ith1,ith2,ith3,i0,j0,k0,i,j,k
    real(DP)  :: f0,term,fac

    ! ---

    fac = 1.0_DP/(pi*volume)
    !$omp  parallel do default(none) &
    !$omp& shared(fr1, fr2, fr3, phi, Q) &
    !$omp& shared(theta1, theta2, theta3) &
    !$omp& firstprivate(fac, nfft1, nfft2, nfft3, order, numatoms) &
    !$omp& private(n, ith1, ith2, ith3, i0, j0, k0, i, j, k) &
    !$omp& private(f0, term)
    do n = 1,numatoms
       f0 = 0.d0
       k0 = int(fr3(n)) - order
       do ith3 = 1,order
          k0 = k0 + 1
          k = k0 + 1 + (nfft3 - isign(nfft3,k0))/2
          j0 = int(fr2(n)) - order
          do ith2 = 1,order
             j0 = j0 + 1
             j = j0 + 1 + (nfft2 - isign(nfft2,j0))/2
             i0 = int(fr1(n)) - order
             do ith1 = 1,order
                i0 = i0 + 1
                i = i0 + 1 + (nfft1 - isign(nfft1,i0))/2
                term = real(Q(i,j,k), DP)
                f0 = f0 + term * theta1(ith1,n) * &
                     theta2(ith2,n) * theta3(ith3,n)
             enddo
          enddo
       enddo
       phi(n) = phi(n) + fac*f0
    enddo

  endsubroutine potential_sum

endmodule pme_kernel
