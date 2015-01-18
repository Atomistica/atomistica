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
! Thread safe (OpenMP) random number generator
!**********************************************************************
module rng
  use system_module
  use units_module

#ifdef _OPENMP
  use omp_lib
#endif

  private

  integer(4)               :: indexf, indexb, buffer(250)
  !$omp threadprivate(indexf, indexb, buffer)

  integer(4), external     :: r250_irand

!  integer(4)               :: MAX_RNG_VAL = (ishl(1, 8*4-1)-1)
  real(DP), parameter      :: MAX_RNG_VAL = 2147483647._DP

  logical                  :: rng_initialized  = .false.

  public :: rng_init, rng_int4, rng_uniform, rng_normal1, rng_initialized
  public :: gaucorr

contains

  !**********************************************************************
  ! Initializes the subroutines rng_uniform by 
  ! creating random seed according to time. Note: This subrou-
  ! tine is called ONLY ONCE per run.
  !**********************************************************************
  subroutine rng_init(sd)
    implicit none

    integer(4), intent(in), optional :: sd

    ! ---

    integer(4)  :: l_sd

    integer     :: now(8)

    ! ---

    if (.not. rng_initialized) then

       if (present(sd)) then
          l_sd  = sd
       else
          call date_and_time(values=now)
          l_sd  = now(7)+1
       endif

#ifdef _OPENMP
       !$omp parallel
       call r250_init(l_sd*(omp_get_thread_num()+1), indexf, indexb, buffer)
       !$omp end parallel
#else
       call r250_init(l_sd, indexf, indexb, buffer)
#endif

       rng_initialized  = .true.

    endif

  endsubroutine rng_init


  !**********************************************************************
  ! Generates a single integer random number
  !**********************************************************************
  function rng_int4() result(rnd)
    implicit none

    integer(4)           :: rnd

    ! ---

    rnd  = r250_irand(indexf, indexb, buffer)

  endfunction rng_int4


  !**********************************************************************
  ! Generates a single random number belonging to [a,b].
  ! On first call, the function calls InitRand to set
  ! the random seed
  !**********************************************************************
  function rng_uniform(a, b) result(rnd)
    implicit none

    real(DP), intent(in), optional :: a, b

    real(DP)                       :: rnd

    ! ---

    integer(4)  :: i

    real(DP)    :: lower, upper

    ! ---

    i = r250_irand(indexf, indexb, buffer)

    !scale number to be between [a,b]
    lower = 0.0_DP
    upper = 1.0_DP
    if (present(a)) lower = a
    if (present(b)) upper = b

    rnd = lower + real((upper-lower)*i, DP)/MAX_RNG_VAL

  endfunction rng_uniform


  !**********************************************************************
  ! Generates a single random number from a normal distribution with
  ! unit variance.
  !**********************************************************************
  function rng_normal1() result(rnd)
    implicit none

    real(DP)             :: rnd

    ! ---

    integer(4)  :: i1, i2
    real        :: x1, x2, w

    ! ---

    w = 1.1_DP
    do while (w >= 1.0_DP)      
       i1  = r250_irand(indexf, indexb, buffer)
       i2  = r250_irand(indexf, indexb, buffer)

       x1  = real(i1, DP)/MAX_RNG_VAL
       x2  = real(i2, DP)/MAX_RNG_VAL

       x1 = 2.0_DP * x1 - 1.0_DP
       x2 = 2.0_DP * x2 - 1.0_DP
       w  = x1*x1 + x2*x2
    enddo

    w    = sqrt( (-2.0_DP * log(w))/w )
    rnd  = x1*w

  endfunction rng_normal1


  !**********************************************************************
  ! Produces two Gaussian random variables that are correlated (cov)
  !**********************************************************************
  subroutine gaucorr(x1, x2, sig1, sig2, cov)
    implicit none

    real(DP), intent(out)  :: x1, x2
    real(DP), intent(in)   :: sig1, sig2, cov

    ! ---

    real(DP)  :: eta1, eta2, a, hilf, b

    ! ---

    if (sig1 < 1d-12 .and. sig2 < 1d-12) then

       x1 = 0.0_DP
       x2 = 0.0_DP

    else

       a = 2.*PI*rng_uniform()
       hilf = rng_uniform()

       !**beschneide die verteilung**
       if (hilf.lt.0.0111) hilf=0.0111
       b    = sqrt(-2.*log(hilf))
       eta1 = b*sin(a)
       eta2 = b*cos(a)
       hilf = sig2**2-cov**2/sig1**2
       x1   = sig1*eta1

       !
       ! Sometimes, hilf can be smaller than 0 (numerical problem)
       !

       if (hilf <= 0) then
          x2   = cov*eta1/sig1
       else
          x2   = cov*eta1/sig1+sqrt(hilf)*eta2
       endif

    endif

  endsubroutine gaucorr

endmodule rng
