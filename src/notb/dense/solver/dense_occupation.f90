!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
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

#include "macros.inc"
#include "filter.inc"

module dense_occupation
  use supplib

  use particles

  use dense_hamiltonian_type

  implicit none

  private

  public :: dense_solver_lapack_t

  type dense_solver_lapack_t

     !
     ! Parameters
     !

     logical(BOOL)  :: enabled      = .false.

     integer   :: solver_type  = 1        ! type of solver to use...
                                          ! 0 = LAPACK standard
                                          ! 1 = LAPACK divide-and-conquer
                                          ! 2 = LAPACK expert

     real(DP)  :: Tele         = 0.01_DP  ! Electronic temperature

     integer   :: norb         = -1       ! Number of orbitals
     integer   :: n_bands      = 10       ! Number of eigenvalues to solve for

     !
     ! The tight-binding object
     !

     type(dense_hamiltonian_t), pointer  :: tb  => NULL()

     !
     ! Solver stuff (eigenvectors and eigenvalues)
     !

     real(DP), allocatable  :: evals(:, :)                  ! eigenvalues
     WF_T(DP), allocatable  :: evecs(:, :, :)               ! eigenvectors

     real(DP), allocatable  :: f(:, :)                      ! occupation, i.e., the Fermi function

     !
     ! Work buffers
     !

     WF_T(DP), allocatable :: S2(:, :)

#ifdef COMPLEX_WF
     WF_T(DP), allocatable :: work(:)
     real(DP), allocatable :: rwork(:)
     integer,  allocatable :: iwork(:)
#else
     WF_T(DP), allocatable :: work(:)
     integer,  allocatable :: iwork(:)
#endif

  endtype dense_solver_lapack_t

  save
  WF_T(DP), allocatable  :: tr_evecs(:, :)

!  public :: construct_density_and_energy_matrix
  public :: tr_evecs, occupy

contains

  !**********************************************************************
  ! Fermi-Dirac distribution
  !**********************************************************************
  function fermi_dirac(e, mu, T)
    implicit none

    real(DP), intent(in)  :: e   ! energy
    real(DP), intent(in)  :: mu  ! chemical potential
    real(DP), intent(in)  :: T   ! temperature

    real(DP)              :: fermi_dirac

    ! ---

    real(DP)  :: x

    ! ---

    if (T < 1e-6) then

       if (e < mu) then
          fermi_dirac = 1.0
       else
          fermi_dirac = 0.0
       endif

    else

       x = (e-mu)/T

       if (x > 100.0) then
          fermi_dirac = 0.0
       else
          fermi_dirac = 1./(1. + exp(x))
       endif

    endif

  endfunction fermi_dirac


  !**********************************************************************
  ! Calculate the occupation numbers for a given number of occupied
  ! orbitals
  !**********************************************************************
  subroutine occupy(s, tb, noc, error)
    implicit none

    type(dense_solver_lapack_t), intent(inout)  :: s
    type(dense_hamiltonian_t),   intent(inout)  :: tb
    real(DP),                    intent(in)     :: noc
    integer,           optional, intent(out)    :: error

    ! ---

    real(DP)  :: mu
    integer   :: i, k

#define F s%f

    ! ---

    call timer_start('occu')

    INIT_ERROR(error)

    !
    ! Find the chemical potential
    !

    mu  = SolveMu(s, tb, s%Tele, 2*noc, error)
    PASS_ERROR(error)

    s%f = 0

    ! Fixme!!! This is O(N^3)!

    !
    ! Construct the occupation
    !

    do k = 1, tb%nk
       !$omp  parallel do default(none) &
       !$omp& shared(k, mu, s, tb)
       do i = 1, tb%norb
          F(i, k) = 2 * fermi_dirac(s%evals(i, k), mu, s%Tele)
       enddo
       !$omp end parallel do
    enddo

    tb%mu = mu

    call timer_stop('occu')

  endsubroutine occupy


  !**********************************************************************
  ! Returns: sum_i 2*f(e_i)-N.
  ! If mu is correct, it returns zero!
  !**********************************************************************
  function fsum(s, tb, T, mu, N)
    implicit none

    type(dense_solver_lapack_t), intent(in)  :: s
    type(dense_hamiltonian_t),   intent(in)  :: tb
    real(DP),                    intent(in)  :: T, mu, N

    real(DP)                                 :: fsum

    ! ---

    integer   :: k, i
    real(DP)  :: r

    ! ---

    r = 0

    do k = 1, tb%nk
       !$omp  parallel do default(none) &
       !$omp& shared(k, mu, s, T, tb) &
       !$omp& reduction(+:r)
       do i = 1, tb%norb
          r = r + 2*fermi_dirac(s%evals(i, k), mu, T)
       enddo
       !$omp end parallel do
    enddo

!    write (*, *)  r

    fsum = r - N

  endfunction fsum


  !>
  !!
  !! function SolveMu
  !
  !! returns the chemical potential mu for a electron system
  !! of N electrons (N/2 lowest energy states occupied on T=0)
  !! (temperature T) occupying energy states e(1:M). It uses
  !! a simple bisection method for solving the nonlinear
  !! equation for mu (e.g. Newton becomes unstable if T is very
  !! small.) Works also for exactly zero temperature.
  !>
  function SolveMu(s, tb, T, N, error) result(res)
    implicit none

    type(dense_solver_lapack_t), intent(in)   :: s
    type(dense_hamiltonian_t),   intent(in)   :: tb
    real(DP),                    intent(in)   :: T, N
    integer,           optional, intent(out)  :: error

    real(DP)                                  :: res

    ! ---

    integer  :: it
    real(DP) :: mu1, mu2, mu3, fmu1, fmu2, fmu3
    
    ! ---

    INIT_ERROR(error)

    !if( N/2>tb%norb ) stop 'noc must be wrong!'

    mu1  = minval(s%evals(1, :))
    mu2  = maxval(s%evals(tb%norb, :))

!    write (*, *)  s%evals(:, :)
!    write (*, *)  tb%norb, s%evals(1, 1)

!    write (*, *)  "mu1, mu2 = ", mu1, mu2

    fmu1 = fsum(s, tb, T, mu1, N)
    fmu2 = fsum(s, tb, T, mu2, N)

!    write (*, *)  "fmu1, fmu2 = ", fmu1, fmu2

    if( fmu1*fmu2 > 0d0 ) then
       RAISE_ERROR("Bisection algorithm could not find root. Did you specify the number of occupied orbitals? State: mu1 = " // mu1 // ", mu2 = " // mu2 // ", N(mu1)-N0 = " // fmu1 // ", N(mu2)-N0 = " // fmu2, error)
    end if

    it = 0
    do 
       it = it+1
       if (it > 10000) then
          RAISE_ERROR("More than 10000 iterations in trying to find the Fermi-level. Something is wrong here.", error)
       endif

       mu3  = 0.5d0*(mu1+mu2)

       fmu3 = fsum(s, tb, T, mu3, N)

!       write (*, *)  mu3, fmu3, N

       if (fmu3 == 0.0) then
          mu1 = mu3
          mu2 = mu3
          exit
       else if( fmu3*fmu1>0d0 ) then
          mu1  = mu3
          mu2  = mu2
          fmu1 = fmu3
       else
          mu1  = mu1
          mu2  = mu3
          fmu2 = fmu3
       end if
       if( abs(mu1-mu2)<1E-12 ) exit
    end do
    res = 0.5d0*(mu1+mu2)
  endfunction SolveMu


  !>
  !! Construct the density matrix
  !<
  subroutine construct_density_matrix(s, tb)
    implicit none

    type(dense_hamiltonian_t),   intent(inout)  :: tb
    type(dense_solver_lapack_t), intent(inout)  :: s

    ! ---

    integer   :: ia, jb, k, o
    WF_T(DP)  :: h1

    WF_T(DP), pointer  :: tb_rho(:, :, :)

    ! ---

    call timer_start('construct_density_matrix')

    call c_f_pointer(tb%rho, tb_rho, [tb%norb, tb%norb, tb%nk])

    call resize(tr_evecs, tb%norb, tb%norb)

    !
    ! Construct the density matrix rho_ll
    !

#define F s%f

    k_loop: do k = 1, tb%nk
       tr_evecs = transpose(s%evecs(:, :, k))

       i_loop: do ia = 1, tb%norb
          j_loop: do jb = 1, tb%norb
             h1 = 0
             do o = 1, tb%norb
                h1 = h1 + F(o, k)*tr_evecs(o, ia)*tr_evecs(o, jb)
             enddo

             tb_rho(ia, jb, k)  = h1

          enddo j_loop
       enddo i_loop
    enddo k_loop

    call timer_stop('construct_density_matrix')
    
  endsubroutine construct_density_matrix

endmodule dense_occupation
