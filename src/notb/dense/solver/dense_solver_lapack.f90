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
!! Determine the density matrix using the LAPACK generalized eigenvalue
!! solver.
!<

#include "macros.inc"
#include "filter.inc"

module dense_solver_lapack
  use, intrinsic :: iso_c_binding

  use supplib

  use particles

  use materials

  use dense_hamiltonian_type

  use dense_occupation

  implicit none

  private

  public :: ST_STANDARD, ST_DIVIDE_AND_CONQUER, ST_EXPERT, n_st, len_st_str
  public :: STR_standard, STR_divide_and_conquer, STR_expert, st_strs

  integer, parameter  :: ST_STANDARD            = 0
  integer, parameter  :: ST_DIVIDE_AND_CONQUER  = 1
  integer, parameter  :: ST_EXPERT              = 2
  integer, parameter  :: n_st                   = 3
  integer, parameter  :: len_st_str             = 25

  character(len_st_str), parameter  :: STR_standard            = CSTR("standard")
  character(len_st_str), parameter  :: STR_divide_and_conquer  = CSTR("divide-and-conquer")
  character(len_st_str), parameter  :: STR_expert              = CSTR("expert")
  character(len_st_str), parameter  :: st_strs(n_st) = &
       (/ STR_standard, STR_divide_and_conquer, STR_expert /)

  !
  ! Solver type
  !

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

#ifdef COMPLEX_WF
     WF_T(DP), allocatable :: work(:)
     real(DP), allocatable :: rwork(:)
     integer,  allocatable :: iwork(:)
#else
     WF_T(DP), allocatable :: work(:)
     integer,  allocatable :: iwork(:)
#endif

  endtype dense_solver_lapack_t

  !
  ! Interface definition
  !

  public :: init
  interface init
     module procedure dense_solver_lapack_init
  endinterface

  public :: del
  interface del
     module procedure dense_solver_lapack_del
  endinterface

  public :: diag_start
  interface diag_start
     module procedure dense_solver_lapack_diag_start
  endinterface

  public :: diag_stop
  interface diag_stop
     module procedure dense_solver_lapack_diag_stop
  endinterface

  public :: diag_HS
  interface diag_HS
     module procedure dense_solver_lapack_diag
  endinterface

  public :: e_bs
  interface e_bs
     module procedure dense_solver_lapack_e_bs
  endinterface

  public :: mulliken
  interface mulliken
     module procedure dense_solver_lapack_mulliken
  endinterface

  public :: get_dict
  interface get_dict
     module procedure dense_solver_lapack_get_dict
  endinterface

  public :: register
  interface register
     module procedure dense_solver_lapack_register
  endinterface

  !
  ! Private interface
  !

  interface set_Hamiltonian
     module procedure dense_solver_lapack_set_Hamiltonian
  endinterface

contains

  !**********************************************************************
  ! Initialize the solver
  !**********************************************************************
  subroutine dense_solver_lapack_init(this, solver_type, Tele, error)
    implicit none
    
    type(dense_solver_lapack_t), intent(inout)  :: this
    integer, intent(in), optional               :: solver_type
    real(DP), intent(in), optional              :: Tele
    integer, intent(out), optional              :: error

    ! ---

    INIT_ERROR(error)

    if (present(solver_type)) then
       this%solver_type  = solver_type
    endif
    if (present(Tele)) then
       this%Tele  = Tele
    endif

    call prlog("- dense_solver_lapack_init -")

    if (this%solver_type == ST_STANDARD) then
       call prlog("     Using standard driver")
    else if (this%solver_type == ST_DIVIDE_AND_CONQUER) then
       call prlog("     Using divide-and-conquer driver")
    else if (this%solver_type == ST_EXPERT) then
       call prlog("     Using expert driver")
    else
       RAISE_ERROR("Solver type " // this%solver_type // " is unknown.", error)
    endif

    call prlog

  endsubroutine dense_solver_lapack_init


  !**********************************************************************
  ! Delete the solver
  !**********************************************************************
  subroutine dense_solver_lapack_del(this)
    implicit none

    type(dense_solver_lapack_t), intent(inout)  :: this

    ! ---

    this%tb => NULL()

    if (allocated(this%f))      deallocate(this%f)
    if (allocated(this%evals))  deallocate(this%evals)
    if (allocated(this%evecs))  deallocate(this%evecs)

    if (allocated(this%work))   deallocate(this%work)
#ifdef COMPLEX_WF
    if (allocated(this%rwork))  deallocate(this%rwork)
#endif
    if (allocated(this%iwork))  deallocate(this%iwork)

  endsubroutine dense_solver_lapack_del


  !**********************************************************************
  ! Initialize the solver
  !**********************************************************************
  subroutine dense_solver_lapack_set_Hamiltonian(this, tb)
    implicit none
    
    type(dense_solver_lapack_t), intent(inout)  :: this
    type(dense_hamiltonian_t), target           :: tb

    ! ---

    call del(this)

    this%tb => tb
    this%norb = tb%norb

    allocate(this%f(tb%norb, tb%nk))
    allocate(this%evals(tb%norb, tb%nk))

    allocate(this%evecs(tb%norb, tb%norb, tb%nk))

  endsubroutine dense_solver_lapack_set_Hamiltonian


  !***************************************************************************
  ! Solve the generalized eigenvalue problem and construct the density matrix
  !***************************************************************************
  subroutine dense_solver_lapack_diag_start(this, tb, error)
    implicit none

    type(dense_solver_lapack_t),           intent(inout)  :: this
    type(dense_hamiltonian_t),   target                   :: tb
    integer,                     optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    if (.not. associated(this%tb, tb) .or. this%norb /= tb%norb) then
       call set_Hamiltonian(this, tb)
    endif

  endsubroutine dense_solver_lapack_diag_start


  !***************************************************************************
  ! Solve the generalized eigenvalue problem and construct the density matrix
  !***************************************************************************
  subroutine dense_solver_lapack_diag_stop(this, tb, error)
    implicit none

    type(dense_solver_lapack_t), intent(inout)        :: this
    type(dense_hamiltonian_t), target, intent(inout)  :: tb
    integer, intent(out), optional                    :: error

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)

  endsubroutine dense_solver_lapack_diag_stop


  !***************************************************************************
  ! Solve the generalized eigenvalue problem
  !***************************************************************************
  subroutine dense_solver_lapack_solve_HS( &
       this, tb, H, S, evals, evecs, phi, error)
    implicit none

    type(dense_solver_lapack_t), intent(inout)  :: this
    type(dense_hamiltonian_t), target           :: tb
    WF_T(DP), intent(in)                        :: H(tb%norb, tb%norb)
    WF_T(DP), intent(in)                        :: S(tb%norb, tb%norb)
    real(DP), intent(out)                       :: evals(tb%norb)
    WF_T(DP), intent(out)                       :: evecs(tb%norb, tb%norb)
    real(DP), intent(in), optional              :: phi(tb%nat)
    integer, intent(out), optional              :: error

    ! ---

    !
    ! Diagonalization workspace
    !

#ifdef COMPLEX_WF
    integer  :: lwork
    integer  :: lrwork
    integer  :: liwork
    WF_T(DP) :: opt_lwork
    real(DP) :: opt_lrwork
    integer  :: opt_liwork
#else
    integer  :: lwork
    integer  :: liwork
    WF_T(DP) :: opt_lwork
    integer  :: opt_liwork
#endif

    ! ---

    integer :: i, j, a, b, ia, jb, info

    real(DP), pointer :: S2(:, :)

    type(particles_t),    pointer :: tb_p
    type(notb_element_t), pointer :: tb_at(:)

    character(14) :: timer_str = "-------------"

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)
    call c_f_pointer(tb%p, tb_p)
    call c_f_pointer(tb%at, tb_at, [tb%nat])

    if (present(phi)) then

       !$omp  parallel do default(none) &
       !$omp& private(a, b, ia, j, jb) &
       !$omp& shared(evecs, H, phi, S, tb, this, tb_p, tb_at)
       do j = 1, tb%nat
          if (IS_EL(tb%f, tb_p, j)) then
             do b = 1, tb_at(j)%no
                do i = 1, tb%nat
                   if (IS_EL(tb%f, tb_p, i)) then
                      do a = 1, tb_at(i)%no
                         ia = tb_at(i)%o1 + a - 1
                         jb = tb_at(j)%o1 + b - 1

                         evecs(ia, jb) = H(ia, jb) &
                              - 0.5_DP*S(ia, jb)*(phi(i) + phi(j))
                      enddo
                   endif
                enddo
             enddo
          endif
       enddo

    else
       evecs(:, :)  = H(:, :)
    endif

    if (this%solver_type == ST_STANDARD) then
       timer_str = "LAPACK dsygv"
    else if (this%solver_type == ST_DIVIDE_AND_CONQUER) then
       timer_str = "LAPACK dsygvd"
    else if (this%solver_type == ST_EXPERT) then
       timer_str = "LAPACK dsygvx"
    else
       RAISE_ERROR("Solver type " // this%solver_type // " is unknown.", error)
    endif

    call timer_start(trim(timer_str))

!    if (this%solver_type == ST_EXPERT) then
!       if (.not. allocated(H)) then 
!          allocate(H(tb%norb, tb%norb))
!       endif
!       if (.not. allocated(ifail)) then
!          allocate(ifail(tb%norb))
!       endif
!    endif

    !
    ! Allocate diagonalization workspace
    !

    ! We abuse rho as a work buffer
    call c_f_pointer(tb%rho, S2, [tb%norb, tb%norb])

#ifdef COMPLEX_WF

    !      write (*, *)  "Allocated workspace."

    if (this%solver_type == ST_STANDARD) then

       call resize(this%rwork, max(1, 3*tb%norb-2))

       lwork = -1
       call zhegv( &
            1, 'V', 'U', tb%norb, &
            evecs(:, :), tb%norb, &
            this%   (:, :), tb%norb, &
            evals(:), &
            opt_lwork, this%lwork, this%rwork, info)

!        write (*, *)  info, opt_lwork

!          if (info /= 0)  stop "[solver_lapack_init] LAPACK workspace query failed."

       lwork = int(opt_lwork)
       call resize(this%work, lwork)

    else if (this%solver_type == ST_DIVIDE_AND_CONQUER) then

       lwork = -1
       lrwork = -1
       liwork = -1
       call zhegvd( &
            1, 'V', 'U', tb%norb, &
            evecs(:, :), tb%norb, &
            S2(:, :), tb%norb, &
            evals(:), &
            opt_lwork, lwork, opt_lrwork, lrwork, opt_liwork, liwork, info)

       if (info /= 0)  stop "[solver_lapack_init] LAPACK workspace query failed."

       lwork = int(opt_lwork)
       lrwork = int(opt_lrwork)
       liwork = int(opt_liwork)
       !            write (*, *)  "lwork, lrwork, liwork = ", lwork, lrwork, liwork
       call resize(this%work, lwork)
       call resize(this%rwork, lrwork)
       call resize(this%iwork, liwork)

!      else if (this%solver_type == ST_EXPERT) then
!
!          allocate(rwork(7*tb%norb))
!          allocate(iwork(5*tb%norb))
!
!          lwork = -1
!          call zhegvx( &
!               1, 'V', 'A', 'U', tb%norb, &
!               H(:, :), tb%norb, &
!               S2(:, :), tb%norb, &
!               0.0d0, 0.0d0, 0.0d0, 0.0d0, -1.0d0, &
!               nev, evals(:), &
!               evecs(:, :), tb%norb, &
!               opt_lwork, lwork, rwork, iwork, ifail, info)
!
!          if (info /= 0)  stop "[solver_lapack_init] LAPACK workspace query failed."
!
!          !            write (*, *)  "lwork = ", lwork
!          lwork = int(opt_lwork)
!          allocate(work(lwork))
!
    else
       RAISE_ERROR_AND_STOP_TIMER("Solver type " // this%solver_type // " is unknown.", trim(timer_str), error)
    endif

#else

    if (this%solver_type == ST_STANDARD) then

       lwork = -1
       call dsygv( &
            1, 'V', 'U', tb%norb, &
            evecs(:, :), tb%norb, &
            S2(:, :), tb%norb, &
            evals(:), &
            opt_lwork, lwork, info)

       if (info /= 0)  stop "[solver_lapack_init] LAPACK workspace query failed."            

       lwork = int(opt_lwork)
       !            write (*, *)  "lwork = ", lwork
       call resize(this%work, lwork)

    else if (this%solver_type == ST_DIVIDE_AND_CONQUER) then

       lwork = -1
       liwork = -1
       call dsygvd( &
            1, 'V', 'U', tb%norb, &
            evecs(:, :), tb%norb, &
            S2(:, :), tb%norb, &
            evals(:), &
            opt_lwork, lwork, opt_liwork, liwork, info)

       if (info /= 0)  stop "[solver_lapack_init] LAPACK workspace query failed."

       lwork = int(opt_lwork)
       liwork = int(opt_liwork)
       !            write (*, *)  "lwork, liwork = ", lwork, liwork
       call resize(this%work, lwork)
       call resize(this%iwork, liwork)

!       else if (this%solver_type == ST_EXPERT) then
!
!          allocate(iwork(5*tb%norb))
!
!          lwork = -1
!          call dsygvx( &
!               1, 'V', 'I', 'U', tb%norb, &
!               H(:, :), tb%norb, &
!               S2(:, :), tb%norb, &
!               0.0d0, 0.0d0, 1, this%n_bands, -1.0d0, &
!               nev, evals(:), &
!               evecs(:, :), tb%norb, &
!               opt_lwork, lwork, iwork, ifail, info)
!
!          if (info /= 0)  stop "[solver_lapack_init] LAPACK workspace query failed."
!
!          lwork = int(opt_lwork)
!          !            write (*, *)  "lwork = ", lwork
!          allocate(work(lwork))
!
    else
       RAISE_ERROR_AND_STOP_TIMER("Solver type " // this%solver_type // " is unknown.", trim(timer_str), error)
    endif

#endif

    ! We have to solve the eigenvalue problem for each k-point
    ! We abuse rho as a work buffer
    S2 = S(:, :)

#ifdef COMPLEX_WF

    if (this%solver_type == ST_STANDARD) then

       call zhegv(1, 'V', 'L', &
            tb%norb, evecs(:, :), &
            tb%norb, S2, &
            tb%norb, evals(:), &
            this%work, lwork, this%rwork, info)

    else if (this%solver_type == ST_DIVIDE_AND_CONQUER) then

       call zhegvd(1, 'V', 'L', &
            tb%norb, evecs(:, :), &
            tb%norb, S2, &
            tb%norb, evals(:), &
            this%work, lwork, this%rwork, lrwork, this%iwork, liwork, info)

!    else if (this%solver_type == ST_EXPERT) then
!
!       H = evecs(:, :)
!
!       call zhegvx(1, 'V', 'A', 'L', &
!            tb%norb, &
!            H, tb%norb, &
!            S2, tb%norb, &
!            0.0d0, 0.0d0, 0, 0, -1.0d0, &
!            n_evecs, evals(:), &
!            evecs(:, :), tb%norb, &
!            work, lwork, rwork, iwork, ifail, info)

    else
       RAISE_ERROR_AND_STOP_TIMER("Solver type " // this%solver_type // " is unknown.", trim(timer_str), error)
    endif

#else
    if (this%solver_type == ST_STANDARD) then

       call dsygv(1, 'V', 'L', &
            tb%norb, evecs(:, :), &
            tb%norb, S2, &
            tb%norb, evals(:), &
            this%work, lwork, info)

    else if (this%solver_type == ST_DIVIDE_AND_CONQUER) then

       call dsygvd(1, 'V', 'L', &
            tb%norb, evecs(:, :), &
            tb%norb, S2, &
            tb%norb, evals(:), &
            this%work, lwork, this%iwork, liwork, info)

!    else if (this%solver_type == ST_EXPERT) then
!
!       H = evecs(:, :)
!
!       call dsygvx(1, 'V', 'I', 'L', &
!            tb%norb, &
!            H, tb%norb, &
!            S2, tb%norb, &
!            0.0d0, 0.0d0, 1, this%n_bands, -1.0d0, &
!            nev, evals(:), &
!            evecs(:, :), tb%norb, &
!            work, lwork, iwork, ifail, info)

    else
       RAISE_ERROR_AND_STOP_TIMER("Solver type " // this%solver_type // " is unknown.", trim(timer_str), error)
    endif

#endif

    if (info /= 0) then
!       if (this%solver_type == ST_EXPERT) then
!
!          if (info <= tb%norb) then
!             write (ilog, '(A)')  "The following eigenvectors did not converge:"
!
!             do i = 1, info
!                write (ilog, '(I10)')  ifail(i)
!             enddo
!
!          endif
!
!       endif

       RAISE_ERROR_AND_STOP_TIMER("Diagonalization failed with error code "//info//".", trim(timer_str), error)
    endif

!    if (this%solver_type == ST_EXPERT) then
!
!       if (nev /= this%n_bands)  then
!          write (*, '(A,I5,I5)')  "[diag] Fatal: Number of eigenvalues not equal number of bands: ", nev, this%n_bands
!          stop
!       else if (this%n_bands > tb%norb) then
!          evals(this%n_bands+1:tb%norb) = evals(this%n_bands) &
!               + 100*(evals(this%n_bands)-evals(1, k))
!       endif
!
!    endif

    call timer_stop(trim(timer_str))

  endsubroutine dense_solver_lapack_solve_HS



  !***************************************************************************
  ! Solve the generalized eigenvalue problem and construct the density matrix
  !***************************************************************************
  subroutine dense_solver_lapack_diag(this, tb, noc, phi, error)
    implicit none

    type(dense_solver_lapack_t), intent(inout)  :: this
    type(dense_hamiltonian_t),   target         :: tb
    real(DP),                    intent(in)     :: noc
    real(DP),          optional, intent(in)     :: phi(tb%nat)
    integer,           optional, intent(out)    :: error

    ! ---

    integer  :: k

    type(particles_t), pointer  :: tb_p
    WF_T(DP),          pointer  :: tb_H(:, :, :), tb_S(:, :, :)

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)
    call c_f_pointer(tb%p, tb_p)
    call c_f_pointer(tb%H, tb_H, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%S, tb_S, [tb%norb, tb%norb, tb%nk])

    call timer_start("solver_lapack_diag")

    !
    ! We have to solve the eigenvalue problem for each k-point
    !

    if (present(phi)) then

       do k = 1, tb%nk
          call dense_solver_lapack_solve_HS( &
               this, tb, &
               tb_H(:, :, k), tb_S(:, :, k), &
               this%evals(:, k), this%evecs(:, :, k), &
               phi(:), error=error)
          PASS_ERROR_AND_STOP_TIMER("solver_lapack_diag", error)
       enddo

    else

       do k = 1, tb%nk
          call dense_solver_lapack_solve_HS( &
               this, tb, &
               tb_H(:, :, k), tb_S(:, :, k), &
               this%evals(:, k), this%evecs(:, :, k), &
               error=error)
          PASS_ERROR_AND_STOP_TIMER("solver_lapack_diag", error)
       enddo

    endif

    call occupy(tb, this%evals, noc, this%Tele, this%f, error=error)
    PASS_ERROR_AND_STOP_TIMER("solver_lapack_diag", error)

    if (present(phi)) then

       call construct_density_and_energy_matrix(this, tb, tb_p, phi)

    else

       call construct_density_and_energy_matrix(this, tb, tb_p)

    endif

    call timer_stop("solver_lapack_diag")

  endsubroutine dense_solver_lapack_diag



  !***************************************************************************
  !
  !                  subroutine E_BS
  !
  ! calculates the electronic structure part of the total energy, i.e. the 
  ! first two terms in equation (19) in Elstner et.al. PRB 58, 7260 (1998)
  !
  !***************************************************************************
  function dense_solver_lapack_e_bs(this, tb, error) result(Ebs)
    implicit none

    type(dense_solver_lapack_t),  intent(in)  :: this
    type(dense_hamiltonian_t), target         :: tb
    integer, intent(out), optional            :: error
    real(DP)                                  :: Ebs

    ! ---

    integer   :: ia, jb, k
    WF_T(DP)  :: Ebs_c !, Ebs_h

    WF_T(DP), pointer  :: tb_H(:, :, :), tb_rho(:, :, :)
    
    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)
    call c_f_pointer(tb%H, tb_H, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%rho, tb_rho, [tb%norb, tb%norb, tb%nk])

    call timer_start('solver_lapack_e_bs')

    Ebs_c = 0.0_DP
!    Ebs_h = 0d0
    do k = 1, tb%nk
!       Ebs_h = Ebs_h + sum(this%f(:, k)*this%evals(:, k))
       !$omp  parallel do default(none) &
       !$omp& private(ia) &
       !$omp& shared(k, tb, this, tb_rho, tb_H) &
       !$omp& reduction(+:Ebs_c)
       do jb = 1, tb%norb
          do ia = 1, tb%norb
             Ebs_c = Ebs_c + tb_rho(ia, jb, k) * tb_H(ia, jb, k)
          enddo
       enddo
    enddo

#ifdef COMPLEX_WF
    Ebs = real(Ebs_c)    
#else
    Ebs = Ebs_c
#endif COMPLEX_WF

    call timer_stop('solver_lapack_e_bs')

  endfunction dense_solver_lapack_e_bs


  subroutine dense_solver_lapack_mulliken(this, tb, q_out, error)
    implicit none

    type(dense_solver_lapack_t),           intent(inout)  :: this
    type(dense_hamiltonian_t),   target                   :: tb
    real(DP),                              intent(out)    :: q_out(:)
    integer,                     optional, intent(out)    :: error

    ! ---

    integer   :: I,J,alpha,beta,Ia,Jb,kp
    WF_T(DP)  :: q(tb%nat) 

    type(particles_t),    pointer  :: tb_p
    type(notb_element_t), pointer  :: tb_at(:)
    WF_T(DP),             pointer  :: tb_S(:, :, :), tb_rho(:, :, :)
    real(DP),             pointer  :: tb_n(:)

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)
    call c_f_pointer(tb%p, tb_p)
    call c_f_pointer(tb%at, tb_at, [tb%nat])
    call c_f_pointer(tb%S, tb_S, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%rho, tb_rho, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%n, tb_n, [tb%nat])

    call timer_start('solver_lapack_mulliken')

    q(:)  = 0.0_DP

    do kp = 1, tb%nk
       !$omp  parallel do default(none) &
       !$omp& private(alpha, ia, j, jb) &
       !$omp& shared(kp, tb, this, tb_p, tb_at, tb_rho, tb_S) &
       !$omp& reduction(+:q)
       i_loop: do I = 1, tb%nat
          if (IS_EL(tb%f, tb_p, I)) then
             do alpha = 0, tb_at(I)%no - 1
                Ia = tb_at(I)%o1 + alpha
                j_loop: do J = 1, tb%nat
                   if (IS_EL(tb%f, tb_p, J)) then
                      do beta = 0, tb_at(J)%no - 1
                         Jb = tb_at(J)%o1 + beta

                         q(I) = q(I) + ( - tb_rho(Ia, Jb, kp) * tb_S(Jb, Ia, kp) )
                      enddo
                   endif
                enddo j_loop
             enddo
          endif
       enddo i_loop
    enddo

#ifdef COMPLEX_WF
    tb_n = real(q)
#else
    tb_n = q
#endif COMPLEX_WF

!    write (*, *)  sum(q(:))

    !$omp  parallel do default(none) &
    !$omp& shared(q_out, tb, this, tb_p, tb_n, tb_at)
    do i = 1, tb%nat
       if (IS_EL(tb%f, tb_p, i)) then
          q_out(i) = (tb_n(i) + tb_at(i)%q0)
       endif
    enddo

    call timer_stop('solver_lapack_mulliken')

  endsubroutine dense_solver_lapack_mulliken


  !**********************************************************************
  ! Construct the density and energy matrix
  !**********************************************************************
  subroutine construct_density_and_energy_matrix(s, tb, p, phi)
    implicit none

    type(dense_hamiltonian_t), intent(inout)  :: tb
    type(dense_solver_lapack_t), intent(inout)  :: s
    type(particles_t), intent(in)         :: p
    real(DP), optional, intent(in)        :: phi(p%nat)

    ! ---

    integer   :: ia1, ia2, jb, i, j, a, b, k
    WF_T(DP)  :: ec1, ec
    WF_T(DP), allocatable, save  :: rho(:, :), e(:, :)

    type(notb_element_t), pointer  :: tb_at(:)
    WF_T(DP),             pointer  :: tb_rho(:, :, :), tb_e(:, :, :)

    ! ---

    call timer_start('construct_density_and_energy_matrix')
    call c_f_pointer(tb%at, tb_at, [tb%nat])
    call c_f_pointer(tb%rho, tb_rho, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%e, tb_e, [tb%norb, tb%norb, tb%nk])
   
    call resize(tr_evecs, tb%norb, tb%norb)
    call resize(rho, tb%norb, tb%norb)
    call resize(e, tb%norb, tb%norb)

    !
    ! Construct the density matrix rho_ll and H_rl * rho_ll (for use in the
    ! forces)
    !

#define F s%f

    k_loop: do k = 1, tb%nk
       tr_evecs(:, :)   = transpose(s%evecs(:, :, k))

       !$omp  parallel do default(none) &
       !$omp& shared(e, k, p, rho, s, tb, tr_evecs, tb_at) &
       !$omp& private(a, ia1)
       do i = 1, p%natloc
          if (IS_EL(tb%f, p, i)) then
             do a = 1, tb_at(i)%no
                ia1 = tb_at(i)%o1 + a - 1

                rho(:, ia1) = F(:, k)*tr_evecs(:, ia1)
                e(:, ia1)   = s%evals(:, k)*rho(:, ia1)

             enddo
          endif
       enddo

       call GEMM('T', 'N', &
            tb%norb, tb%norb, tb%norb, &
            1.0_DP, &
            tr_evecs(:, :), tb%norb, &
            rho(:, :), tb%norb, &
            0.0_DP, &
            tb_rho(:, :, k), tb%norb)

       call GEMM('T', 'N', &
            tb%norb, tb%norb, tb%norb, &
            1.0_DP, &
            tr_evecs(:, :), tb%norb, &
            e(:, :), tb%norb, &
            0.0_DP, &
            tb_e(:, :, k), tb%norb)


       if (present(phi)) then

          !$omp  parallel do default(none) &
          !$omp& shared(k, p, phi, tb, tb_at, tb_e, tb_rho) &
          !$omp& private(ec, ec1, ia1, ia2, jb)
          do j = 1, p%nat
             if (IS_EL(tb%f, p, j)) then
                ec1 = -0.5 * phi(j)
             
                do b = 1, tb_at(j)%no
                   jb = tb_at(j)%o1 + b - 1

                   do i = 1, p%nat
                      if (IS_EL(tb%f, p, i)) then
                         ec = ec1 - 0.5 * phi(i)
                         ia1 = tb_at(i)%o1
                         ia2 = tb_at(i)%o1+tb_at(i)%no-1
                         tb_e(ia1:ia2, jb, k) = tb_e(ia1:ia2, jb, k) - tb_rho(ia1:ia2, jb, k)*ec
                      endif
                   enddo
                enddo
             endif
          enddo

       endif

    enddo k_loop

    call timer_stop('construct_density_and_energy_matrix')
    
  endsubroutine construct_density_and_energy_matrix


  !>
  !! Return dictionary object containing pointers to internal data
  !<
  subroutine dense_solver_lapack_get_dict(this, dict, error)
    implicit none

    type(dense_solver_lapack_t), intent(inout) :: this        !< NOTB object
    type(ptrdict_t),             intent(inout) :: dict
    integer,           optional, intent(out)   :: error       !< Error signals

    ! ---

    integer :: nk

    ! ---

    INIT_ERROR(error)

    nk = size(this%evals, 2)
    if (nk == 1) then
       if (allocated(this%evals)) then
          call register(dict, this%evals(:, 1), "eigenvalues")
       endif
       if (allocated(this%evecs)) then
          call register(dict, this%evecs(:, :, 1), "eigenvectors")
       endif
       if (allocated(this%f)) then
          call register(dict, this%f(:, 1), "occupation")
       endif
    else
       if (allocated(this%evals)) then
          call register(dict, this%evals, "eigenvalues")
       endif
       if (allocated(this%evecs)) then
          call register(dict, this%evecs, "eigenvectors")
       endif
       if (allocated(this%f)) then
          call register(dict, this%f, "occupation")
       endif
    endif

  endsubroutine dense_solver_lapack_get_dict


  !>
  !! Register this solver object
  !<
  subroutine dense_solver_lapack_register(this, cfg)
    implicit none

    type(dense_solver_lapack_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                        :: cfg

    ! ---

    type(c_ptr)  :: cfg2
    
    ! ---

    cfg2 = ptrdict_register_module(cfg, c_loc(this%enabled), &
         CSTR("SolverLAPACK"), &
         CSTR("Use the standard LAPACK routines for diagonalization and determination of the density matrix."))

    call ptrdict_register_enum_property(cfg2, c_loc(this%solver_type), &
         n_st, len_st_str, st_strs, &
         CSTR("solver_type"), &
         CSTR("Solver to use: LAPACK 'standard' or 'divide-and-conquer'"))

    call ptrdict_register_real_property(cfg2, c_loc(this%Tele), &
         CSTR("electronic_T"), &
         CSTR("Electronic temperature"))

  endsubroutine dense_solver_lapack_register

endmodule dense_solver_lapack
