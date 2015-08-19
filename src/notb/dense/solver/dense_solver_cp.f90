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
!! Canonical purification method
!!
!! See:
!! A. H. R. Palser and D. E. Manolopulos, Phys. Rev. B 58, 12704 (1998)
!<

#include "macros.inc"
#include "filter.inc"

module dense_solver_cp
  use, intrinsic :: iso_c_binding

  use supplib

  use particles
  use materials
  use dense_hamiltonian_type

  implicit none

  private

  character(*), parameter  :: MODULE_STR  = "SolverCP"

  public :: dense_solver_cp_t
  type dense_solver_cp_t

     logical(BOOL)  :: enabled  = .false.

     integer   :: it
     integer   :: it_diag

     !
     ! Tight-binding object
     !

     type(dense_hamiltonian_t), pointer  :: tb  => NULL()

     !
     ! Trace convergence?
     !

     logical(BOOL)  :: trace  = .false.

     !
     ! Solver stuff
     !

     real(DP)               :: emin, emax       ! Eigenvalue bounds (of H_rl)
     real(DP)               :: E_BS             ! Bandstructure energy

     WF_T(DP), allocatable  :: S_rr(:, :, :)    ! the inverse of the overlap matrix

     WF_T(DP), allocatable  :: H0_rl(:, :, :)   ! the (unshifted) Hamiltonian matrix (contra-covariant)
     WF_T(DP), allocatable  :: H_rl(:, :, :)    ! the Hamiltonian matrix
     WF_T(DP), allocatable  :: rho_rl(:, :, :)  ! the density matrix (contra-covariant)

     !
     ! Convergence criteria
     !

     real(DP)           :: epsilon = 1d-6

     !
     ! Template helper matrices
     !

     WF_T(DP), allocatable  :: help1(:, :)
     WF_T(DP), allocatable  :: help2(:, :)

  endtype dense_solver_cp_t


  !
  ! Interface definition
  !

  public :: init
  interface init
     module procedure dense_solver_cp_init
  endinterface

  public :: del
  interface del
     module procedure dense_solver_cp_del
  endinterface

  public :: diag_start
  interface diag_start
     module procedure dense_solver_cp_diag_start
  endinterface

  public :: diag_stop
  interface diag_stop
     module procedure dense_solver_cp_diag_stop
  endinterface

  public :: diag_HS
  interface diag_HS
     module procedure dense_solver_cp_diag
  endinterface

  public :: e_bs
  interface e_bs
     module procedure dense_solver_cp_e_bs
  endinterface

  public :: mulliken
  interface mulliken
     module procedure dense_solver_cp_mulliken
  endinterface

  public :: register
  interface register
     module procedure dense_solver_cp_register
  endinterface

  !
  ! Private interface
  !

  interface set_Hamiltonian
     module procedure dense_solver_cp_set_Hamiltonian
  endinterface

contains

  !**********************************************************************
  ! Initialize the solver
  !**********************************************************************
  subroutine dense_solver_cp_init(this, epsilon, error)
    implicit none
    
    type(dense_solver_cp_t), intent(inout)  :: this
    real(DP), intent(in), optional          :: epsilon
    integer, intent(out), optional          :: error


    ! ---

    INIT_ERROR(error)

    if (present(epsilon)) then
       this%epsilon  = epsilon
    endif

    this%it       = 0
    this%it_diag  = 0

  endsubroutine dense_solver_cp_init


  !**********************************************************************
  ! Delete the solver
  !**********************************************************************
  subroutine dense_solver_cp_del(this)
    implicit none

    type(dense_solver_cp_t), intent(inout)  :: this

    ! ---

    this%tb => NULL()

    if (allocated(this%S_rr))    deallocate(this%S_rr)

    if (allocated(this%H0_rl))   deallocate(this%H0_rl)
    if (allocated(this%H_rl))    deallocate(this%H_rl)
    if (allocated(this%rho_rl))  deallocate(this%rho_rl)

    if (allocated(this%help1))   deallocate(this%help1)
    if (allocated(this%help2))   deallocate(this%help2)

  endsubroutine dense_solver_cp_del


  !**********************************************************************
  ! Set the tight-binding object
  !**********************************************************************
  subroutine dense_solver_cp_set_Hamiltonian(this, tb)
    implicit none
    
    type(dense_solver_cp_t), intent(inout)            :: this
    type(dense_hamiltonian_t), target, intent(inout)  :: tb

    ! ---

    call del(this)

    this%tb  => tb

    allocate(this%S_rr(tb%norb, tb%norb, tb%nk))

    allocate(this%H0_rl(tb%norb, tb%norb, tb%nk))
    allocate(this%H_rl(tb%norb, tb%norb, tb%nk))
    allocate(this%rho_rl(tb%norb, tb%norb, tb%nk))

    allocate(this%help1(tb%norb, tb%norb))
    allocate(this%help2(tb%norb, tb%norb))

  endsubroutine dense_solver_cp_set_Hamiltonian


  !**********************************************************************
  ! Calculate the inverse of the overlap matrix
  !**********************************************************************
  subroutine dense_solver_cp_diag_start(this, tb, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(dense_solver_cp_t), intent(inout)  :: this
    type(dense_hamiltonian_t), target       :: tb
    integer, intent(out), optional          :: error

    ! ---

    integer :: k
#ifdef NO_BIND_C_OPTIONAL
    integer :: nit
#endif

    WF_T(DP), pointer :: tb_H(:, :, :), tb_S(:, :, :)
    logical(C_BOOL) :: l

    ! ---

    INIT_ERROR(error)

    if (.not. associated(this%tb, tb)) then
       call set_Hamiltonian(this, tb)
    endif

    call c_f_pointer(tb%S, tb_S, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%H, tb_H, [tb%norb, tb%norb, tb%nk])

    call timer_start("dense_solver_cp_diag_start")

    !
    ! Solve for S^-1
    !

    do k = 1, tb%nk
     
       l = this%it > 0
       call iterative_matrix_inverse( &
            tb_S(:, :, k), &
            this%S_rr(:, :, k), &
            tb%norb, &
            l, &
            this%epsilon, &
#ifdef NO_BIND_C_OPTIONAL
            cublas_handle = C_NULL_PTR, &
            nit = nit, &
#endif
            work1 = this%help1, &
            work2 = this%help2, &
            error = error)
       PASS_ERROR(error)
       call MM(tb%norb, 1.0d0, this%S_rr(:, :, k), tb_H(:, :, k), 0.0d0, this%H0_rl(:, :, k))

       this%H0_rl(:, :, k)  = transpose(this%H0_rl(:, :, k))

    enddo

    call timer_stop("dense_solver_cp_diag_start")

  endsubroutine dense_solver_cp_diag_start


  !**********************************************************************
  ! Finalize
  !**********************************************************************
  subroutine dense_solver_cp_diag_stop(this, tb, error)
    implicit none

    type(dense_solver_cp_t), intent(inout)  :: this
    type(dense_hamiltonian_t), target       :: tb
    integer, intent(out), optional          :: error

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)

    this%it = this%it + 1

  endsubroutine dense_solver_cp_diag_stop


  !**********************************************************************
  ! Calculate the density matrix using the CP method
  !**********************************************************************
  subroutine dense_solver_cp_diag(this, tb, N0, phi, error)
!#ifdef MKL
!    use mkl95_blas
!#endif

    implicit none

    type(dense_solver_cp_t), intent(inout)  :: this
    type(dense_hamiltonian_t), target       :: tb
    real(DP), intent(in)                    :: N0    ! The number of electrons in the system
    real(DP), intent(in), optional          :: phi(tb%nat)
    integer, intent(out), optional          :: error

    ! ---

    real(DP)       :: emin, emax

    real(DP)       :: lambda, mu, E, E_old, c, tr_r, tr_r2, tr_r3, ec

    integer        :: i, j, a, b, ia, jb, k, nit, un

    character(100) :: fn

    type(particles_t),    pointer  :: tb_p
    type(notb_element_t), pointer  :: tb_at(:)
    WF_T(DP),             pointer  :: tb_H(:, :, :), tb_S(:, :, :)
    WF_T(DP),             pointer  :: tb_rho(:, :, :), tb_e(:, :, :)

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)
    call c_f_pointer(tb%p, tb_p)
    call c_f_pointer(tb%at, tb_at, [tb%nat])
    call c_f_pointer(tb%H, tb_H, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%S, tb_S, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%rho, tb_rho, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%e, tb_e, [tb%norb, tb%norb, tb%nk])

    call timer_start("dense_solver_cp_diag")

    !
    ! Construct the shifted Hamiltonian matrix
    !

    do k = 1, tb%nk

       if (present(phi)) then

          do j = 1, tb%nat
             if (IS_EL(tb%f, tb_p, j)) then
                do b = 1, tb_at(j)%no
                   do i = 1, tb%nat
                      if (IS_EL(tb%f, tb_p, i)) then
                         do a = 1, tb_at(i)%no
                            ia = tb_at(i)%o1 + a - 1
                            jb = tb_at(j)%o1 + b - 1

                            this%help1(ia, jb) = tb_H(ia, jb, k) &
                                 - 0.5_DP*tb_S(ia, jb, k)*(phi(i) + phi(j))

                         enddo
                      endif
                   enddo
                enddo
             endif
          enddo

          call MM(tb%norb, 1.0d0, this%S_rr(:, :, k), this%help1, 0.0d0, this%H_rl(:, :, k))

       else

          this%H_rl(:, :, k)  = transpose(this%H0_rl(:, :, k))

       endif

    enddo

    nit        = 0
    this%E_BS  = 0.0_DP
    do k = 1, tb%nk

       !
       ! Guess a lower and an upper bound for the eigenvalue spectrum.
       !

       call ev_bounds(tb%norb, this%H_rl(:, :, k), emin, emax)

       mu      = tr(tb%norb, this%H_rl(:, :, k))/tb%norb
       lambda  = min(N0/(emax-mu), (tb%norb-N0)/(mu-emin))/tb%norb

       this%rho_rl(:, :, k) = -lambda*this%H_rl(:, :, k)
       do i = 1, tb%norb
          this%rho_rl(i, i, k) = this%rho_rl(i, i, k)+lambda*mu+real(N0, DP)/tb%norb
       enddo

!       E     = multr(tb%norb, this%rho_rl(:, :, k), this%H0_rl(:, :, k))
       E     = sum(this%rho_rl(:, :, k)*this%H0_rl(:, :, k))
       E_old = E+1e6

       c = 0.5

       this%it_diag  = this%it_diag + 1

       if (this%trace) then
          write (fn, '(A,I6.6,A)')  "dense_solver_cp_convergence_", this%it_diag, ".out"
          un = fopen(fn, F_WRITE)
          write (un, '(ES20.10)')  E
       endif

       do while (abs(E - E_old) > tb%norb*this%epsilon .and. (c >= 0 .and. c <= 1))

          tr_r = tr(tb%norb, this%rho_rl(:, :, k))

          call MM(tb%norb, 1.0d0, this%rho_rl(:, :, k), this%rho_rl(:, :, k), 0.0d0, this%help1)

          tr_r2 = tr(tb%norb, this%help1)

          call MM(tb%norb, 1.0d0, this%rho_rl(:, :, k), this%help1, 0.0d0, this%help2)

          tr_r3 = tr(tb%norb, this%help2)

          E_old = E
          if (tr_r /= tr_r2) then
             c = (tr_r2-tr_r3)/(tr_r-tr_r2)

             if (c >= 0 .and. c <= 1) then
                if (c < 0.5) then
                   this%rho_rl(:, :, k) = ((1-2*c)*this%rho_rl(:, :, k)+(1+c)*this%help1-this%help2)/(1-c)
                else
                   this%rho_rl(:, :, k) = ((1+c)*this%help1-this%help2)/c
                endif

!                E     = multr(tb%norb, this%rho_rl(:, :, k), this%H0_rl(:, :, k))
                E     = sum(this%rho_rl(:, :, k)*this%H0_rl(:, :, k))
             endif
          endif

          if (this%trace) then
             write (un, '(ES20.10)')  E
          endif

          nit = nit+1

          if (mod(nit, 100) == 0) then
             WARN("No convergence after " // nit // " iterations (density matrix).")
          endif

       enddo

       if (this%trace) then
          call fclose(un)
       endif

       this%E_BS = this%E_BS + 2*E
    enddo

    !
    ! Now calculate the lowered density matrix
    !

    do k = 1, tb%nk
       call MM(tb%norb, 2.0d0, this%rho_rl(:, :, k), this%S_rr(:, :, k), 0.0d0, tb_rho(:, :, k))

       call MM(tb%norb, 1.0d0, this%H_rl(:, :, k), tb_rho(:, :, k), 0.0d0, tb_e(:, :, k))

       if (present(phi)) then
       
          i_loop: do i = 1, tb%nat
             if (IS_EL(tb%f, tb_p, i)) then
                j_loop: do j = 1, tb%nat
                   if (IS_EL(tb%f, tb_p, j)) then
                      ec = -0.5 * ( phi(i) + phi(j) )

                      a_loop: do a = 1, tb_at(i)%no
                         b_loop: do b = 1, tb_at(j)%no
                            ia = tb_at(i)%o1 + a - 1
                            jb = tb_at(j)%o1 + b - 1
       
                            tb_e(ia, jb, k) = tb_e(ia, jb, k) - tb_rho(ia, jb, k)*ec
                         enddo b_loop
                      enddo a_loop
                   endif
                enddo j_loop
             endif
          enddo i_loop

       endif

    enddo

    tb%mu      = mu
    this%emin  = emin
    this%emax  = emax

    call timer_stop("dense_solver_cp_diag")
    
  endsubroutine dense_solver_cp_diag


  !***************************************************************************
  ! Mulliken charge analysis
  !***************************************************************************
  subroutine dense_solver_cp_mulliken(this, tb, q, error)
    implicit none

    type(dense_solver_cp_t), intent(inout)  :: this
    type(dense_hamiltonian_t), target       :: tb
    real(DP), intent(out)                   :: q(:)
    integer, intent(out), optional          :: error

    ! ---

    integer   :: k, i, a, ia

    type(particles_t),    pointer  :: tb_p
    type(notb_element_t), pointer  :: tb_at(:)
    real(DP),             pointer  :: tb_n(:)

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)
    call c_f_pointer(tb%p, tb_p)
    call c_f_pointer(tb%at, tb_at, [tb%nat])
    call c_f_pointer(tb%n, tb_n, [tb%nat])

    call timer_start('dense_solver_cp_mulliken')

    tb_n  = 0.0_DP

    do k = 1, tb%nk
       do i = 1, tb%nat
          if (IS_EL(tb%f, tb_p, i)) then
             do a = 0, tb_at(i)%no - 1
                ia = tb_at(i)%o1 + a
                tb_n(i) = tb_n(i) - 2*this%rho_rl(ia, ia, k)
             enddo
          endif
       enddo
    enddo

    do i = 1, tb%nat
       if (IS_EL(tb%f, tb_p, i)) then
          q(i) = (tb_n(i) + tb_at(i)%q0)
       endif
    enddo

    call timer_stop('dense_solver_cp_mulliken')

  endsubroutine dense_solver_cp_mulliken


  !**********************************************************************
  ! Calculate the band-structure energy
  !**********************************************************************
  function dense_solver_cp_e_bs(this, tb, error)
    implicit none

    type(dense_solver_cp_t), intent(in)  :: this
    type(dense_hamiltonian_t), target    :: tb
    integer, intent(out), optional       :: error
    real(DP)                             :: dense_solver_cp_e_bs

    ! ---

    INIT_ERROR(error)
    ASSERT_ASSOCIATION(this%tb, tb, error)

    dense_solver_cp_e_bs = this%E_BS

  endfunction dense_solver_cp_e_bs


  !>
  !! Register the solver object
  !<
  subroutine dense_solver_cp_register(this, cfg)
    implicit none

    type(dense_solver_cp_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                    :: cfg

    ! ---

    type(c_ptr)  :: cfg2
    
    ! ---

    cfg2 = ptrdict_register_module(cfg, c_loc(this%enabled), CSTR("SolverCP"), &
         CSTR("Canonical purification solver."))

    call ptrdict_register_real_property(cfg2, c_loc(this%epsilon), &
         CSTR("epsilon"), &
         CSTR("Convergence criterion (change in energy)"))

    call ptrdict_register_boolean_property(cfg2, c_loc(this%trace), &
         CSTR("trace"), &
         CSTR("Convergence criterion (change in energy)"))

  endsubroutine dense_solver_cp_register

endmodule dense_solver_cp
