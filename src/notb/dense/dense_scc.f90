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

#include "macros.inc"

!#define DEBUG_ENERGY_MINIMUM


!>
!! Charge self-consistant tight-binding
!!
!! Contains the self-consistent looping. Calculate charges, form new Hamiltonian, diagonalize, occupy, calculate
!! charges, check convergence ... Uses Anderson mixing to speed up 
!! convergence.
!!
!! To adjust there are 4 parameters:
!!  - mixing parameter beta
!!  - Mmax (~5), memory used in Anderson mixing
!!  - itmax; maximum number of iterations (in one cycle; when one cycle
!!    is done and no convergence, make a new cycle with itmax = 2*itmax
!!    and beta(new) = beta(current) / 2. 
!!  - itmax2: Maximum number itmax
!!
!! See
!! Elstner et al., Phys. Rev. B 58, 7260 (1998)
!<
module dense_scc
  use, intrinsic :: iso_c_binding

  use supplib

  use anderson_mixer

  use particles
  use neighbors
  use filter

  use coulomb

  use materials

  use dense_hamiltonian_type
  use dense_hamiltonian
  use dense_solver

#ifdef DEBUG_ENERGY_MINIMUM
  use rng
#endif

  implicit none

  private

  character(*), parameter, private  :: MODULE_STR = "SCC"

  public :: dense_scc_t
  type dense_scc_t

     logical(BOOL)  :: enabled = .false.  !< For standalone code

     !
     ! SCC stuff
     !

     integer       :: max_nit         = 200         !< max nr of iterations in self-consistency
     real(DP)      :: dq_crit         = 0.0001_DP   !< limit for convergence = max(dq_in-qd_out) 
     real(DP)      :: beta            = 0.2_DP      !< density mixing parameter in iteration
     integer       :: andersen_memory = 3           !< M in Anderson mixing in iteration

     integer       :: warn            = 20          !< warn after 20 iterations
     logical(BOOL) :: log             = .false.     !< write a status report for each SCC step

     logical       :: charges_only    = .false.     !< Only calculate Mulliken charges

     real(DP), allocatable  :: phi(:)

     !
     ! Position and charge history
     !

     integer               :: extrapolation_memory = 3  !< Number of past steps to keep
     integer               :: history_counter = 0
     real(DP), allocatable :: r(:, :, :)
     real(DP), allocatable :: q(:, :)

     !
     ! Statistics/diagnostic
     !

     integer :: niterations = 0
     integer :: nsteps      = 0
     integer :: nfail       = 0

     !
     ! Configuration objects
     !

     type(dense_solver_t), pointer :: solv  => NULL()
     type(C_PTR)                   :: coul  = C_NULL_PTR

     !
     ! Associated objects, control internal buffer size
     !

     type(particles_t), pointer          :: p     => NULL()
     type(dense_hamiltonian_t), pointer  :: tb    => NULL()

  endtype dense_scc_t


  !
  ! Interface
  !

  public :: init
  interface init
     module procedure dense_scc_init
  endinterface

  public :: set
  interface set
     module procedure dense_scc_set
  endinterface

  public :: del
  interface del
     module procedure dense_scc_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure dense_scc_bind_to
  endinterface

  public :: set_Coulomb
  interface set_Coulomb
     module procedure dense_scc_set_Coulomb
  endinterface set_Coulomb

  public :: set_solver
  interface set_solver
     module procedure dense_scc_set_solver
  endinterface set_solver

  public :: establish_self_consistency
  interface establish_self_consistency
     module procedure dense_scc_establish_self_consistency
  endinterface

  interface extrapolate_charges
     module procedure dense_scc_extrapolate_charges
  endinterface

  public :: register
  interface register
     module procedure dense_scc_register
  endinterface

  !
  ! Internal interface
  !

  interface internal_init
     module procedure dense_scc_internal_init
  endinterface internal_init

contains

  !>
  !! Constructor
  !!
  !! Initialize self-consistent charge calculation and allocate memory.
  !<
  subroutine dense_scc_init(this, solv, coul, dq_crit, beta, est, max_nit, andersen_memory, warn, charges_only, log, error)
    use, intrinsic :: iso_c_binding

    implicit none
    
    type(dense_scc_t),              intent(inout) :: this
    type(dense_solver_t), optional, target        :: solv
    type(C_PTR),          optional, intent(in)    :: coul
    real(DP),             optional, intent(in)    :: dq_crit
    real(DP),             optional, intent(in)    :: beta
    real(DP),             optional, intent(in)    :: est
    integer,              optional, intent(in)    :: max_nit
    integer,              optional, intent(in)    :: andersen_memory
    integer,              optional, intent(in)    :: warn
    logical,              optional, intent(in)    :: charges_only
    logical,              optional, intent(in)    :: log
    integer,              optional, intent(inout) :: error

    ! ---

    INIT_ERROR(error)

    this%enabled = .true.

    this%niterations = 0
    this%nsteps      = 0
    this%nfail       = 0

    ! - set defaults and requested parameters
    call set(this, dq_crit, beta, est, max_nit, andersen_memory, warn, log, charges_only)

    ! - set pointers
    if (present(solv)) then
       call set_solver(this, solv, error)
       PASS_ERROR(error)
    endif
    if (present(coul)) then
       call set_Coulomb(this, coul, error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_scc_init


  !>
  !! Constructor
  !!
  !! Set SCC parameters
  !<
  subroutine dense_scc_set(this, dq_crit, beta, est, max_nit, andersen_memory, warn, log, charges_only)
    implicit none
    
    type(dense_scc_t), intent(inout)  :: this
    real(DP), intent(in), optional         :: dq_crit
    real(DP), intent(in), optional         :: beta
    real(DP), intent(in), optional         :: est
    integer, intent(in), optional          :: max_nit
    integer, intent(in), optional          :: andersen_memory
    integer, intent(in), optional          :: warn
    logical, intent(in), optional          :: log
    logical, intent(in), optional          :: charges_only

    ! ---

    ! - set defaults and requested parameters

    if (present(dq_crit)) then
       this%dq_crit  = dq_crit
    endif

    if (present(beta)) then
       this%beta  = beta
    endif

    if (present(max_nit)) then
       this%max_nit  = max_nit
    endif

    if (present(andersen_memory)) then
       this%andersen_memory  = andersen_memory
    endif

    if (present(warn)) then
       this%warn  = warn
    endif

    if (present(log)) then
       this%log  = log
    endif

    if (present(charges_only)) then
       this%charges_only  = charges_only
    endif

  endsubroutine dense_scc_set


  !>
  !! Destructor
  !!
  !! Remove all buffers from memory
  !<
  subroutine dense_scc_del(this)
    implicit none

    type(dense_scc_t), intent(inout)  :: this

    ! ---

    if (this%niterations > 0) then
       call prlog("- dense_scc_del -")
       call prlog("Solver was called "//this%nsteps//" times.")
       call prlog("Average number of iterations for self-consistency = "//((1.0_DP*this%nsteps)/this%niterations))
       call prlog("Convergence failed "//this%nfail//" times.")
       call prlog
    endif

    if (allocated(this%phi))  deallocate(this%phi)
    if (allocated(this%r))    deallocate(this%r)
    if (allocated(this%q))    deallocate(this%q)

    this%p   => NULL()
    this%tb  => NULL()

  endsubroutine dense_scc_del


  !>
  !! Set the associated Coulomb solver object (internal)
  !!
  !! Set the associated Coulomb solver object (internal)
  !<
  subroutine dense_scc_set_Coulomb(this, coul, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(dense_scc_t), intent(inout) :: this
    type(C_PTR),       intent(in)    :: coul  !< Coulomb object to be associated
    integer, optional, intent(out)   :: error !< Error signals

    ! ---

    this%coul = coul

    call internal_init(this, error)
    PASS_ERROR(error)

  endsubroutine dense_scc_set_Coulomb


  !>
  !! Set the associated eigenvalue solver object (internal)
  !!
  !! Set the associated eigenvalue solver object (internal)
  !<
  subroutine dense_scc_set_solver(this, solv, error)
    implicit none

    type(dense_scc_t),              intent(inout)  :: this
    type(dense_solver_t), target                   :: solv
    integer,              optional, intent(out)    :: error   !< Error signals

    ! ---

    this%solv => solv

    call internal_init(this, error)
    PASS_ERROR(error)

  endsubroutine dense_scc_set_solver


  !>
  !! Set the associated particles and Hamiltonian object
  !!
  !! Set the associated particles and Hamiltonian object
  !<
  subroutine dense_scc_bind_to(this, p, tb, error)
    implicit none

    type(dense_scc_t),                   intent(inout)  :: this
    type(particles_t),         target                   :: p        !< Particles to be associated
    type(dense_hamiltonian_t), target                   :: tb       !< Hamiltonian to be associated
    integer,                   optional, intent(out)    :: error   !< Error signals

    ! ---

    INIT_ERROR(error)

    this%p  => p

    ! allocate phi array
    if (allocated(this%phi)) then
       deallocate(this%phi)
    endif
    if (p%maxnatloc <= 0) then
       RAISE_ERROR("scc_init: Particles doesn't seem to contain any atoms.", error)
    endif
    allocate(this%phi(p%maxnatloc))

    ! allocate history
    this%history_counter = 0
    if (allocated(this%r)) then
       deallocate(this%r)
    endif
    if (allocated(this%q)) then
       deallocate(this%q)
    endif
    if (this%extrapolation_memory >= 2) then
       allocate(this%r(3, p%maxnatloc, this%extrapolation_memory))
       allocate(this%q(p%maxnatloc, this%extrapolation_memory))
    endif

    this%tb    => tb

    call internal_init(this, error)
    PASS_ERROR(error)

  end subroutine dense_scc_bind_to


  !>
  !! Set the associated particles and Hamiltonian object
  !!
  !! Set the associated particles and Hamiltonian object
  !<
  subroutine dense_scc_internal_init(this, error)
    implicit none

    type(dense_scc_t),           intent(inout)  :: this
    integer,           optional, intent(out)    :: error   !< Error signals

    ! ---

    INIT_ERROR(error)

    if (c_associated(this%coul) .and. &
        associated(this%p) .and. &
        associated(this%tb)) then

       ! Report
       call prlog("- dense_scc_internal_init -")
       call prlog("dq_crit              = "//this%dq_crit)
       call prlog("miximg               = "//this%beta)
       call prlog("andersen_memory      = "//this%andersen_memory)
       call prlog("max_nit              = "//this%max_nit)
       call prlog("extrapolation_memory = "//this%extrapolation_memory)
       call prlog

       call dense_scc_copy_Hubbard_U(this%coul, this%p, this%tb, error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_scc_internal_init


  !>
  !! Copy Hubbard-U to ChargeOverlap module (internal)
  !!
  !! Copy Hubbard-U to ChargeOverlap module (internal)
  !<
  subroutine dense_scc_copy_Hubbard_U(coul, p, tb, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(C_PTR),               intent(inout) :: coul
    type(particles_t),         intent(in)    :: p
    type(dense_hamiltonian_t), target        :: tb
    integer,         optional, intent(inout) :: error

    ! ---

    integer               :: i
    real(DP)              :: U(p%nel)
    type(notb_element_t)  :: el
    integer               :: z         ! temp

    type(materials_t), pointer  :: tb_mat

    ! ---

    INIT_ERROR(error)
    call c_f_pointer(tb%mat, tb_mat)

    ! - checks

    if(.not. associated(tb_mat)) then
       RAISE_ERROR("scc_copy_Hubbard_U: The materials database does not seem to be associated to the Hamiltonian yet.", error)
    end if

    ! - copy U's

    do i = 1, p%nel
       if (element_by_Z(tb_mat, p%el2Z(i), el)) then
          U(i)  = el%U
       else
          !RAISE_ERROR("scc_copy_Hubbard_U: Element with atomic number '" // p%el2Z(i) // "' encountered, but not found in the materials database.", error)
          z = p%el2Z(i)
          WARN("scc_copy_Hubbard_U: Element with atomic number '" // z // "' not found in materials database, setting Hubbard U to zero.")
          U(i) = 0.0_DP
       endif
    enddo

    call coulomb_set_Hubbard_U(coul, p, U, error)
    PASS_ERROR(error)

  endsubroutine dense_scc_copy_Hubbard_U


  !>
  !! Run the self-consistency loop
  !!
  !! Run the self-consistency loop
  !<
  subroutine dense_scc_establish_self_consistency(this, p, nl, tb, q, noc, f, error)
    implicit none

    type(dense_scc_t), intent(inout)   :: this      !< SCC object
    type(particles_t), target          :: p         !< Particles
    type(neighbors_t), target          :: nl        !< Neighbor list
    type(dense_hamiltonian_t), target  :: tb
    real(DP), intent(inout)            :: q(p%nat)  !< Charges
    real(DP), intent(in)               :: noc       !< Number of occupied orbitals
    integer, intent(in), optional      :: f         !< Filter for atom types
    integer, intent(inout), optional   :: error     !< Error signals

    ! ---

    integer   :: it                                 ! SCC loop iteration
    logical   :: done                               ! charges converged?
    integer   :: nf                                 ! number of atoms within filter
    integer   :: filter                             ! filter used

    type(anderson_mixer_t) :: mixer                 ! Anderson mixer

    real(DP)  :: f_q_prev(p%natloc), f_q_new(p%natloc)    ! filtered charge arrays, previous and new
    real(DP)  :: prev_mu

#ifdef DEBUG_ENERGY_MINIMUM
    integer :: M
#endif

    ! ---

    INIT_ERROR(error)

    ASSERT_ASSOCIATION(this%p, p, error)
    ASSERT_ASSOCIATION(this%tb, tb, error)

    if (.not. associated(this%solv)) then
       RAISE_ERROR("dense_scc_establish_self_consistency: No eigenvalue solver specified.", error)
    endif

    if (.not. c_associated(this%coul)) then
       RAISE_ERROR("dense_scc_establish_self_consistency: No Coulomb solver specified.", error)
    endif

    !
    ! XXX: Test. Is this correct?
    !
    if(this%charges_only) then
       this%phi = 0.0_DP
       call diag_start(this%solv, tb, error=error)
       PASS_ERROR(error)
       call diag_HS(this%solv, tb, noc, error=error)
       PASS_ERROR(error)
       call diag_stop(this%solv, tb, error=error)
       PASS_ERROR(error)
       call mulliken(this%solv, tb, q, error=error)
       PASS_ERROR(error)
#ifndef LAMMPS
       call I_changed_other(p)
#endif
       return
    end if

    this%niterations = this%niterations + 1

    !
    ! Extrapolate charges
    !

    if (this%extrapolation_memory >= 2) then
       call extrapolate_charges(this, p, q, error=error)
       PASS_ERROR(error)
    endif

    !
    ! Init
    !

    ! - Init mixed and solver
    call init(mixer, this%andersen_memory)
    call diag_start(this%solv, tb)

    !
    ! Pack charges within filter, which are the ones we want to be changing
    !

    ! init filter
    if(present(f)) then
       filter = f
    else
       filter = filter_from_string("*", p, ierror=error)
       PASS_ERROR(error)
    end if

    call timer_start("scc_establish_self_consistency")

    ! pack charges from q to f_q_prev
    call filter_pack(filter, p, q, f_q_prev)

    ! number of atoms
    nf = filter_count(filter, p)

    !
    ! Logging
    !

    if (this%log) then
       write (ilog, '(1X,A10,1X,A4,4A12)')  "scc|", "it", "sum(q)", "max(dq)", "mu[eV]", "dmu[eV]"
    endif

    !
    ! Charge self-consistency loop
    !

    done = .false.
    it = 0
    prev_mu = this%tb%mu
    do while(.not. done .and. it < this%max_nit)
       it          = it + 1
       this%nsteps = this%nsteps + 1

       ! solve: calculate potential -> diagonalize -> calculate new charges
       call solve(error=error)
       PASS_ERROR_AND_STOP_TIMER("scc_establish_self_consistency", error)

       ! new charges from q to f_q_new
       call filter_pack(filter, p, q, f_q_new)

       ! mix, new charges to f_q_prev
       call mix(mixer, it, nf, f_q_prev, f_q_new, this%beta, this%dq_crit, done, 0.05d0)

       ! unpack new charges from f_q_prev to q
       call filter_unpack(filter, p, f_q_prev, q)

#ifndef LAMMPS
       ! notify p that charges changed
       call I_changed_other(p)
#endif

       ! output and logging
       if( mod(it, this%warn)==0 .or. (done .and. it>this%warn) ) call prscrlog("Warning: Charge self-consistency at iteration "//it//".")

       if (this%log) then
          if (it > 1) then
             write (ilog, '(12X,I4,F12.3,3ES12.3)')  it, sum(q), maxval( abs(f_q_prev(1:nf) - f_q_new(1:nf)) ), this%tb%mu, this%tb%mu-prev_mu
          else
             write (ilog, '(12X,I4,F12.3,2ES12.3)')  it, sum(q), maxval( abs(f_q_prev(1:nf) - f_q_new(1:nf)) ), this%tb%mu
          endif
       endif

       prev_mu = this%tb%mu

    enddo  ! end of charge self-consistency loop

    !
    ! Do the rest
    !

    ! delete mixer
    call del(mixer)

    ! warn of problems with convergence
    if (it >= this%max_nit) then
       call prscrlog("Warning: Maximum number of SCC iterations (= "//this%max_nit//") exceeded.")
       this%nfail = this%nfail + 1
    endif

#ifdef DEBUG_ENERGY_MINIMUM

    !
    ! Probe if this is really the minimum energy configuration by
    ! adding random perturbations to the charge.
    !

    ! XXX: This may not work!

    call coulomb_charge_changed(part, nl)
    phi  = 0.0_DP
    call coulomb_potential(part, nl, phi_in)
    call diag_HS(this%solv, tb, part, noc, phi_in, error=error)
    PASS_ERROR_AND_STOP_TIMER("scc_establish_self_consistency", error)
    e1 = e_bs(solver)
    e2 = 0.0_DP
    call coulomb_force(part, nl, e2)

    max = e1+e2

    do it = 1, 100

       do M = 1, p%natloc-1
          q(M)  = q(M) + rng_uniform(-0.001_DP, 0.001_DP)
       enddo
       q(p%natloc)  = -sum(q(1:p%natloc-1))

       call coulomb_charge_changed(part, nl)
       phi_in  = 0.0_DP
       call coulomb_potential(part, nl, phi_in)
       call diag_HS(solver, tb, part, noc, phi_in, error=error)
       PASS_ERROR_AND_STOP_TIMER("scc_establish_self_consistency", error)
       e1 = e_bs(solver)
       e2 = 0.0_DP
       call coulomb_force(part, nl, e2)

       write (*, '(4ES20.10)')  e1, e2, e1+e2, e1+e2-max

       q  = dq1

    enddo

#endif

    call timer_stop("scc_establish_self_consistency")

  contains

    subroutine solve(error)
      implicit none

      integer, intent(inout), optional  :: error

      ! ---

      this%phi  = 0.0_DP
      call coulomb_potential(this%coul, p, nl, q, this%phi, ierror=error)
      PASS_ERROR(error)

      call diag_HS(this%solv, tb, noc, this%phi, error=error)
      PASS_ERROR(error)

      call mulliken(this%solv, tb, q)

    endsubroutine solve

  endsubroutine dense_scc_establish_self_consistency


  !>
  !! Extrapolate charges from past time steps
  !<
  subroutine dense_scc_extrapolate_charges(this, p, q, error)
    implicit none

    type(dense_scc_t), intent(inout) :: this      !< SCC object
    type(particles_t), intent(in)    :: p         !< Particles
    real(DP),          intent(inout) :: q(p%nat)  !< Charges
    integer, optional, intent(out)   :: error     !< Error status

    ! ---

    integer :: i, k, l

    real(DP) :: a(this%extrapolation_memory-1, this%extrapolation_memory-1)
    real(DP) :: b(this%extrapolation_memory-1)
    real(DP) :: alpha(this%extrapolation_memory-1)
    real(DP) :: q0(p%nat)

    real(DP) :: drk(3, p%natloc), drl(3, p%natloc)

    ! ---

    INIT_ERROR(error)

    q0 = q

    if (this%history_counter >= this%extrapolation_memory) then
       do k = 1, this%extrapolation_memory-1
          drk = this%r(1:3, 1:p%natloc, modulo(this%history_counter-k, this%extrapolation_memory)+1) - &
                this%r(1:3, 1:p%natloc, modulo(this%history_counter-k-1, this%extrapolation_memory)+1)
          do l = 1, this%extrapolation_memory-1
             drl = this%r(1:3, 1:p%natloc, modulo(this%history_counter-l, this%extrapolation_memory)+1) - &
                   this%r(1:3, 1:p%natloc, modulo(this%history_counter-l-1, this%extrapolation_memory)+1)
             a(k, l) = dot_product(reshape(drk, [3*p%natloc]), reshape(drl, [3*p%natloc]))
          enddo
          b(k) = dot_product( &
              reshape(PCN3(p, 1:p%natloc) - &
                      this%r(1:3, 1:p%natloc, modulo(this%history_counter-1, this%extrapolation_memory)+1), &
                      [3*p%natloc]), &
              reshape(this%r(1:3, 1:p%natloc, modulo(this%history_counter-k, this%extrapolation_memory)+1) - &
                      this%r(1:3, 1:p%natloc, modulo(this%history_counter-k-1, this%extrapolation_memory)+1), &
                      [3*p%natloc]) &
              )
       enddo

       alpha = matmul(inverse(a, error=error), b)
       if (error == ERROR_NONE) then
          !                q(t) - q(t-dt)
          q = q0 + alpha(1)*(q0 - this%q(1:p%nat, modulo(this%history_counter-1, this%extrapolation_memory)+1))
          do i = 2, this%extrapolation_memory-1
             q = q + alpha(i)*(this%q(1:p%nat, modulo(this%history_counter-i+1, this%extrapolation_memory)+1) - &
                               this%q(1:p%nat, modulo(this%history_counter-i, this%extrapolation_memory)+1))
          enddo
       else
          call prlog("Warning: Unable to extrapolate charges. Resetting history.")
          call clear_error(error)
          this%history_counter = 0
       endif
    endif

    this%history_counter = this%history_counter+1
    i = modulo(this%history_counter-1, this%extrapolation_memory)+1
    ! This is current r(t+dt)
    this%r(1:3, 1:p%nat, i) = PCN3(p, 1:p%nat)
    ! This is last q(t)
    this%q(1:p%nat, i)      = q0(1:p%nat)
    
  endsubroutine dense_scc_extrapolate_charges


  !>
  !! Register the SCC object
  !<
  subroutine dense_scc_register(this, cfg)
    implicit none

    type(dense_scc_t), target, intent(inout) :: this
    type(c_ptr),               intent(in) :: cfg

    ! ---

    type(c_ptr) :: m

    ! ---

    this%max_nit              = 200         !< max nr of iterations in self-consistency
    this%dq_crit              = 0.0001_DP   !< limit for convergence = max(dq_in-qd_out)
    this%beta                 = 0.2_DP      !< density mixing parameter in iteration
    this%andersen_memory      = 3           !< M in Anderson mixing in iteration

    this%warn                 = 20          !< warn after 20 iterations
    this%log                  = .false.     !< write a status report for each SCC step

    this%charges_only         = .false.     !< Only calculate Mulliken charges

    this%extrapolation_memory = 3

    m = ptrdict_register_module(cfg, c_loc(this%enabled), CSTR("SCC"), &
         CSTR("Use charge self-consistency in the tight-binding calculation."))

    call ptrdict_register_real_property(m, c_loc(this%dq_crit), &
         CSTR("dq_crit"), &
         CSTR("Convergence criterium for the self-consistent determination of the charges."))
    call ptrdict_register_real_property(m, c_loc(this%beta), CSTR("mixing"), &
         CSTR("Mixing parameter for charge self-consistency."))
    call ptrdict_register_integer_property(m, c_loc(this%max_nit), &
         CSTR("maximum_iterations"), &
         CSTR("Maximum number of SCC iterations."))

    call ptrdict_register_integer_property(m, c_loc(this%andersen_memory), &
         CSTR("andersen_memory"), &
         CSTR("Andersen mixing memory."))

    call ptrdict_register_integer_property(m, c_loc(this%extrapolation_memory), &
         CSTR("extrapolation_memory"), &
         CSTR("Number of past time steps to consider for charge extrapolation (min 2, extrapolation disabled if less)."))

    call ptrdict_register_integer_property(m, c_loc(this%warn), CSTR("warn"), &
         CSTR("Warn after a number of iterations without self-consistency."))

    call ptrdict_register_boolean_property(m, c_loc(this%log), CSTR("log"), &
         CSTR("Print a status for each iteration step to the log file."))

  endsubroutine dense_scc_register

endmodule dense_scc

