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
! @meta
!   shared
!   classtype:variable_charge_t classname:VariableCharge interface:potentials
! @endmeta

!> Charge transfer
!!
!! Module for the evaluation of charge transfer. Does only work
!! in combination with a Coulomb (and possibly charge overlap) module. Uses
!! either a conjugate gradient algorithm, anderson mixing or
!! Car-Parrinello to optimize the charges.
!!
!! The method is described in:
!!
!!   A. K. Rappe and W. A. Goddard III, J. Phys. Chem. 95, 3358 (1991)
!!   F. H. Streitz and J. W. Mintmire, Phys. Rev. B 50, 11996 (1994)
!!   Pastewka, Jaervi, Mayrhofer, Moseler, Phys. Rev. B 83, 165418 (2011)
!!
!! The charges are optimized by minimizing the energy
!!
!! \f[
!!   E = \frac{1}{p} V \sum_i |q_i|^p + \frac{1}{2} \sum_{i,j \neq i} \gamma_{ij} q_i q_j + \frac{1}{2} U \sum_i q_i^2 + \sum_i q_i \left( \phi_i + \chi_i \right).
!! \f]
!!
!! The three first terms are: band structure energy, charge overlap + Coulomb tail,
!! self-energy. In the last term, \f$\phi_i\f$ is the potential at
!! charge \f$i\f$ from other charges and possible external potential
!! and \f$\chi\f$ can be used to adjust the chemical potential of
!! the electron reservoir of atom \f$i\f$.
!!
!! A Coulomb module is used to calculate the contribution from
!! \f$\phi_i\f$, in which the charge overlap can also be
!! included by the appropriate module.
!!
!! The potential energy contribution from the charge transfer module
!! only includes the band structure energy and the contribution from
!! \f$\chi\f$. Although the others are calculated during minimization
!! they are assumed to be added by the Coulomb and charge overlap modules
!! to the total energy.
!!
!! The parameters are read in from a file 'variable_charge.dat' of the
!! following format:
!!
!! 1                         [number of atom types]
!!
!!---
!!
!!C                          [element for first atom type]
!!
!!1                          [group of first atom type (6th column of atoms.dat)]
!!
!!X = 0.0                    [\f$ \chi \f$ of GROUP]
!!
!!U = 9.93   ! eV/e**2       [U of ELEMENT]
!!
!!V = 17.0   ! eV            [V of ELEMENT]
!!
!!p = 2.0                    [p of ELEMENT]
!!
!!---
!!
!! ADD BOTH INSTRUCTIONS HOW TO CALL FROM OWN MAIN PROGRAM AND INPUT ENTRIES
!! FOR STANDALONE CODE!
!!
!<

#include "macros.inc"
#include "filter.inc"

module variable_charge
  use supplib

  use particles
  use neighbors

  use anderson_mixer
  use extrapolation

  use filter

  use coulomb

  implicit none

  private

  integer, parameter  :: ST_CONJUGATE_GRADIENTS  = 0   !< Conjugate-gradients minimization
  integer, parameter  :: ST_ANDERSON             = 1   !< Anderson mixing
  integer, parameter  :: ST_CAR_PARRINELLO       = 2   !< Car-Parrinello type fictious dynamics
  integer, parameter  :: ST_DISABLE              = 3

  character(*), parameter  :: QV_STR  = "q_velocities"
  character(*), parameter  :: QA_STR  = "q_accelerations"

  !>
  !! Element / group data for reading input file
  !<
  type ct_element_t

     character(2)  :: name   !< Atom type
     integer       :: group  !< Group id

     integer       :: Z      !< Atomic number of type

     real(DP)      :: X      !< Chemical potential (see energy expression above) (per group)
     real(DP)      :: U      !< Hubbard U (per element)
     real(DP)      :: V      !< Band structure prefactor (E includes term (V/2)*|q|^p (per element)
     real(DP)      :: Vp     !< Band structure exponent (per element)

  endtype ct_element_t
  public :: ct_element_t

  
  !>
  !! Variable charge object
  !!
  !! Variable charge object
  !<
  type variable_charge_t

     !
     ! Stuff read from the input file
     !

     character(MAX_EL_STR)  :: elements      = "*"              !< Which elements to include?

     integer                :: solver_type   = ST_ANDERSON      !< Solver type

     real(DP)               :: total_charge  = 0.0_DP           !< Total charge (of subsystem included in equilibration)
     
     logical(BOOL)          :: log           = .false.          !< log? (only very little logging)
     logical(BOOL)          :: trace         = .false.          !< detailed logging? (turns log true as well)


     !
     ! For the Anderson solver
     !

     real(DP)               :: convergence      = 0.001_DP      !< Anderson: converge the charge equilibration up to this point
     real(DP)               :: freq             = -1.0_DP       !< Anderson: charge update frequency
     integer                :: anderson_memory  = 3             !< Anderson: history of mixer
     real(DP)               :: beta             = 0.5           !< Anderson: mixing parameter
     real(DP)               :: beta_max         = 0.5           !< Anderson: maximum value of the mixing parameter
     real(DP)               :: beta_mul         = 1.0           !< Anderson: increase mixing parameter by this factor up to beta_max
     integer                :: max_iterations   = 100           !< Anderson: maximum number of iterations

     !
     ! For the Car-Parrinello solver
     !

     real(DP)               :: mq                = 10.0_DP      !< C-P: fictitious mass
     real(DP)               :: gamma             = 1.0_DP       !< C-P: damping constant
     real(DP)               :: dE_max            = 0.01_DP      !< C-P: maximum |dE| at which an Anderson mixing step is called
     integer                :: fail_max          = 100          !< C-P: number of failures before an Anderson mixing step is called

     !
     ! First iteration?
     !

     logical                :: first_iteration
     integer                :: n_fail

     !
     ! Background charge
     !

     real(DP)               :: q0

     !
     ! Counter
     !

     real(DP)               :: ti

     !
     ! Filter
     !

     integer                :: f                    ! Filter on elements for charge transfer

     !
     ! Output file
     !

     integer                :: dE_un  = -1

     !
     ! Chemical potential per group
     !

     real(DP), allocatable  :: X(:)

     !
     ! Band structure energy and Hubbard U per element
     !

     real(DP), allocatable  :: U(:)
     real(DP), allocatable  :: V(:)
     real(DP), allocatable  :: Vp(:)

     !
     ! Mixer
     !

     type(anderson_mixer_t)  :: mixer

     !
     ! Internal buffers
     !

     real(DP), allocatable   :: phi(:)
     real(DP), allocatable   :: r_in(:)
     real(DP), allocatable   :: r_out(:)

     !
     ! Charge velocities
     !

     real(DP), pointer  :: qv(:)    !< C-P: Velocities for the charges
     real(DP), pointer  :: qa(:)    !< C-P: Accelerations for the charges

     !
     ! Chemical potential, Hubbard-U's etc 
     !
     
     type(ct_element_t), allocatable  :: at(:)

     !
     ! Position and charge history
     !

     integer               :: extrapolation_memory = 0  !< Number of past steps to keep
     type(extrapolation_t) :: extrapolation

  endtype variable_charge_t
  public :: variable_charge_t

  public :: init
  interface init
     module procedure variable_charge_init
  endinterface

  public :: set
  interface set
     module procedure variable_charge_init
  endinterface

  public :: del
  interface del
     module procedure variable_charge_del
  endinterface

  public :: bind_to_with_coul
  interface bind_to_with_coul
     module procedure variable_charge_bind_to_with_coul
  endinterface

  public :: energy_and_forces_with_charges_and_coul
  interface energy_and_forces_with_charges_and_coul
     module procedure variable_charge_energy_and_forces_with_charges_and_coul
  endinterface

  public :: get_epot
  interface get_epot
     module procedure variable_charge_get_epot
  endinterface

  public :: register
  interface register
     module procedure variable_charge_register
  endinterface register

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine variable_charge_init(this, &
       elements, at, solver_type, total_charge, convergence, &
       anderson_memory, extrapolation_memory, beta, beta_max, beta_mul, &
       mq, gamma, dE_max, fail_max, &
       max_iterations, &
       log, trace)
    implicit none

    type(variable_charge_t),      intent(inout) :: this

    character(*),       optional, intent(in)    :: elements
    type(ct_element_t), optional, intent(in)    :: at(:)
    integer,            optional, intent(in)    :: solver_type
    real(DP),           optional, intent(in)    :: total_charge
    real(DP),           optional, intent(in)    :: convergence
    integer,            optional, intent(in)    :: anderson_memory
    integer,            optional, intent(in)    :: extrapolation_memory
    real(DP),           optional, intent(in)    :: beta
    real(DP),           optional, intent(in)    :: beta_max
    real(DP),           optional, intent(in)    :: beta_mul
    real(DP),           optional, intent(in)    :: mq
    real(DP),           optional, intent(in)    :: gamma
    real(DP),           optional, intent(in)    :: dE_max
    real(DP),           optional, intent(in)    :: fail_max
    integer,            optional, intent(in)    :: max_iterations
    logical,            optional, intent(in)    :: log
    logical,            optional, intent(in)    :: trace

    ! ---

    ASSIGN_PROPERTY(elements)
    ASSIGN_PROPERTY(solver_type)
    ASSIGN_PROPERTY(total_charge)
    ASSIGN_PROPERTY(convergence)
    ASSIGN_PROPERTY(anderson_memory)
    ASSIGN_PROPERTY(extrapolation_memory)
    ASSIGN_PROPERTY(beta)
    ASSIGN_PROPERTY(beta_max)
    ASSIGN_PROPERTY(beta_mul)
    ASSIGN_PROPERTY(mq)
    ASSIGN_PROPERTY(gamma)
    ASSIGN_PROPERTY(dE_max)
    ASSIGN_PROPERTY(fail_max)
    ASSIGN_PROPERTY(max_iterations)
    ASSIGN_PROPERTY(log)
    ASSIGN_PROPERTY(trace)

    if (present(at)) then
       if (allocated(this%at)) then
          deallocate(this%at)
       endif
       allocate(this%at(lbound(at, 1):ubound(at, 1)))
       this%at = at
    endif

  endsubroutine variable_charge_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine variable_charge_del(this)
    implicit none

    type(variable_charge_t), intent(inout)  :: this

    ! ---

    if (allocated(this%at)) then
       deallocate(this%at)
    endif

    if (allocated(this%X)) then
       deallocate(this%X)
    endif
    if (allocated(this%U)) then
       deallocate(this%U)
    endif
    if (allocated(this%V)) then
       deallocate(this%V)
    endif
    if (allocated(this%Vp)) then
       deallocate(this%Vp)
    endif

    if (allocated(this%r_in)) then
       deallocate(this%r_in)
    endif
    if (allocated(this%r_out)) then
       deallocate(this%r_out)
    endif
    if (allocated(this%phi)) then
       deallocate(this%phi)
    endif

    call del(this%mixer)
    call del(this%extrapolation)

    if (this%dE_un > 0) then
       call fclose(this%dE_un)
       this%dE_un  = -1
    endif

  endsubroutine variable_charge_del


  !>
  !! Initialize the variable charge module once all information is available
  !!
  !! Initialize the variable charge module once all information is available
  !<
  subroutine variable_charge_bind_to_with_coul(this, p, nl, coul, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(variable_charge_t), intent(inout) :: this
    type(particles_t),       intent(in)    :: p
    type(neighbors_t),       intent(inout) :: nl
    type(c_ptr),             intent(in)    :: coul
    integer,       optional, intent(inout) :: ierror

    ! ---

    integer                          :: un, nel, i, k, min_group, max_group, Z
    character(200)                   :: line, act, dat
    logical                          :: done
    logical                          :: gotU, gotV, gotVp, gotX

    ! ---

    call prlog("- variable_charge_bind_to_with_coul -")

    this%f = filter_from_string(this%elements, p, ierror=ierror)
    PASS_ERROR(ierror)

    !
    ! Read parameters from variable_charge.dat
    !

    if (.not. allocated(this%at)) then
       call prlog("     Reading 'variable_charge.dat'")
       un = fopen("variable_charge.dat")
       read (un, *)  nel
       read (un, *)

       allocate(this%at(nel))

       ! loop over elements
       do i = 1, nel
          read (un, *)  this%at(i)%name
          read (un, *)  this%at(i)%group

          Z = atomic_number(this%at(i)%name)
          if (Z <= 0 .or. Z > MAX_Z) then
             RAISE_ERROR("Unknown element '" // trim(this%at(i)%name) // "'.", ierror)
          endif

          this%at(i)%Z  = Z

          this%at(i)%X  = 0.0_DP
          this%at(i)%U  = 0.0_DP
          this%at(i)%V  = 0.0_DP
          this%at(i)%Vp  = 0.0_DP

          ! loop over parameters
          gotU = .false.
          gotX = .false.
          gotV = .false.
          gotVp = .false.
          done = .false.
          do while (.not. done)
             read (un, '(A)')  line

             k = scan(line,'=')
             if (k /= 0) then
                act = adjustl(line(1:k-1))
                dat = line(k+1:)
             else
                act = line
             endif

             select case(trim(act))
             case("X")
                read (dat, *)  this%at(i)%X
                gotX = .true.
             case("U")
                read (dat, *)  this%at(i)%U
                gotU = .true.
             case("V")
                read (dat, *)  this%at(i)%V
                gotV = .true.
             case("p")
                read (dat, *)  this%at(i)%Vp
                gotVp = .true.
             case default
                done = .true.
             endselect
          enddo  ! end of loop over parameters

          if(.not. (gotU .and. gotX .and. gotV .and. gotVp)) then
             RAISE_ERROR("Did not find all of X, U, V, p for atom " // this%at(i)%name // ", group " // this %at(i)%group, ierror)
          end if

       enddo  ! end of loop over elements

       call fclose(un)
    endif

    !
    ! Fill data according to group and element
    !

    min_group = 0
    max_group = 0

    do i = 1, p%natloc
       min_group = min(min_group, p%g(i))
       max_group = max(max_group, p%g(i))
    enddo
#ifdef _MP
    min_group = min(mod_parallel_3d%mpi, min_group)
    max_group = max(mod_parallel_3d%mpi, max_group)
#endif

    ! reallocate arrays
    if (allocated(this%X)) deallocate(this%X)
    if (allocated(this%U)) deallocate(this%U)
    if (allocated(this%V)) deallocate(this%V)
    if (allocated(this%Vp)) deallocate(this%Vp)
    allocate(this%X(min_group:max_group))
    allocate(this%U(p%nel))
    allocate(this%V(p%nel))
    allocate(this%Vp(p%nel))
    call log_memory_estimate(this%X, "X")
    call log_memory_estimate(this%U, "U")
    call log_memory_estimate(this%V, "V")
    call log_memory_estimate(this%Vp, "p")

    ! Set variables
    this%X  = 0.0_DP
    this%U  = 0.0_DP
    this%V  = 0.0_DP
    this%Vp  = 0.0_DP
    do i = lbound(this%at, 1), ubound(this%at, 1)
       if (.not. IS_EL2(this%f, p%Z2el(this%at(i)%Z))) then
          RAISE_ERROR("Element '" // trim(this%at(i)%name) // "' defined in the charge-transfer file (variable_charge.dat), but not in the simulation control file (md.dat).", ierror)
       endif

       if (this%at(i)%Z <= 0 .or. this%at(i)%Z > MAX_Z) then
          call prlog("     Could not find atom with atomic number " // this%at(i)%Z // ".")
       else
          call prlog("     " // ElementName(this%at(i)%Z) // "(" // this%at(i)%group // ")")
          call prlog("     - X = " // this%at(i)%X // ", U = " // this%at(i)%U // ", V = " // this%at(i)%V // ", p = " // this%at(i)%Vp)
       endif

       if (this%at(i)%group < min_group .or. this%at(i)%group > max_group) then
          RAISE_ERROR("Atoms with group number " // this%at(i)%group // " do not exists in this simulation.", ierror)
       endif

       this%X(this%at(i)%group)      = this%at(i)%X
       this%U(p%Z2el(this%at(i)%Z))  = this%at(i)%U
       this%V(p%Z2el(this%at(i)%Z))  = this%at(i)%V
       this%Vp(p%Z2el(this%at(i)%Z))  = this%at(i)%Vp
    enddo

    call coulomb_set_Hubbard_U(coul, p, this%U)

    if (this%trace) this%log = .true.

    this%ti = this%freq

    k = filter_count(this%f, p)
#ifdef _MP
    call sum_in_place(mod_parallel_3d%mpi, k)
#endif
    this%q0 = this%total_charge / k

    write (ilog, '(5X,A,F20.10)')  "total_charge         = ", this%total_charge
    write (ilog, '(5X,A,F20.10)')  "* q0                 = ", this%q0

    call init(this%mixer, this%anderson_memory)

    ! allocate extrapolation history
    call prlog("extrapolation_memory = "//this%extrapolation_memory)
    call init(this%extrapolation, p, this%extrapolation_memory)

    if (this%solver_type == ST_CAR_PARRINELLO) then
       call ptr_by_name(p%data, QV_STR, this%qv, ierror=ierror)
       PASS_ERROR(ierror)
       call ptr_by_name(p%data, QA_STR, this%qa, ierror=ierror)
       PASS_ERROR(ierror)

       if (this%trace) then
          this%dE_un  = fopen("variable_charge_max_dE.out", F_WRITE)
       endif
    endif

    this%first_iteration  = .true.
    this%n_fail           = 0

    write (ilog, *)

  endsubroutine variable_charge_bind_to_with_coul


  !>
  !! Energy from band structure and chemical potential (internal)
  !!
  !! Energy of a single atom (charge) from the band structure (V)
  !!  term in the energy expression.
  !!
  !! So, here
  !!
  !! \f[
  !!   E = \frac{V}{p} |q|^p
  !! \f]
  !!
  !<
  pure function E_V(q, V, p)
    implicit none

    real(DP), intent(in)                   :: q    !< charge
    real(DP), intent(in)                   :: V    !< prefactor of charge in band structure energy
    real(DP), intent(in)                   :: p    !< exponent of charge in band structure energy
    real(DP)                               :: E_V

    if (q == 0.0_DP) then
       E_V = 0.0_DP
    else
       E_V = (V/p)*abs(q)**p
    endif

  end function E_V

  ! XXX
  pure function dE_V(q, V, p)
    implicit none

    real(DP), intent(in)                   :: q    !< charge
    real(DP), intent(in)                   :: V    !< prefactor of charge in band structure energy
    real(DP), intent(in)                   :: p    !< exponent of charge in band structure energy
    real(DP)                               :: dE_V

!    write (*, *)  q, V, p

    if (q == 0.0_DP) then
       dE_V = 0.0_DP
    else
       dE_V = sign(1.0_DP, q)*V*abs(q)**(p-1.0_DP)
    endif

  end function dE_V

  pure function d2E_V(q, V, p)
    implicit none

    real(DP), intent(in)                   :: q    !< charge
    real(DP), intent(in)                   :: V    !< prefactor of charge in band structure energy
    real(DP), intent(in)                   :: p    !< exponent of charge in band structure energy
    real(DP)                               :: d2E_V

!    write (*, *)  q, V, p

    if (q == 0.0_DP) then
       d2E_V = 0.0_DP
    else
       d2E_V = V*(p-1.0_DP)*abs(q)**(p-2.0_DP)
    endif

  end function d2E_V


  !>
  !! Evaluate the new charges
  !!
  !! Evaluate the new charges
  !<
  subroutine variable_charge_energy_and_forces_with_charges_and_coul(this, p, &
       nl, q, coul, epot, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(variable_charge_t), intent(inout) :: this
    type(particles_t),       intent(inout) :: p
    type(neighbors_t),       intent(inout) :: nl
    real(DP),                intent(inout) :: q(p%maxnatloc)
    type(C_PTR),             intent(in)    :: coul
    real(DP),                intent(inout) :: epot
    integer,       optional, intent(inout) :: ierror

    ! ---

    integer   :: i

    ! ---

    this%ti  = this%ti + 1.0_DP

    if (this%ti > this%freq .and. this%convergence > 0) then

       if (this%trace) then
          call prlog("- variable_charge_energy_and_forces -")
       endif

       if (this%first_iteration) then
          ! First iteration is always anderson (much faster)
          select case (this%solver_type)

          case (ST_CONJUGATE_GRADIENTS)
             call variable_charge_conjugate_gradients(this, p, nl, q, coul, epot, ierror)
             PASS_ERROR(ierror)
          case (ST_ANDERSON)
             call variable_charge_anderson(this, p, nl, q, coul, epot, ierror)
             PASS_ERROR(ierror)
          case (ST_CAR_PARRINELLO)
             call variable_charge_anderson(this, p, nl, q, coul, epot, ierror)
             PASS_ERROR(ierror)
          case default
             RAISE_ERROR("Internal error: Unknown solver type encountered.", ierror)

          endselect

          this%first_iteration  = .false.
       else

          select case (this%solver_type)

          case (ST_CONJUGATE_GRADIENTS)
             call variable_charge_conjugate_gradients(this, p, nl, q, coul, epot, ierror)
             PASS_ERROR(ierror)
          case (ST_ANDERSON)
             call variable_charge_anderson(this, p, nl, q, coul, epot, ierror)
             PASS_ERROR(ierror)
          case (ST_CAR_PARRINELLO)
             call variable_charge_car_parrinello(this, p, nl, q, coul, epot, ierror)
             PASS_ERROR(ierror)
          case default
             RAISE_ERROR("Internal error: Unknown solver type encountered.", ierror)

          endselect

       endif

    else

       do i = 1, p%natloc
          if (IS_EL(this%f, p, i)) then
             epot = epot + E_V(q(i), this%V(p%el(i)), this%Vp(p%el(i)))
          endif
       enddo

    endif

  endsubroutine variable_charge_energy_and_forces_with_charges_and_coul


  !>
  !! Transformation to the auxilliary space
  !!
  !! Transformation to the auxilliary space
  !<
  subroutine q2aux(p, f, q0, q, naux, r, error)
    implicit none

    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: f
    real(DP), intent(in)           :: q0
    real(DP), intent(in)           :: q(p%natloc)
    integer, intent(out)           :: naux
    real(DP), intent(out)          :: r(p%natloc)
    integer, intent(out), optional :: error

    ! ---

    integer   :: i, j, k

#ifdef _MP
    real(DP)  :: last_q(mod_parallel_3d%mpi%n_procs)
#endif

    ! ---

    INIT_ERROR(error)

    k = 1
    do while (k <= p%natloc .and. .not. IS_EL(f, p, k))
       k = k+1
    enddo

    if (k > p%natloc) then
       RAISE_ERROR("No element matching the filter found.", error)
    endif

    j       = 2
    r(1)    = -q(k) + q0

    do i = k+1, p%natloc

       if (IS_EL(f, p, i)) then
          r(j)  = r(j-1) - q(i) + q0
          j     = j+1
       endif

    enddo

#ifdef _MP

    call mpi_allgather( &
         r(j-1), 1, MPI_DOUBLE_PRECISION, &
         last_q(:), 1, MPI_DOUBLE_PRECISION, &
         mod_parallel_3d%mpi%communicator, i)
    PASS_MPI_ERROR(i, error)

    if (mod_parallel_3d%mpi%my_proc == mod_parallel_3d%mpi%n_procs-1) then
       naux = j-2
    else
       naux = j-1
    endif

    if (mod_parallel_3d%mpi%my_proc /= 0) then
       r(1:naux) = r(1:naux) + sum(last_q(1:mod_parallel_3d%mpi%my_proc))
    endif

#else
    naux = j-2
#endif

  endsubroutine q2aux


  !>
  !! Transformation of gradients to the auxilliary space
  !!
  !! Transformation of gradients to the auxilliary space
  !<
  subroutine dq2aux(p, f, dq, naux, dr, error)
    implicit none

    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: f
    real(DP), intent(in)           :: dq(p%natloc)
    integer, intent(in)            :: naux
    real(DP), intent(out)          :: dr(p%natloc)
    integer, intent(out), optional :: error

    ! ---

    integer  :: i, last_i, j, k

#ifdef _MP
    real(DP)  :: last_q
    integer   :: left, right
#endif

    ! ---

    INIT_ERROR(error)

    k = 1
    do while (k <= p%natloc .and. .not. IS_EL(f, p, k))
       k = k+1
    enddo

    if (k > p%natloc) then
       RAISE_ERROR("No element matching the filter found.", error)
    endif

    last_i  = k
    j       = 1

#ifdef _MP

    left   = modulo(mod_parallel_3d%mpi%my_proc-1, mod_parallel_3d%mpi%n_procs)
    right  = modulo(mod_parallel_3d%mpi%my_proc+1, mod_parallel_3d%mpi%n_procs)

    call sendrecv(mod_parallel_3d%mpi, &
         dq(k), &
         left, 0, &
         last_q, &
         right, 0, &
         i)

#endif

    do i = k+1, p%natloc

       if (IS_EL(f, p, i)) then
          dr(j)   = -dq(last_i) + dq(i)
          last_i  = i
          j       = j+1
       endif

    enddo

#ifdef _MP
    if (mod_parallel_3d%mpi%my_proc /= mod_parallel_3d%mpi%n_procs-1) then
       dr(j)  = -dq(last_i) + last_q
       j      = j+1
    endif
#endif

!    AASSERT(naux == j-1, "naux == j-1", error)

  endsubroutine dq2aux


  !>
  !! Back-transformation
  !!
  !! Back-transformation
  !<
  subroutine aux2q(p, f, q0, naux, r, q, error)
    implicit none

    type(particles_t), intent(in)  :: p
    integer, intent(in)            :: f
    real(DP), intent(in)           :: q0
    integer, intent(in)            :: naux
    real(DP), intent(in)           :: r(naux)
    real(DP), intent(inout)        :: q(p%natloc)
    integer, intent(out), optional :: error

    ! ---

    integer   :: i, j, k

#ifdef _MP
    integer   :: left, right
    real(DP)  :: first_q
#endif

    ! ---

    INIT_ERROR(error)

    k = 1
    do while (k <= p%natloc .and. .not. IS_EL(f, p, k))
       k = k+1
    enddo

    if (k > p%natloc) then
       RAISE_ERROR("No element matching the filter found.", error)
    endif

#ifdef _MP

    left   = modulo(mod_parallel_3d%mpi%my_proc-1, mod_parallel_3d%mpi%n_procs)
    right  = modulo(mod_parallel_3d%mpi%my_proc+1, mod_parallel_3d%mpi%n_procs)

    call sendrecv(mod_parallel_3d%mpi,  &
         r(naux), &
         right, 0, &
         first_q, &
         left, 0, &
         i)

    if (mod_parallel_3d%mpi%my_proc == 0) then
       q(k) = -r(1)
    else
       q(k) = first_q - r(1)
    endif

#else

    q(k)    = -r(1)

#endif

    j       = 2

    do i = k+1, p%natloc

       if (IS_EL(f, p, i)) then
!#ifdef _MP
!          AASSERT((dmp_rank == dmp_numprocs-1 .and. j <= naux+1) .or. j <= naux)
!#else
!          AASSERT(j <= naux+1)
!#endif

#ifdef _MP
          if (mod_parallel_3d%mpi%my_proc == mod_parallel_3d%mpi%n_procs-1 .and. j == naux+1) then
#else
          if (j == naux+1) then
#endif
             q(i)  = r(j-1)
          else
             q(i)  = r(j-1) - r(j)
          endif
          j     = j+1
       endif

    enddo

    if (q0 /= 0.0_DP) then

       do i = 1, p%natloc
          if (IS_EL(f, p, i)) then
             q(i)  = q(i) + q0
          endif
       enddo

    endif

  endsubroutine aux2q


  !>
  !! Conjugate-gradients solver
  !!
  !! Compute equilibrated charges using a conjugate-gradients algorithm
  !<
  subroutine variable_charge_conjugate_gradients(this, p, nl, q, coul, epot, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(variable_charge_t), intent(inout) :: this
    type(particles_t),       intent(inout) :: p
    type(neighbors_t),       intent(inout) :: nl
    real(DP),                intent(inout) :: q(p%maxnatloc)
    type(C_PTR),             intent(in)    :: coul
    real(DP),                intent(inout) :: epot
    integer,       optional, intent(inout) :: error

    ! ---

    real(DP), allocatable  :: g(:), h(:), xi(:), dq(:), phi1(:), phi2(:), r(:)
    real(DP)               :: lambda, gg, dgg, gamma, E, previous_E, tot_q
    real(DP)               :: dmu, f, df

    integer                :: i, nit, naux, eli

    character(1)           :: updown_str

    ! ---

    call timer_start("variable_charge_conjugate_gradients")

    if (this%log) then
       call prlog("- variable_charge_conjugate_gradients -")
    endif

    !
    ! Extrapolate charges from previous time steps
    !

    call extrapolate(this%extrapolation, p, q, error=error)
    PASS_ERROR(error)

    !
    ! FIXME!!! Memory leak... Deallocate again
    !

    if (.not. allocated(g)) then
       allocate(g(p%maxnatloc))
       allocate(h(p%maxnatloc))
       allocate(xi(p%maxnatloc))
       allocate(r(p%maxnatloc))
       allocate(dq(p%maxnatloc))
       allocate(phi1(p%maxnatloc))
       allocate(phi2(p%maxnatloc))
    endif

    !
    ! Transform the charge distribution to the auxilliary space
    !

    call q2aux(p, this%f, this%q0, q, naux, r)

    if (this%trace) then
       write (ilog, '(5X,A5,A2,5A20)')  &
            "it", " ", "E", "d(mu)", "d(q)", "sum(q)", "d(E)"
    endif

    tot_q = sum(q(1:p%natloc))
#ifdef _MP
    call sum_in_place(mod_parallel_3d%mpi, tot_q)
#endif

    ! Back transformation, now sum(q) = this%q0
    call aux2q(p, this%f, this%q0, naux, r, q)
    !ASSERT(sum(q(1:p%natloc)) == 0, "sum(q(1:p%natloc)) == 0", error)
    call I_changed_other(p)

    dq  = 0.0_DP

    ! Get electrostatic potential phi1 and steepest descent direction dq
    phi1(1:p%natloc)  = 0.0_DP
    call coulomb_potential(coul, p, nl, q, phi1, error)
    PASS_ERROR(error)
    call add_X(this, p, phi1)

    do i = 1, p%natloc
       if (IS_EL(this%f, p, i)) then
          eli = p%el(i)
          dq(i) = phi1(i) + dE_V(q(i), this%V(eli), this%Vp(eli))
       endif
    enddo

    ! Some debug output
    if (this%trace) then
       E = 0.0_DP
       do i = 1, p%natloc
          if (IS_EL(this%f, p, i)) then
             eli = p%el(i)
             E = E + q(i)*phi1(i) + 2*E_V(q(i), this%V(eli), this%Vp(eli))
          else
             E = E + q(i)*phi1(i)
          endif
       enddo
       E = E/2
       previous_E = E

       tot_q = sum(q(1:p%natloc))
#ifdef _MP
       call sum_in_place(mod_parallel_3d%mpi, E)
       call sum_in_place(mod_parallel_3d%mpi, tot_q)
#endif

       write (ilog, '(12X,ES20.10,20X,2ES20.10)')  &
            E, maxval(abs(dq)), tot_q
    endif

    ! Transform the gradient to the plane in which sum q_i = 0 (sum xi_i = 0)
    call dq2aux(p, this%f, dq, naux, xi)
    g(1:naux)   = -xi(1:naux)
    h(1:naux)   = g(1:naux)
    xi(1:naux)  = h(1:naux)

    ! Initialize variables
    nit     = 1

    dmu  = this%convergence + 1.0_DP
    gg   = 1.0_DP
    do while (dmu > this%convergence .and. gg /= 0.0_DP)

       !**
       !* BEGIN LINE SEARCH IN DIRECTION xi (dq)
       !**

       ! Compute the pseudo electrostatic potential for the gradient dq E
       ! Back transformation -> sum(dq) = 0
       call aux2q(p, this%f, this%q0, naux, xi, dq)
       call I_changed_other(p)
       phi2(1:p%natloc)    = 0.0_DP
       call coulomb_potential(coul, p, nl, dq, phi2, error)
       PASS_ERROR(error)

       dmu = maxval(abs(dq))
#ifdef _MP
       dmu = max(mod_parallel_3d%mpi, dmu)
#endif

#if 0
          ! Nonlinear line search in direction phi1 (Newton's method)
          call f_and_df(p, dq, phi1, phi2, lambda, f, df)
   !       call print("Begin lambda iteration...")
          lambda = 0.0_DP ! Precondition?
   !       write (*, '(6ES20.10)')  lambda, f, df, maxval(abs(dq)), dot_product(dq, dq), dot_product(xi(1:naux), xi(1:naux))
          do while (abs(f) > 1d-12)
             if (df == 0.0_DP) then
                ! Some heuristics if the second derivative becomes zero
                lambda = 0.9*lambda
             else
                lambda = lambda - f/df
             endif
             call f_and_df(p, dq, phi1, phi2, lambda, f, df)
   !          write (*, '(3ES20.10)')  lambda, f, df
          enddo
       endif
#endif

       lambda = - sum(phi1*dq)/sum(phi2*dq)

       ! Update charges
       q(1:p%natloc)  = q(1:p%natloc) + lambda*dq(1:p%natloc)
       call I_changed_other(p)

       !**
       !* END LINE SEARCH
       !**

       ! Get electrostatic potential phi1 and new steepest descent direction dq
       phi1  = 0.0_DP
       call coulomb_potential(coul, p, nl, q, phi1, error)
       PASS_ERROR(error)
       call add_X(this, p, phi1)

!      This should be equal to 0
!       call f_and_df(p, dq, phi1, phi2, 0.0_DP, f, df)
!       write (*, '(ES20.10)')  f

       do i = 1, p%natloc
          if (IS_EL(this%f, p, i)) then
             eli = p%el(i)
             dq(i) = phi1(i) + dE_V(q(i), this%V(eli), this%Vp(eli))
          endif
       enddo

       ! Some debug output
       if (this%trace) then
          E = 0.0_DP
          do i = 1, p%natloc
             if (IS_EL(this%f, p, i)) then
                eli = p%el(i)
                E = E + q(i)*phi1(i) + 2*E_V(q(i), this%V(eli), this%Vp(eli))
             else
                E = E + q(i)*phi1(i)
             endif
          enddo
          E = E/2

          tot_q = sum(q(1:p%natloc))
#ifdef _MP
          call sum_in_place(mod_parallel_3d%mpi, E)
          call sum_in_place(mod_parallel_3d%mpi, tot_q)
#endif

          if (E < previous_E) then
             updown_str = 'v'
          else
             if (E > previous_E) then
                updown_str = '^'
             else
                updown_str = '='
             endif
          endif

          write (ilog, '(5X,I5,A2,5ES20.10)')  &
               nit, updown_str, E, dmu, maxval(abs(dq)), tot_q, E - previous_E

          previous_E = E
       endif

       ! Transform to the plane in which sum q_i = 0
       call dq2aux(p, this%f, dq, naux, xi)

       ! Conjugate gradients
       gg   = dot_product(g(1:naux), g(1:naux))
       dgg  = dot_product(xi(1:naux), xi(1:naux)) + dot_product(g(1:naux), xi(1:naux))  ! Polak-Ribiere

#ifdef _MP
       call sum_in_place(mod_parallel_3d%mpi, gg)
       call sum_in_place(mod_parallel_3d%mpi, dgg)
#endif

       if (gg /= 0.0_DP) then

          gamma = dgg/gg

          ! Get the gradient of the energy with respect to q
          g(1:naux)   = -xi(1:naux)
          h(1:naux)   = g(1:naux) + gamma*h(1:naux)   ! <--- conjugate gradients
          !             h(1:naux)   = g(1:naux)  ! <--- steepest descent
          xi(1:naux)  = h(1:naux)

       endif

       nit = nit + 1

    enddo

    if (this%trace) then

       dmu = maxval(abs(dq))
#ifdef _MP
       dmu = max(mod_parallel_3d%mpi, dmu)
#endif

       write (ilog, '(5X,5X,20X,2ES20.10)')  dmu, maxval(abs(dq))

    endif

    epot = epot + get_epot(this, p, q)

    if (this%log) then
       write (ilog, '(5X,I5,A)')  nit, " iterations to convergence."
       write (ilog, *)
    endif

    call timer_stop("variable_charge_conjugate_gradients")

  contains

    subroutine f_and_df(p, dq, phi1, phi2, lambda, f, df)
      implicit none

      type(particles_t), intent(in)  :: p
      real(DP), intent(in)           :: dq(p%natloc)     ! Gradient
      real(DP), intent(in)           :: phi1(p%natloc)   ! Electrostatic potential
      real(DP), intent(in)           :: phi2(p%natloc)   ! Pseudo electrostatic potential of the gradient
      real(DP), intent(in)           :: lambda
      real(DP), intent(out)          :: f
      real(DP), intent(out)          :: df

      ! ---

      integer   :: i, eli
      real(DP)  :: q2

      ! ---
      
      f   = 0.0_DP
      df  = 0.0_DP
      do i = 1, p%natloc
         if (IS_EL(this%f, p, i)) then
            eli  = p%el(i)
            q2   = q(i) + lambda*dq(i)

            f   = f + phi1(i)*dq(i) + lambda*phi2(i)*dq(i) + &
                 dq(i)*dE_V(q2, this%V(eli), this%Vp(eli))
            df  = df + phi2(i)*dq(i) + &
                 dq(i)**2*d2E_V(q2, this%V(eli), this%Vp(eli))
         endif
      enddo

    endsubroutine f_and_df

  endsubroutine variable_charge_conjugate_gradients


  !>
  !! Anderson solver
  !!
  !! Compute equilibrated charges by a self-consistent loop with
  !! Anderson mixing.
  !<
  subroutine variable_charge_anderson(this, p, nl, q, coul, epot, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(variable_charge_t), intent(inout) :: this     !< Variable charge object
    type(particles_t),       intent(inout) :: p                !< Particles
    type(neighbors_t),       intent(inout) :: nl               !< Neighbor list
    real(DP),                intent(inout) :: q(p%maxnatloc)   !< Charges
    type(C_PTR),             intent(in)    :: coul             !< Coulomb module
    real(DP),                intent(inout) :: epot           !< Potential energy
    integer,       optional, intent(inout) :: error           !< Error signals

    ! ---

    real(DP)  :: max_dmu, max_dq, lambda, E, previous_E, tot_q, beta, f, df

    integer   :: i, n, naux, nlambda
    integer   :: nit

    ! ---

    this%trace = .true.
    this%log = .true.

    call timer_start("variable_charge_anderson")

    if (this%log) then
       call prlog("- variable_charge_anderson -")
    endif

    !
    ! Allocate arrays
    !

    call log_memory_start("variable_charge_anderson")

    naux  = filter_count(this%f, p)  ! number of atoms matching filter and thus included in charge transfer

    if (.not. allocated(this%phi)) then
       if (this%log) then
          write (ilog, '(A)')  "- variable_charge_anderson -"
       endif

       call log_memory_start("variable_charge_anderson")

       allocate(this%phi(p%maxnatloc))
       call log_memory_estimate(this%phi, "phi")
    endif

    if (allocated(this%r_in)) then
       if (size(this%r_in) < naux) then
          deallocate(this%r_in)
          deallocate(this%r_out)
          
          if (this%log) then
             write (ilog, '(A)')  "- variable_charge_anderson -"
          endif

          call log_memory_start("variable_charge_anderson")
       endif
    endif

    if (.not. allocated(this%r_in)) then
       allocate(this%r_in(naux))
       allocate(this%r_out(naux))

       call log_memory_estimate(this%r_in, "r_in")
       call log_memory_estimate(this%r_out, "r_out")

       call log_memory_stop("variable_charge_anderson")

       if (this%log) then
          write (ilog, *)
       endif
    endif

    !
    ! Set total charge to requested value
    !

    n      = 0
    tot_q  = 0.0_DP
    do i = 1, p%natloc
       if (IS_EL(this%f, p, i)) then
          n      = n + 1
          tot_q  = tot_q + q(i)
       endif
    enddo

#ifdef _MP
    call sum_in_place(mod_parallel_3d%mpi, n)
    call sum_in_place(mod_parallel_3d%mpi, tot_q)
#endif

    if (abs(tot_q - this%total_charge) > 1e-12) then
       write (ilog, '(5X,A)')         "Adjusting charge to give the requested total charge."
       write (ilog, '(5X,A,F20.10)')  "tot_q  = ", tot_q
       write (ilog, '(5X,A,I20)')     "n      = ", n

       ! This need to run to nat, because the fast multipole solver uses the charges in this
       ! domain, independent of whether it's a ghost or not
       do i = 1, p%nat
          if (IS_EL(this%f, p, i)) then
             q(i)  = q(i) + (this%total_charge - tot_q)/n
          endif
       enddo

       call I_changed_other(p)
    endif

    !
    ! Get the potential
    !

    this%phi(1:p%nat)  = 0.0_DP   ! Note: This needs to be 1:p%nat, NOT 1:p%natloc for the fast-multipole 
    call coulomb_potential(coul, p, nl, q, this%phi, error)
    PASS_ERROR(error)
    call add_X(this, p, this%phi)

    !
    ! Take only treated atoms to r_in
    !

    call get_dmu(p, q, this%phi, max_dmu)

    call filter_pack(this%f, p, q, this%r_in)

    E        = 1.0_DP
    previous_E   = 0.0_DP                                  ! energy of previous step
    beta     = this%beta                               ! starting value for beta
#ifdef _MP
    lambda   = -filter_average(this%f, p, this%phi, &
         mpi=mod_parallel_3d%mpi)                      ! starting value for lambda
#else
    lambda   = -filter_average(this%f, p, this%phi)
#endif
    nlambda  = 0
    nit      = 0                                       ! current iteration
    max_dq   = 0.0_DP
    do while (abs(max_dmu) > this%convergence .and. nit < this%max_iterations)
       nit = nit+1

       ! Output energy and other stuff
       if (this%trace) then
          ! Energy
          E      = 0.0_DP
          tot_q  = 0.0_DP
          do i = 1, p%natloc
             if (IS_EL(this%f, p, i)) then
                E = E + q(i)*this%phi(i)/2 + E_V(q(i), this%V(p%el(i)), this%Vp(p%el(i)))
             else
                E = E + q(i)*this%phi(i)/2
             endif
             tot_q = tot_q + q(i)
          enddo

#ifdef _MP
          call sum_in_place(mod_parallel_3d%mpi, E)
          call sum_in_place(mod_parallel_3d%mpi, tot_q)
#endif

          if (nit == 1) then
             write (ilog, '(5X,A5,A6,2X,4A20,A12)')  &
                  "it", "beta", "E", "d(mu)", "d(q)", "sum(q)", "n(lambda)"
             write (ilog, '(5X,I5,F6.3,2X,ES20.10,20X,2ES20.10)')  &
                  nit, beta, E, max_dmu, tot_q
          else
             if (E > previous_E) then
                write (ilog, '(5X,I5,F6.3,A2,4ES20.10,I12)')  &
                     nit, beta, '^', E, max_dmu, max_dq, tot_q, nlambda
             else
                write (ilog, '(5X,I5,F6.3,A2,4ES20.10,I12)')  &
                     nit, beta, 'v', E, max_dmu, max_dq, tot_q, nlambda
             endif
          endif

          previous_E  = E

          flush(ilog)
       endif

       !
       ! Main contents of main loop
       !

       ! Compute Lagrange multiplier, Newton's method
       f = 1.0_DP
       nlambda = 0
       do while (abs(f) > 1d-6)
          call f_and_df(p, this%total_charge, q, this%V, this%Vp, this%phi, lambda, f, df)
          lambda  = lambda - f/df
          nlambda = nlambda + 1
       enddo

       ! New charges based on iterative formula
       max_dq  = 0.0_DP
       n       = 0
       do i = 1, p%natloc
          if (IS_EL(this%f, p, i)) then
             n  = n + 1
             ! Lambda is the Lagrange-multiplier for the total charge conservation constraint
             this%r_out(n) = next_q(q(i), this%V(p%el(i)), this%Vp(p%el(i)), this%phi(i), lambda)
             max_dq = max(max_dq, abs(q(i)-this%r_out(n)))
          endif
       enddo

       ! Anderson mixer
       call mix(this%mixer, nit, naux, this%r_in, this%r_out, beta)

       ! Charges back to particles object
       call filter_unpack(this%f, p, this%r_in, q)

#ifdef _MP
       call communicate_ghosts(mod_parallel_3d, p, .false.)
#endif

       ! Notify particles that charges changed
       call I_changed_other(p)

       ! Calculate new potential
       this%phi(1:p%nat)  = 0.0_DP   ! Note: This needs to be 1:p%nat, NOT 1:p%natloc for the fast-multipole solver
       call coulomb_potential(coul, p, nl, q, this%phi, error)
       PASS_ERROR(error)
       call add_X(this, p, this%phi)

       ! Update beta
       beta  = min(beta*this%beta_mul, this%beta_max)

       ! Get maximum variable in chemical potential
       call get_dmu(p, q, this%phi, max_dmu)
    enddo

    if (this%trace) then
       write (ilog, '(38X,ES20.10)') max_dmu
    endif

    !
    ! Add our contribution to potential energy
    !
    epot = epot + get_epot(this, p, q)

    if (this%log) then
       write (ilog, '(5X,I5,A)')  nit, " iterations to convergence."
       write (ilog, *)
    endif

    if (nit >= this%max_iterations) then
       call prscrlog("Charge equilibration did not converge: " // nit // " iterations exceeded.")
    endif

    call timer_stop("variable_charge_anderson")

  contains

    ! Compute the Lagrange multiplier and maximum derivative of E wrt. q_i
    subroutine get_dmu(p, q, phi, max_dE)
      implicit none

      type(particles_t), intent(in)  :: p
      real(DP), intent(in)           :: q(p%natloc)
      real(DP), intent(in)           :: phi(p%natloc)
      real(DP), intent(out)          :: max_dE

      ! ---

      real(DP)  :: mu(p%natloc), mu1, mu2
      integer   :: i

      ! ---

      max_dE  = 0.0_DP
      mu      = 0.0_DP
      do i = 1, p%natloc
         if (IS_EL(this%f, p, i)) then
            mu(i)  = phi(i) + dE_V(q(i), this%V(p%el(i)), this%Vp(p%el(i)))  ! dE/dq_i
!            write (*, '(I5,3ES20.10)') i, mu(i), q(i), phi(i)
         endif
      enddo

      mu1 = minval(mu, filter_mask(this%f, p))
      mu2 = maxval(mu, filter_mask(this%f, p))

#ifdef _MP
      mu1 = min(mod_parallel_3d%mpi, mu1)
      mu2 = max(mod_parallel_3d%mpi, mu2)
#endif

      max_dE = mu2 - mu1

    end subroutine get_dmu


    subroutine f_and_df(p, total_charge, q, V, Vp, phi, lambda, f, df)
      implicit none
 
      type(particles_t), intent(in)  :: p
      real(DP), intent(in)           :: total_charge
      real(DP), intent(in)           :: q(p%natloc)
      real(DP), intent(in)           :: V(p%nel)
      real(DP), intent(in)           :: Vp(p%nel)
      real(DP), intent(in)           :: phi(p%natloc)
      real(DP), intent(in)           :: lambda
      real(DP), intent(out)          :: f
      real(DP), intent(out)          :: df

      ! ---

      integer   :: i, eli
      real(DP)  :: xp, h
!      real(DP)  :: maxh

      ! ---

      f   = 0.0_DP
      df  = 0.0_DP

!      maxh = 0.0_DP

      do i = 1, p%natloc
         if (IS_EL(this%f, p, i)) then
            eli = p%el(i)

            xp  = (2.0_DP-Vp(eli))/(Vp(eli)-1.0_DP)

            h   = next_q(q(i), V(eli), Vp(eli), phi(i), lambda)
!            maxh  = max(maxh, abs(q(i)-h))
            f   = f + h
            if (abs(phi(i) + lambda) > 1d-12) then
               h   = 1.0_DP/(Vp(eli)-1.0_DP)
               h   = h / (V(eli)**h)
               df  = df + ( abs(phi(i) + lambda)**xp ) * h
            endif
         endif
      enddo

#ifdef _MP
      call sum_in_place(mod_parallel_3d%mpi, f)
      call sum_in_place(mod_parallel_3d%mpi, df)
#endif

!      write (*, *)  "maxh = ", maxh

      f   = total_charge - f

    endsubroutine f_and_df


    ! Next q for iteration
    pure function next_q(q, V, p, phi, lambda)
      implicit none

      real(DP), intent(in)  :: q
      real(DP), intent(in)  :: V
      real(DP), intent(in)  :: p
      real(DP), intent(in)  :: phi
      real(DP), intent(in)  :: lambda
      real(DP)              :: next_q

      ! ---

      real(DP)  :: h, xp

      ! ---

      h  = (phi + lambda)/V
      xp = 1.0_DP/(p-1.0_DP)

      if (h > 0.0_DP) then
         next_q = -(h**xp)
      else if (h < 0.0_DP) then
         next_q = ((-h)**xp)
      else
         next_q = 0.0_DP
      endif

    end function next_q

  endsubroutine variable_charge_anderson


  !>
  !! Car-Parrinello solver
  !!
  !! Propagate the charges using a Car-Parrinello type metadynamics
  !<
  subroutine variable_charge_car_parrinello(this, p, nl, q, coul, epot, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(variable_charge_t), intent(inout) :: this
    type(particles_t),       intent(inout) :: p
    type(neighbors_t),       intent(inout) :: nl
    real(DP),                intent(inout) :: q(p%maxnatloc)
    type(C_PTR),             intent(in)    :: coul
    real(DP),                intent(inout) :: epot
    integer,       optional, intent(inout) :: error

    ! ---

    real(DP)  :: lambda   !> Lagrange multiplier
    real(DP)  :: epot_ct, max_dE

    integer   :: i, n

    ! ---

    call timer_start("variable_charge_car_parrinello")

    if (.not. associated(this%qa)) then
       ! XXX FIXME!! This does not work if *p* changes
       call ptr_by_name(p%data, QA_STR, this%qa, ierror=error)
       PASS_ERROR(error)
    endif

    call coulomb_potential(coul, p, nl, q, this%qa, ierror=error)
    PASS_ERROR(error)
    call add_X(this, p, this%qa)

    epot_ct  = 0.0_DP
    lambda   = 0.0_DP
    n        = 0
    do i = 1, p%natloc
       if (IS_EL(this%f, p, i)) then
          epot_ct     = epot_ct + this%X(p%el(i))*q(i) + E_V(q(i), this%V(p%el(i)), this%Vp(p%el(i)))

          this%qa(i)  = this%qa(i) + dE_V(q(i), this%V(p%el(i)), this%Vp(p%el(i)))
          lambda      = lambda + this%qa(i)
          n           = n + 1
       endif
    enddo

#ifdef _MP
    call sum_in_place(mod_parallel_3d%mpi, lambda)
    call sum_in_place(mod_parallel_3d%mpi, n)
#endif

    lambda  = lambda/n
    max_dE  = 0.0_DP
    do i = 1, p%natloc
       if (IS_EL(this%f, p, i)) then
          this%qa(i)  = -( this%qa(i) - lambda )/this%mq
          if (abs(this%qa(i)) > max_dE) then
             max_dE  = abs(this%qa(i))
          endif
          this%qa(i)  = this%qa(i) - this%gamma*this%qv(i)
       endif
    enddo

#ifdef _MP
    call sum_in_place(mod_parallel_3d%mpi, max_dE)
#endif

    if (this%trace) then
       lambda  = sum(q(1:p%natloc))

#ifdef _MP
       call sum_in_place(mod_parallel_3d%mpi, lambda)
#endif

       write (ilog, '(A)')             "- variable_charge_car_parrinello -"
       write (ilog, '(5X,A,ES20.10)')  "max|dE|  = ", max_dE
       write (ilog, '(5X,A,ES20.10)')  "sum(q)   = ", lambda
       write (ilog, *)

       write (this%dE_un, '(2ES20.10)')  max_dE, lambda
    endif

    if (max_dE > this%dE_max) then
       this%n_fail  = this%n_fail + 1
    else
       epot         = epot + epot_ct
       this%n_fail  = 0
    endif

    if (this%n_fail >= this%fail_max) then
       ! Anderson mixing step if the criterion is exceeded
       call variable_charge_anderson(this, p, nl, q, coul, epot, error)
       PASS_ERROR(error)

       this%qv      = 0.0_DP
       this%qa      = 0.0_DP

       this%n_fail  = 0
    endif

    call timer_stop("variable_charge_car_parrinello")

  endsubroutine variable_charge_car_parrinello


  !>
  !! "Band-structure" energy
  !!
  !! Energy due to penalty term
  !<
  function variable_charge_get_epot(this, p, q) result(epot)
    implicit none

    type(variable_charge_t), intent(in)  :: this
    type(particles_t),       intent(in)  :: p
    real(DP),                intent(in)  :: q(p%maxnatloc)   !< Charges

    real(DP)                             :: epot

    ! ---

    integer  :: i, eli

    ! ---

    epot = 0.0_DP
    do i = 1, p%natloc
       if (IS_EL(this%f, p, i)) then
          eli = p%el(i)
          epot = epot + E_V(q(i), this%V(eli), this%Vp(eli))
       endif
    enddo

  endfunction variable_charge_get_epot



  subroutine add_X(this, p, phi)
    implicit none

    type(variable_charge_t), intent(in)  :: this
    type(particles_t), intent(in)        :: p
    real(DP), intent(inout)              :: phi(p%natloc)

    ! ---

    integer  :: i

    ! ---

    do i = 1, p%natloc
       if (IS_EL(this%f, p, i)) then
          phi(i)  = phi(i) + this%X(p%g(i))
       endif
    enddo

  endsubroutine add_X


  subroutine variable_charge_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(variable_charge_t), target, intent(inout) :: this
    type(c_ptr),                     intent(in)    :: cfg
    type(c_ptr),                     intent(out)   :: m

    ! ---

    integer, parameter  :: n_st                    = 4
    integer, parameter  :: len_solver_type_str     = 25

    character(len_solver_type_str), parameter  :: STR_conjugate_gradients  = CSTR("conjugate-gradients")
    character(len_solver_type_str), parameter  :: STR_anderson             = CSTR("anderson")
    character(len_solver_type_str), parameter  :: STR_car_parrinello       = CSTR("car-parrinello")
    character(len_solver_type_str), parameter  :: STR_disable              = CSTR("disable")
    character(len_solver_type_str), parameter  :: solver_type_strs(n_st)  = &
         (/ STR_conjugate_gradients, STR_anderson, STR_car_parrinello, STR_disable /)

    ! ---

    m = ptrdict_register_section(cfg, CSTR("VariableCharge"), &
         CSTR("A simple charge transfer model (see: A. K. Rappe and W. A. Goddard III, J. Phys. Chem. 95, 3358 (1991), F. H. Streitz and J. W. Mintmire, Phys. Rev. B 50, 11996 (1994))"))

    call ptrdict_register_string_property(m, c_loc(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements for which to enable charge transfer."))

    call ptrdict_register_enum_property(m, c_loc(this%solver_type), &
         n_st, len_solver_type_str, solver_type_strs(:), &
         CSTR("solver_type"), &
         CSTR("Type of solver to use ('conjugate-gradients', 'anderson' or 'disable'."))

    call ptrdict_register_real_property(m, c_loc(this%total_charge), &
         CSTR("total_charge"), &
         CSTR("Total charge on the system for which to enable charge transfer."))
         
    call ptrdict_register_real_property(m, c_loc(this%convergence), &
         CSTR("convergence"), &
         CSTR("Convergence criterium (on the maximum chemical potential, i.e. the gradient of the total energy)."))

    call ptrdict_register_integer_property(m, c_loc(this%extrapolation_memory), &
         CSTR("extrapolation_memory"), &
         CSTR("Number of past time steps to consider for charge extrapolation (minimum of 2, extrapolation is disabled if less)."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Frequency of charge equilibration."))

    call ptrdict_register_boolean_property(m, c_loc(this%log), CSTR("log"), &
         CSTR("Write number of iterations to log-file."))
    call ptrdict_register_boolean_property(m, c_loc(this%trace), CSTR("trace"), &
         CSTR("Trace self-consistency convergence."))

    call ptrdict_register_integer_property(m, c_loc(this%anderson_memory), &
         CSTR("anderson_memory"), &
         CSTR("Anderson mixing memory."))
    call ptrdict_register_real_property(m, c_loc(this%beta), CSTR("beta"), &
         CSTR("Mixing parameter."))
    call ptrdict_register_real_property(m, c_loc(this%beta_max), &
         CSTR("beta_max"), &
         CSTR("Maximum mixing parameter."))
    call ptrdict_register_real_property(m, c_loc(this%beta_mul), &
         CSTR("beta_mul"), &
         CSTR("Multiplicator for mixer increase."))

    call ptrdict_register_real_property(m, c_loc(this%mq), CSTR("mq"), &
         CSTR("Fictuous mass for Car-Parrinello type dynamics."))
    call ptrdict_register_real_property(m, c_loc(this%gamma), CSTR("gamma"), &
         CSTR("Dampening constant for Car-Parrinello type dynamics."))
    call ptrdict_register_real_property(m, c_loc(this%dE_max), CSTR("dE_max"), &
         CSTR("Upper limit of the derivative at which an Anderson mixing step is called."))
    call ptrdict_register_integer_property(m, c_loc(this%fail_max), &
         CSTR("fail_max"), &
         CSTR("Number of times the limit may be broken before and Anderson mixing step is called."))

    call ptrdict_register_integer_property(m, c_loc(this%max_iterations), &
         CSTR("max_iterations"), &
         CSTR("Maximum number of iterations."))

  endsubroutine variable_charge_register

endmodule variable_charge
