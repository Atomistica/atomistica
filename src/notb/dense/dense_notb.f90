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

! @meta
!   shared:directory
!   dependencies:../materials.f90,dense_hamiltonian_type.f90,dense_hamiltonian.f90,c_dense_hamiltonian.cpp,dense_forces.f90,dense_hs.f90,dense_repulsion.f90,solver/dense_occupation.f90,solver/dense_solver_lapack.f90,solver/dense_solver_cp.f90,solver/dense_solver_dispatch.f90,dense_scc.f90,analysis/dense_bonds.f90
!   classtype:dense_notb_t classname:TightBinding interface:potentials
! @endmeta

#include "macros.inc"
#include "filter.inc"

!>
!! The non-orthogonal tight-binding potential
!!
!! The non-orthogonal tight-binding potential
!!
!! See for example
!! M. Finnis, Interatomic Forces in Condensed Matter (2004).
!!
!<
module dense_notb
  use, intrinsic :: iso_c_binding

  use supplib

  use particles
  use neighbors
  use filter

#ifdef _MP
  use communicator
#endif

  use coulomb

  use materials

  use dense_hamiltonian

  use dense_force_calc
  use dense_hs
  use dense_repulsion

  use dense_solver

  use dense_scc

  use dense_bonds

  implicit none

  private

  public :: DENSE_NOTB_MAX_FOLDER_STRING

  character(*), parameter, private  :: MODULE_STR = "NOTB"

  integer, parameter                :: DENSE_NOTB_MAX_FOLDER_STRING = 1000

  public :: dense_notb_t
  type dense_notb_t

     !
     ! # of iterations
     !

     logical(BOOL)             :: enabled = .false.
     integer                   :: it

     !>
     !! Where to find the tight-binding database
     !<
     character(DENSE_NOTB_MAX_FOLDER_STRING)  :: database_folder = "*"

     !
     ! General stuff
     !

     !>
     !! Number of orbitals to populate
     !!
     !! Set this before bind_to is run to 
     !! set the number of orbitals to populate.
     !! If this is specified on input, qtot should not be,
     !! and vice versa.
     !<
     real(DP)                  :: noc  = 0.0_DP

     !>
     !! Total charge of TB subsystem
     !!
     !! Set this before bind_to is run to 
     !! set the total charge of the tight-binding subsystem.
     !! If this is specified on input, noc should not be,
     !! and vice versa.
     !!
     !! Negative qtot implies more electrons.
     !!
     !! The charges read in from the input file will be ignored if
     !! their sum doesn't match the desired value. If the charges
     !! are changed, a warning is printed. In this case, a homogeneous
     !! charge distribution is created.
     !!
     !! Note that the NOTB total charge can be set so that the total
     !! charge of the system is non-zero, but then things
     !! will probably not work always, so be careful!
     !<
     real(DP)                  :: qtot = 0.0_DP

     character(MAX_EL_STR)     :: elements  = "*"  !< Elements included in tight-binding
     integer                   :: el

#ifdef _MPI
     !> Cutoff for divide-and-conquer buffer zone
     !!
     !! Determine the cutoff for the divide-and-conquer buffer zone in
     !! distance units. If negative, the Hamiltonian cutoff is used as the
     !! buffer zone cutoff.
     !<
     real(DP)                  :: buffer_cutoff = -1.0_DP
#endif _MPI

     !
     ! Debug
     !

     logical(BOOL)             :: output_tables  = .false.

     !
     ! Tight-binding stuff
     !

     type(dense_hamiltonian_t)  :: tb                      !< Hamiltonian
     type(materials_t)          :: mat                     !< Materials database

     !
     ! Solver and SCC
     !

     type(dense_scc_t), pointer      :: scc     => NULL()
     type(dense_solver_t), pointer   :: solver  => NULL()

  endtype dense_notb_t


  public :: init
  interface init
     module procedure dense_notb_init
  endinterface

  public :: del
  interface del
     module procedure dense_notb_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure dense_notb_bind_to
  endinterface

  public :: energy_and_forces_with_charges
  interface energy_and_forces_with_charges
     module procedure dense_notb_energy_and_forces
  endinterface

  public :: set_solver
  interface set_solver
     module procedure dense_notb_set_solver
  endinterface

  public :: set_scc
  interface set_scc
     module procedure dense_notb_set_scc
  endinterface

  public :: set_Coulomb
  interface set_Coulomb
     module procedure dense_notb_set_Coulomb
  endinterface

  public :: get_dict
  interface get_dict
     module procedure dense_notb_get_dict
  endinterface

  public :: get_per_bond_property
  interface get_per_bond_property
     module procedure dense_notb_get_per_bond_property
  endinterface

  public :: register, dense_notb_register
  interface register
     module procedure dense_notb_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize the tight-binding potential
  !!   i.e.  read materials database
  !!
  !! Note that \a nl, \a solver and \a scc
  !! can be set a later stage using the \a set_neighbors, \a set_solver
  !! and \a set_scc methods.
  !! \sa set_neighbors
  !! \sa set_solver
  !! \sa set_scc
  !<
  subroutine dense_notb_init(this, solver, scc, elements, qtot, ierror)
    implicit none

    type(dense_notb_t), intent(inout)       :: this         !< NOTB object
    type(dense_solver_t), target, optional  :: solver       !< Solver object
    type(dense_scc_t), target, optional     :: scc          !< SCC object, if set, charge self-consistency is enabled
    character(*), intent(in), optional      :: elements
    real(DP), intent(in), optional          :: qtot
    integer, intent(out), optional          :: ierror       !< Error handling

    ! ---

    INIT_ERROR(ierror)

    ASSIGN_PROPERTY(elements)
    if (present(qtot)) then
       this%qtot = qtot
    endif

    ! init materials database
    call dense_notb_init_materials(this, ierror)
    PASS_ERROR(ierror)

    ! Standalone code. If scc is associated, but not enabled, erase it
    if (associated(this%scc)) then
       if (.not. this%scc%enabled) then
          deallocate(this%scc)
          this%scc => NULL()
       endif
    endif    

    ! Init dependent
    if (associated(this%scc)) then
       call init(this%scc, error=ierror)
       PASS_ERROR(ierror)
    endif

    if (present(solver)) then
       call set_solver(this, solver)
    else
       if (associated(this%solver)) then
          call init(this%solver, error=ierror)
          PASS_ERROR(ierror)
       else
          RAISE_ERROR("Please specify a solver.", ierror)
       endif
    endif

    if (present(scc)) then
       call set_scc(this, scc, ierror)
       PASS_ERROR(ierror)
    endif

  endsubroutine dense_notb_init


  !>
  !! Destructor
  !!
  !! Remove the tight-binding potential from memory, cleanup all allocated data
  !! structures.
  !<
  subroutine dense_notb_del(this)
    implicit none

    type(dense_notb_t), intent(inout)  :: this

    ! ---

    call del(this%tb)

    if (associated(this%solver)) then
       call del(this%solver)
       this%solver  => NULL()
    endif

    if (associated(this%scc)) then
       call del(this%scc)
       this%scc  => NULL()
    endif
    
  endsubroutine dense_notb_del


  !>
  !! Initialize materials database (internal)
  !!
  !! Initialize materials database (internal)
  !<
  subroutine dense_notb_init_materials(this, ierror)
    implicit none

    type(dense_notb_t), intent(inout)  :: this    !< NOTB object
    integer, optional, intent(out)     :: ierror  !< Error signals

    ! ---

    INIT_ERROR(ierror)

    if (trim(this%database_folder) == "*") then
       call read_database(this%mat, Hartree, Bohr, error=ierror)
    else
       call read_database(this%mat, Hartree, Bohr, folder=this%database_folder, error=ierror)
    endif
    PASS_ERROR(ierror)

    if (this%output_tables) then
       call write_tables(this%mat)
    endif

  end subroutine dense_notb_init_materials


  !>
  !! Initialize
  !!
  !! Initialize this NOTB object if all information is available
  !<
  subroutine dense_notb_bind_to(this, p, nl, ierror)
    implicit none

    type(dense_notb_t), target                   :: this
    type(particles_t),            intent(inout)  :: p
    type(neighbors_t),            intent(inout)  :: nl
    integer,            optional, intent(out)    :: ierror
    
    ! ---

    integer          :: i

    ! ---

    INIT_ERROR(ierror)
    call del(this%tb)

    ! Report
    call prlog("- dense_notb_bind_to -")
#ifdef COMPLEX_WF
    call prlog("The tight-binding module has been compiled for complex arithmetics.")
#else
    call prlog("The tight-binding module has been compiled for real arithmetics.")
#endif

    ! Init

    this%el  = filter_from_string(this%elements, p, ierror=ierror)
    PASS_ERROR(ierror)
    call init(this%tb, &
         this%mat, &
         p = p, &
         f = this%el, &
         error = ierror)
    PASS_ERROR(ierror)

    ! Request interaction range
    call request_interaction_range(nl, this%tb%cutoff)

#ifdef _MP
    ! Request border for buffer zone
    if (this%buffer_cutoff > 0.0_DP) then
       call request_border(mod_communicator, p, this%buffer_cutoff)
    endif
#endif

    ! Init other stuff

    this%it  = 0

    if (associated(this%scc)) then
       call bind_to(this%scc, p, this%tb, error=ierror)
       PASS_ERROR(ierror)

       call set_solver(this%scc, this%solver, error=ierror)
       PASS_ERROR(ierror)
    endif

  endsubroutine dense_notb_bind_to


  !>
  !! Set the solver object
  !!
  !! Set the solver object
  !<
  subroutine dense_notb_set_solver(this, solver)
    implicit none

    type(dense_notb_t), intent(inout)  :: this
    type(dense_solver_t), target       :: solver

    ! ---

    this%solver  => solver

    if (associated(this%scc)) then
       call set_solver(this%scc, this%solver)
    endif

  endsubroutine dense_notb_set_solver


  !>
  !! Set the SCC object
  !!
  !! SCC is only enabled after this function has been called.
  !<
  subroutine dense_notb_set_scc(this, scc, ierror)
    implicit none

    type(dense_notb_t), intent(inout)  :: this    !< NOTB object
    type(dense_scc_t), target          :: scc     !< SCC object
    integer, intent(out), optional     :: ierror  !< Error signals

    ! ---

    INIT_ERROR(ierror)

    this%scc  => scc

    if (associated(this%solver)) then
       call set_solver(this%scc, this%solver)
    endif

  endsubroutine dense_notb_set_scc


  !>
  !! Set the Coulomb solver
  !!
  !! Set the Coulomb solver, i.e. pass it to SCC
  !<
  subroutine dense_notb_set_Coulomb(this, coul, ierror)
    use, intrinsic :: iso_c_binding

    implicit none

    type(dense_notb_t), intent(inout) :: this    !< NOTB object
    type(C_PTR),        intent(in)    :: coul    !< SCC object
    integer,  optional, intent(out)   :: ierror  !< Error signals

    ! ---

    INIT_ERROR(ierror)

    if (associated(this%scc)) then
       call set_Coulomb(this%scc, coul, ierror)
       PASS_ERROR(ierror)
    endif

  endsubroutine dense_notb_set_Coulomb


  !>
  !! Return additional dependent per bond properties
  !!
  !! Return additional dependent per bond properties
  !<
  subroutine dense_notb_get_per_bond_property(this, p, nl, propstr, propout, error)
    implicit none

    type(dense_notb_t), intent(inout) :: this        !< NOTB object
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    character(*),       intent(in)    :: propstr     !< Name of property to return
    real(DP),           intent(out)   :: propout(*)  !< Return buffer, needs to have same length as neighbor list object
    integer,  optional, intent(out)   :: error       !< Error signals

    ! ---

    INIT_ERROR(error)

    if (trim(propstr) == "overlap_population") then
       call bond_analysis(this%tb, p, nl, overlap_population=propout)
    else
    if (trim(propstr) == "Loewdin_bond_order") then
       call bond_analysis(this%tb, p, nl, Loewdin_bond_order=propout)
    else
    if (trim(propstr) == "covalent_bond_energy") then
       call bond_analysis(this%tb, p, nl, e_cov=propout)
    else
       RAISE_ERROR("Unknown bond property '" // propstr // "'.", error)
    endif
    endif
    endif

  endsubroutine dense_notb_get_per_bond_property


  !>
  !! Return dictionary object containing pointers to internal data
  !<
  subroutine dense_notb_get_dict(this, dict, error)
    implicit none

    type(dense_notb_t), intent(inout) :: this        !< NOTB object
    type(ptrdict_t),    intent(inout) :: dict
    integer,  optional, intent(out)   :: error       !< Error signals

    ! ---

    INIT_ERROR(error)

    call get_dict(this%tb, dict, error)
    PASS_ERROR(error)
    if (associated(this%solver)) then
       call get_dict(this%solver, dict, error)
       PASS_ERROR(error)
    endif

  endsubroutine dense_notb_get_dict


  !>
  !! Compute the energies and forces
  !!
  !! Compute the energies and forces
  !<
  subroutine dense_notb_energy_and_forces(this, p, nl, epot, f, wpot, q, &
       epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
       ierror)
    implicit none

    type(dense_notb_t), intent(inout) :: this
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(inout) :: epot
    real(DP),           intent(inout) :: f(3, p%maxnatloc)
    real(DP),           intent(inout) :: wpot(3, 3)
    real(DP), optional, intent(inout) :: q(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_at(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_bond(nl%neighbors_size)
    real(DP), optional, intent(inout) :: f_per_bond(3, nl%neighbors_size)
#ifdef LAMMPS
    real(DP), optional, intent(inout) :: wpot_per_at(6, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot_per_bond(6, nl%neighbors_size)
#else
    real(DP), optional, intent(inout) :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot_per_bond(3, 3, nl%neighbors_size)
#endif
    integer,  optional, intent(out)   :: ierror

    ! ---

    real(DP) :: noc

    ! ---

    INIT_ERROR(ierror)

    ! Verify pointers
    if(.not. associated(this%solver)) then
       RAISE_ERROR("No solver object associated to NOTB.", ierror)
    end if

    call timer_start("dense_notb_energy_and_forces")

    ! Update neighbor list
    call update(nl, p, ierror)
    PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)

    ! Update the atom to H/S matrix-block match
    call update_orbitals(this%tb, p, error=ierror)
    PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)

    call hs_setup(this%tb, this%mat, p, nl, error=ierror)
    PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)

    ! Compute number of occupied orbitals
    noc = get_occupied_orbitals(p, this%el, this%tb, this%qtot)

    ! Solve either with or without charge self-consistency
    if (associated(this%scc)) then

       if (.not. present(q)) then
          RAISE_ERROR_AND_STOP_TIMER("Please provide a charge-array for the self-consistent solution.", "dense_notb_energy_and_forces", ierror)
       endif

       this%it  = this%it + 1
       call establish_self_consistency( &
            this%scc, p, nl, this%tb, q, &
            noc, &
            f     = this%tb%f, &
            error = ierror)
       PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)
    else
       call diag_start(this%solver, this%tb, error=ierror)
       PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)
       call diag_HS(this%solver, this%tb, noc, error=ierror)
       PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)
       call diag_stop(this%solver, this%tb, error=ierror)
       PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)
    endif

    ! Calculate forces and energies
    !       Note: The Coulomb energy and forces are calculated from
    !             within the Coulomb module, which *has* to be loaded
    !             for SCC calculations and called separately.
    !   e_bs     = Band-structure energy
    !   e_rep    = Repulsive energy
    !   e_atomic = Energy of the individual charge neutral atoms
    call forces(p, nl, this%tb, this%mat, f, wpot, wpot_per_bond, error=ierror)
    PASS_ERROR_AND_STOP_TIMER("dense_notb_energy_and_forces", ierror)
    epot  = epot + e_bs(this%solver, this%tb) + &
         e_rep(this%tb, this%mat, p, nl) - e_atomic(this%tb, p)

    call timer_stop("dense_notb_energy_and_forces")

  endsubroutine dense_notb_energy_and_forces


  !>
  !!
  !!
  !!
  !<
  subroutine dense_notb_adjust_total_charge(p, f, q, qtot)
    implicit none

    type(particles_t), intent(in)  :: p               !< Particles
    integer, intent(in)            :: f               !< Filter for atom types
    real(DP), intent(inout)        :: q(p%maxnatloc)  !< Charges
    real(DP), intent(in)           :: qtot            !< Total charge to be set

    ! ---

    integer                        :: i               ! loops
    real(DP)                       :: qc              ! current charge
    integer                        :: n               ! atoms matching filter

    ! ---

    ! check total charge
    qc = 0.0_DP
    n = 0
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          n = n + 1
          qc = qc + q(i)
       end if
    end do

#ifdef _MPI
    call sum_in_place(mod_communicator%mpi, qc)
#endif

    if(abs(qc-qtot) > 1e-10) then
       WARN("Adjusting charge of TB (sub)system from " // qc // " to " // qtot // ". Using homogeneous charge distribution.")
       do i = 1, p%natloc
          if (IS_EL(f, p, i)) then
             q(i) = qtot/n
          end if
       end do
    end if

  end subroutine dense_notb_adjust_total_charge



  !>
  !! Number of occupied orbitals
  !!
  !! Calculates the total number of occupied orbitals from the elements in the system
  !! (to form a neutral molecule (occ=0), anion (occ=-1), or cation (occ=-2) )
  !<
  function get_occupied_orbitals(p, f, tb, q) result(noc)
    implicit none

    type(particles_t),         intent(in) :: p    !< Particles object
    integer,                   intent(in) :: f    !< Filter for atom types
    type(dense_hamiltonian_t), intent(in) :: tb   !< Atom type data from materials database
    real(DP),                  intent(in) :: q    !< Total charge

    real(DP)                              :: noc  !< Number of occupied orbitals

    ! ---

    real(DP)  :: occ
    integer   :: i

    ! ---

    type(notb_element_t), pointer :: at(:)

    ! ---

    call c_f_pointer(tb%at, at, [p%nat])

    occ = 0.0_DP
    do i = 1, p%nat
       if (IS_EL(f, p, i)) then
          occ = occ + at(i)%q0
       endif
    enddo

!#ifdef _MPI
!    call sum_in_place(mod_communicator%mpi, occ)
!#endif

    noc = (occ - q)/2

  endfunction get_occupied_orbitals


  !>
  !! Total charge
  !!
  !! Calculates the total charge from the elements in the system  
  !<
  function get_total_charge(p, f, at, noc) result(q)
    implicit none

    type(particles_t), intent(in)     :: p          !< Particles object
    integer, intent(in)               :: f          !< Filter for atom types
    type(notb_element_t), intent(in)  :: at(p%nat)  !< Atom type data from materials database
    real(DP), intent(in)              :: noc        !< Number of occupied orbitals

    real(DP)                          :: q

    ! ---

    real(DP)  :: occ
    integer   :: i

    ! ---

    occ = 0.0_DP
    do i = 1, p%natloc
       if (IS_EL(f, p, i)) then
          occ = occ + at(i)%q0
       endif
    enddo

    q = -(2*noc - occ)

  endfunction get_total_charge


  ! --- REGISTRY ---

  subroutine dense_notb_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(dense_notb_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)               :: cfg
    type(c_ptr), intent(out)              :: m

    ! ---

    m = ptrdict_register_module(cfg, c_loc(this%enabled), CSTR("TightBinding"), &
         CSTR("Non-orthogonal tight-binding potential."))

    call ptrdict_register_string_property(m, c_locs(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements for which to activate this module."))

    call ptrdict_register_real_property(m, c_loc(this%noc), &
         CSTR("occupied_orbitals"), &
         CSTR("Number of occupied orbitals (automatically determined if set to zero)."))

    call ptrdict_register_boolean_property(m, c_loc(this%output_tables), &
         CSTR("output_tables"), &
         CSTR("Debug: Output Slater-Koster and repulsion tables to file."))

    call ptrdict_register_string_property(m, c_locs(this%database_folder), &
         DENSE_NOTB_MAX_FOLDER_STRING, CSTR("database_folder"), &
         CSTR("Folder containing the NOTB parametrization."))

#ifdef _MPI
    call ptrdict_register_real_property(m, c_loc(this%buffer_cutoff), &
         CSTR("buffer_cutoff"), &
         CSTR("Cutoff (width) of the divide-and-conquer buffer zone. Default: Same as Hamiltonian cutoff"))
#endif


    allocate(this%solver)
    call register(this%solver, m)
    allocate(this%scc)
    call register(this%scc, m)

  endsubroutine dense_notb_register

endmodule dense_notb
