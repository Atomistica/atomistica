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
! @endmeta

#include "macros.inc"

program main
  use supplib

  use particles
  use neighbors
  use dynamics

  use potentials
  use integrators
  use coulomb
  use callables

  use native_io

  use atomistica

  use signal_handler

#ifdef _MP
  use communicator
#endif

  implicit none

  character(*), parameter :: MODULE_STR = "main"

  ! Registries from *_factory_c.c
  interface

    subroutine integrators_factory_register(ptrdict_file, this) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value  :: ptrdict_file
      type(c_ptr), value  :: this
    endsubroutine integrators_factory_register

    subroutine potentials_factory_register(ptrdict_file, this) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value  :: ptrdict_file
      type(c_ptr), value  :: this
    endsubroutine potentials_factory_register

    subroutine coulomb_factory_register(ptrdict_file, this) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value  :: ptrdict_file
      type(c_ptr), value  :: this
    endsubroutine coulomb_factory_register

    subroutine callables_factory_register(ptrdict_file, this) bind(C)
      use, intrinsic :: iso_c_binding
      type(c_ptr), value  :: ptrdict_file
      type(c_ptr), value  :: this
    endsubroutine callables_factory_register

  endinterface

  !
  ! Stuff read from the simulation control file
  !

  integer, target            :: scr_freq  = 10
  integer, target            :: file_freq = 10

  integer, target            :: n_iterations = -1
  real(DP), target           :: max_time = -1
  real(DP), target           :: max_dt = 0.01
  real(DP)                   :: t0   = 0d0
  real(DP), target           :: max_dr = -1
  integer, target            :: max_nat = -1
  integer, target            :: avgn = 20

  real(DP), target           :: bin_size = -1

  real(DP), target           :: cutoff_add = -1.0_DP
  real(DP), target           :: abs_cutoff = -1.0_DP
  logical(BOOL), target      :: force_binning = .false.

  logical(BOOL), target      :: in_clear_velocities  = .false.
  logical(BOOL), target      :: in_clear_charges     = .false.

  integer, target            :: pbc(3) = (/ 1, 1, 1 /)

  logical(BOOL), target      :: in_sort = .false.

  logical(BOOL), target      :: center_molecule = .false.

  logical(BOOL), target      :: per_atom_virial = .false.
  
  real(DP), target           :: shear_gamma(3)
  real(DP), target           :: shear_tau
  real(DP), target           :: shear_amplitude(3)

  real(DP), pointer          :: q(:)  => NULL()

  !
  ! Initial transformations
  !

  real(DP), target           :: translate_by(3)

  ! ---

  call main_loop

contains

  !**********************************************************************
  ! The main subroutine
  !**********************************************************************    
  subroutine main_loop()
    use, intrinsic :: iso_c_binding

#ifdef HAVE_IFPORT
    use ifport
#endif

    implicit none

    character(3), parameter   :: month(12) = &
         (/ "Jan", "Feb", "Mar", "Apr", &
            "May", "Jun", "Jul", "Aug", &
            "Sep", "Oct", "Nov", "Dec" /)

    character(10), parameter :: out_filename(2) = &
         [ "atomsA.out", "atomsB.out" ]

    !
    ! General stuff: Particle information, binning, neighbor list
    !

    type(particles_t),   target :: p                 ! Particle information

    type(dynamics_t)            :: dyn

    type(neighbors_t)           :: nl                ! Pair information

    integer                     :: it, i, done_file
    integer                     :: ierror  = ERROR_NONE

    logical                     :: user_exit

    ! ---

#ifdef _MP
    type(MPI_Context) :: mpi
#endif

    type(integrators_t), target :: ints
    type(potentials_t),  target :: pots
    type(coulomb_t),     target :: coul
    type(callables_t),   target :: cals

    ! ---

    !
    ! Read input files
    !

    call read_ptrdict_file(p, c_loc(ints), c_loc(pots), c_loc(coul), c_loc(cals))

#ifdef _MP
    call initialise(mpi)
#endif

    call atomistica_startup

    !
    ! Initialize units
    !

    call units_init(system_of_units)

    !
    ! Initialize particles object and callables
    !

    call init(p)
    call init(dyn, p, dt = max_dt, mymaxtime = max_time)
    p%pbc          = pbc
    p%locally_pbc  = p%pbc /= 0

    if (coulomb_is_enabled(c_loc(coul))) then
       call add_real(p%data, Q_STR, Q_TAG, ierror=ierror)
       HANDLE_ERROR(ierror)
       call coulomb_register_data(c_loc(coul), p)
    endif

    if (per_atom_virial) then
       call add_real3x3(p%data, "virial", F_TO_TRAJ, ierror=ierror)
       HANDLE_ERROR(ierror)
    endif

    DEBUG_WRITE("- integrator_init -")
    call integrators_init(c_loc(ints), p)

    call potentials_init(c_loc(pots))
    call potentials_register_data(c_loc(pots), p, ierror)
    HANDLE_ERROR(ierror)
    call coulomb_init(c_loc(coul))

    call potentials_set_Coulomb(c_loc(pots), c_loc(coul), ierror)
    HANDLE_ERROR(ierror)

    call callables_init(c_loc(cals))
    call callables_register_data(c_loc(cals), p)

    if (max_nat > 0) then
       call allocate(p, max_nat, allow_def=.true., error=ierror)
       HANDLE_ERROR(ierror)
    endif

    ! Cell information is need for initialization of the domain decomposition module
    DEBUG_WRITE("- read_cell_from_atoms -")
    call read_cell_from_atoms(p, "atoms.dat", allow_def=.true.)
    if (centered_box) then
       p%upper  = p%upper/2
       p%lower  = -p%upper

       p%lower_with_border  = p%lower
       p%upper_with_border  = p%upper
    endif

#ifdef _MP
    if (abs_cutoff > 0.0_DP) then
       stop "abs_cutoff > 0.0 and parallel computation. Does not work yet. Please implement."
    endif
    call init(mod_communicator, p, verlet_shell = cutoff_add, error = ierror)
    HANDLE_ERROR(ierror)
#endif

    !
    ! Read particle positions
    !

    DEBUG_WRITE("- read_atoms -")
    if (allocated(ints%settle)) then
       call read_atoms(p, "atoms.dat", &
            mol        = ints%settle%molecules, &
            skip_cell  = .true., &
            allow_def  = .true.)
    else
       call read_atoms(p, "atoms.dat", &
            skip_cell  = .true., &
            allow_def  = .true.)
    endif

    p%dof  = 3*p%natloc
    do i = 1, p%natloc
       if (p%g(i) <= 0) then
          p%dof  = p%dof - 3
       endif
    enddo

#ifdef _MP
    DEBUG_WRITE("- computing dof -")
    call sum_in_place(mod_communicator%mpi, p%dof)
#endif

    p%dof = p%dof-3

    call print_to_log(p%data)

    !
    ! Initial transformations
    !

    if (center_molecule) then
       call center(p, cell=p%Abox)
    endif

    if (any(translate_by /= 0.0_DP)) then
       call prlog
       call prlog("     translate_by  = ( "//translate_by//" )")
       call prlog

!       call group_rigid_objects(p)
       call translate(p, translate_by)

       PNC3(p, :)  = POS3(p, :)
       PCN3(p, :)  = POS3(p, :)

       translate_by  = 0.0_DP
    endif

    call inbox(p)

    !
    ! Initialize integrators
    !

    call prlog("     degrees-of-freedom  = "//p%dof)
    call prlog

    !
    ! Initialize binning/neighbor lists
    !

    DEBUG_WRITE("- init(nl) -")
    call init(nl, &
         avgn          = avgn, &
         cutoff        = abs_cutoff, &
         verlet_shell  = cutoff_add, &
         sort          = logical(in_sort))
    if (bin_size > 0.0_DP) then
       call set(nl, bin_size = bin_size)
    endif

#ifdef _MP
    call allocate(mod_communicator, p)
#endif

    !
    ! Bind potentials to the Particles and Neighbor list objects
    !

    ! Note: potentials_bind_to calls Coulomb modules!
    call potentials_bind_to(c_loc(pots), p, nl, c_loc(coul), ierror=ierror)
    HANDLE_ERROR(ierror)
    call callables_bind_to(c_loc(cals), p, nl, c_loc(pots), c_loc(coul), ierror=ierror)
    HANDLE_ERROR(ierror)

#ifdef _MP
    call update(dyn, advance_time=.false., mpi=mod_communicator%mpi)
#else
    call update(dyn, advance_time=.false.)
#endif
    t0  = dyn%ti

    !
    ! Main loop
    !
    
#ifdef _MP
    if (mod_communicator%mpi%my_proc == ROOT) then
#endif
    write (*, *)
#ifdef _MP
    endif
#endif

    call prlog
    call prlog
    call prlog("====> ENTERING MAIN LOOP <====")
    call prlog

    it = 0
    done = .false.

#ifdef HAVE_IFPORT
    if (signal(SIGTERM, handle_signal, -1) == SIG$ERR) then
       WARN("Failed to trap SIGTERM.")
    endif
    if (signal(12, handle_signal, -1) == SIG$ERR) then  ! SIGUSR2
       WARN("Failed to trap SIGUSR2.")
    endif
#endif

    if (n_iterations == 0) then
       it = -1
       dyn%v  = 0.0_DP
       dyn%f  = 0.0_DP
    endif

    if (in_clear_velocities) then
       dyn%v  = 0.0_DP
    endif

!    if (exists(p%data, Q_STR)) then
    if (coulomb_is_enabled(c_loc(coul))) then
       call ptr_by_name(p%data, Q_STR, q)

       if (in_clear_charges) then
          q  = 0.0_DP
       endif
    endif

    if (per_atom_virial) then
       call ptr_by_name(p%data, "virial", pots%wpot_per_at)
    endif

    call timer_start("main loop")

    ! Compute forces for initial Verlet step
    ! Note: potentials_energy_and_forces calls Coulomb modules!
    DEBUG_WRITE("- potentials_force -")
    if (associated(q)) then
       call potentials_energy_and_forces(c_loc(pots), dyn, nl, coul=c_loc(coul), q=q, &
            ierror=ierror)
    else
       call potentials_energy_and_forces(c_loc(pots), dyn, nl, ierror=ierror)
    endif
    HANDLE_ERROR(ierror)

#ifdef _MP
    call update(dyn, advance_time=.false., mpi=mod_communicator%mpi)
#else
    call update(dyn, advance_time=.false.)
#endif

    DEBUG_WRITE("- callables_invoke -")
    call callables_invoke(c_loc(cals), dyn, nl, c_loc(pots), c_loc(coul), ierror)
    HANDLE_ERROR(ierror)

    do while ((n_iterations < 0 .or. dyn%it < n_iterations) .and. (max_time < 0.0_DP .or. dyn%ti < max_time) .and. .not. done)
       it = it + 1


       !
       ! Check if we need to checkpoint our run
       !

       if (it == 1 .or. mod(it, file_freq) == 0) then
          call write_atoms(p, out_filename(mod(it/file_freq, 2)+1), error=ierror)
          HANDLE_ERROR(ierror)
       endif


       !
       ! Integrator: Step 1
       !

       DEBUG_WRITE("- Integrator (1) -")
       call integrators_step1(c_loc(ints), c_loc(pots), dyn, max_dt, max_dr, &
            ierror=ierror)
       HANDLE_ERROR(ierror)
       if (shear_tau > 0.0_DP) then
          p%shear_dv  = shear_amplitude * sin(2*PI*(dyn%ti-t0)/shear_tau)
       endif
       if (any(shear_gamma /= 0.0_DP)) then
          call set_lees_edwards(p, p%shear_dx + shear_gamma*p%Abox(3, 3)*dyn%dt, dv=shear_gamma*p%Abox(3, 3))
       else if (any(p%shear_dv /= 0.0_DP)) then
          call set_lees_edwards(p, p%shear_dx + p%shear_dv*dyn%dt, dv=p%shear_dv)
       endif


       !
       ! Compute potential energy and forces
       !

       DEBUG_WRITE("- potentials_force -")
       if (associated(q)) then
          call potentials_energy_and_forces(c_loc(pots), dyn, nl, coul=c_loc(coul), &
               q=q, ierror=ierror)
       else
          call potentials_energy_and_forces(c_loc(pots), dyn, nl, ierror=ierror)
       endif
       HANDLE_ERROR(ierror)


       !
       ! Integrator: Step 2
       !

       DEBUG_WRITE("- Integrator (2) -")
       call integrators_step2(c_loc(ints), dyn, ierror=ierror)
       HANDLE_ERROR(ierror)

#ifdef _MP
       call update(dyn, it, mpi=mod_communicator%mpi)
#else
       call update(dyn, it)
#endif

       DEBUG_WRITE("- callables_invoke -")
       call callables_invoke(c_loc(cals), dyn, nl, c_loc(pots), c_loc(coul), ierror)
       HANDLE_ERROR(ierror)

       !
       ! Are we optimizing? Abort if convergence criterion has been reached.
       !

       if (allocated(ints%fire)) then
          done = .true.
          do i = 1, p%natloc
             if (p%g(i) > 0) then
                if (sqrt(dot_product(VEC3(dyn%f, i), VEC3(dyn%f, i))) > &
                     ints%fire%fmax) then
                   done = .false.
                endif
             endif
          enddo
       endif

       !
       ! Print to screen
       !

       if (mod(it, scr_freq) == 0) then
          
          call print_status(dyn)

          inquire(file="EXIT", exist=user_exit)
          if (user_exit) then
             call prscrlog("FOUND *EXIT*; USER REQUESTED ABORT")
             done = .true.
          endif
       endif

#ifdef _MP
       !
       ! This barrier is necessary for user requested aborts
       !

       call barrier(mod_communicator%mpi, ierror)
       HANDLE_ERROR(ierror)
#endif

    enddo

    call timer_stop("main loop")

    !
    ! Cleanup
    !

    call print_status(dyn)

    call prlog

    call write_atoms(p, "atoms.out", error=ierror)
    HANDLE_ERROR(ierror)

    call integrators_del(c_loc(ints))
    call potentials_del(c_loc(pots))
    call coulomb_del(c_loc(coul))
    call callables_del(c_loc(cals))

    call del(dyn)
    call del(nl)
#ifdef _MP
    call del(mod_communicator)
#endif
    call del(p)

#ifdef HAVE_IFPORT
    if (signal(SIGTERM, handle_signal, 0) == SIG$ERR) then
       WARN("Failed to untrap SIGTERM.")
    endif
    if (signal(12, handle_signal, 0) == SIG$ERR) then  ! SIGUSR2
       WARN("Failed to untrap SIGUSR2.")
    endif
#endif

    call atomistica_shutdown

#ifdef _MP
    call finalise(mpi)
#endif

    !
    ! Create the (empty) file DONE to let everyone know we finished properly
    !

    done_file = fopen("DONE", mode=F_WRITE)
    call fclose(done_file)

  endsubroutine main_loop


  !**********************************************************************
  ! Read the configuration file
  !**********************************************************************    
  subroutine read_ptrdict_file(p, ints, pots, coul, cals)
    use, intrinsic :: iso_c_binding

    implicit none

    type(particles_t), target     :: p
    type(C_PTR),       intent(in) :: ints
    type(C_PTR),       intent(in) :: pots
    type(C_PTR),       intent(in) :: coul
    type(C_PTR),       intent(in) :: cals


    ! ---

    type(c_ptr)  :: ptrdict_file, m

    ! ---

    ptrdict_file = ptrdict_register_section(C_NULL_PTR, CSTR("Simulation"), &
         CSTR("MD simulation control file."))

    call ptrdict_register_enum_property(ptrdict_file, c_loc(system_of_units), &
         n_units, len_unit_str, unit_strs, &
         CSTR("system_of_units"), &
         CSTR("'eV/A' or 'H/Bohr'"))

    call ptrdict_register_boolean_property(ptrdict_file, c_loc(centered_box), &
         CSTR("centered_box"), &
         CSTR("If true, the simulation box will be centered around the origin. Otherwise the origin will be the lower left corner."))

    call ptrdict_register_boolean_property(ptrdict_file, c_loc(center_molecule), &
         CSTR("center_molecule"), &
         CSTR("Center the initial structure within the simulation cell."))

    call ptrdict_register_boolean_property(ptrdict_file, c_loc(per_atom_virial), &
         CSTR("per_atom_virial"), &
         CSTR("Compute per atom virial information to be used with e.g. the Slicing module."))

    call ptrdict_register_intpoint_property(ptrdict_file, c_loc(pbc(1)), &
         CSTR("pbc"), &
         CSTR("Periodicity in x-, y- and z-direction."))

    call ptrdict_register_integer_property(ptrdict_file, c_loc(scr_freq), &
         CSTR("scr_freq"), &
         CSTR("Screen output interval."))
    call ptrdict_register_integer_property(ptrdict_file, c_loc(file_freq), &
         CSTR("file_freq"), &
         CSTR("File output interval (=> traj.xyz)."))

    call ptrdict_register_real_property(ptrdict_file, c_loc(max_dt), CSTR("dt"), &
         CSTR("Time step."))
    call ptrdict_register_real_property(ptrdict_file, c_loc(max_dr), CSTR("max_dr"), &
         CSTR("Maximum displacement (adaptive time stepping enabled if greater than 0)."))
    call ptrdict_register_integer_property(ptrdict_file, c_loc(n_iterations), &
         CSTR("n_iterations"), &
         CSTR("Maximum number of iterations."))
    call ptrdict_register_real_property(ptrdict_file, c_loc(max_time), &
         CSTR("max_time"), &
         CSTR("When to stop simulation (time)."))

    call ptrdict_register_integer_property(ptrdict_file, c_loc(max_nat), &
         CSTR("max_nat"), &
         CSTR("Size of internal particles object."))

    call ptrdict_register_real_property(ptrdict_file, c_loc(bin_size), &
         CSTR("bin_size"), &
         CSTR("Cutoff used for the binning routing (same as the interaction cut-off if < 0)."))
    call ptrdict_register_integer_property(ptrdict_file, c_loc(avgn), CSTR("avgn"), &
         CSTR("Average number of neighbors (for neighbor list initialization)."))

    call ptrdict_register_real_property(ptrdict_file, c_loc(abs_cutoff), &
         CSTR("cutoff"), &
         CSTR("Specify cutoff. Shell larger than the interaction range will be treated as a skin depth/Verlet shell."))
    call ptrdict_register_real_property(ptrdict_file, c_loc(cutoff_add), &
         CSTR("cutoff_add"), &
         CSTR("Additional shell added to the cutoff (skin depth/Verlet shell)."))

    call ptrdict_register_boolean_property(ptrdict_file, c_loc(in_sort), CSTR("sort"), &
         CSTR("Sort particles before binning."))

    call ptrdict_register_boolean_property(ptrdict_file, c_loc(force_binning), &
         CSTR("force_binning"), &
         CSTR("Force binning in every time step even when Verlet shells are enabled (cutoff_add > 0)."))   

    call ptrdict_register_boolean_property(ptrdict_file, c_loc(in_clear_velocities), &
         CSTR("clear_velocities"), &
         CSTR("Reset velocities to zero at the beginning of the simulation (overrides input file)."))
    call ptrdict_register_boolean_property(ptrdict_file, c_loc(in_clear_charges), &
         CSTR("clear_charges"), &
         CSTR("Reset charges to zero at the beginning of the simulation (overrides input file)."))

    !
    ! Initial transformations
    !

    translate_by  = 0.0_DP
    call ptrdict_register_point_property(ptrdict_file, c_loc(translate_by(1)), &
         CSTR("translate_by"), &
         CSTR("Initially translate the system by these values."))


    !
    ! Lees-Edwards stuff
    !

!    p%shear_dx(:)     = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
    p%shear_dv        = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
    shear_gamma       = (/ 0.0_DP, 0.0_DP, 0.0_DP /)
    shear_tau         = -1.0_DP
    shear_amplitude   = (/ 0.0_DP, 0.0_DP, 0.0_DP /)

    m = ptrdict_register_section(ptrdict_file, CSTR("LeesEdwards"), &
         CSTR("LeesEdwards shearing boundary conditions."))

!    call ptrdict_register_point_property(m, p%shear_dx(:), CSTR("dx"), &
!         CSTR("Constant (initial) displacement of the boxes."))
    call ptrdict_register_point_property(m, c_loc(p%shear_dv(1)), CSTR("dv"), &
         CSTR("Constant velocity."))
    call ptrdict_register_point_property(m, c_loc(shear_gamma(1)), &
         CSTR("gamma"), &
         CSTR("Shear rate."))
    call ptrdict_register_real_property(m, c_loc(shear_tau), CSTR("tau"), &
         CSTR("Oscillation period (no oscillation < 0.0)."))
    call ptrdict_register_point_property(m, c_loc(shear_amplitude(1)), CSTR("amplitude"), &
         CSTR("Amplitude of the oscillation (velocity)."))

    !
    ! Factories. Upon call to ptrdict_read the members of the *ints*, *pots*,
    ! *coul* and *cals* structures will be selectively allocated.
    !

    ! Implementation is in integrators_factory_c.c
    call integrators_factory_register(ptrdict_file, ints)

    ! Implementation is in potentials_factory_c.c
    call potentials_factory_register(ptrdict_file, pots)

    ! Implementation is in coulomb_factory_c.c
    call coulomb_factory_register(ptrdict_file, coul)

    ! Implementation is in callables_factory_c.c
    call callables_factory_register(ptrdict_file, cals)

#ifdef _MP
    call register(mod_communicator, ptrdict_file)
    m = ptrdict_register_section(ptrdict_file, CSTR("Parallel3D"), &
         CSTR("Domain decomposition module."))

    call ptrdict_register_intpoint_property(m, &
         c_loc(mod_communicator%decomposition(1)), CSTR("decomposition"), &
         CSTR("Number of domains in each direction, i.e. type of the decomposition."))
#endif

    ! Read the config file. After this call, *ints*, *pots*, *coul* and *cals*
    ! will contain instantiated objects of whatever is defined in md.dat.
    call ptrdict_read(ptrdict_file, CSTR("md.dat"))

    call ptrdict_write(ptrdict_file, CSTR("md.out"))

  endsubroutine read_ptrdict_file

endprogram main
