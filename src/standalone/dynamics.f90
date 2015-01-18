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

!>
!! Simulation dynamics
!!
!! Information on dynamics like velocities, energies, pressures etc.
!<
module dynamics
  use supplib
  use rng

  use particles

  implicit none

  private

  public :: dynamics_t
  type dynamics_t

     integer            :: it                !< Current iteration
     real(DP), pointer  :: ti                !< Current time
     real(DP)           :: dt                !< Time step
     real(DP)           :: maxtime          !< The maximum time for the simulation that was specified by the user via "md.dat"

     real(DP), pointer  :: v(:, :)           !< Velocities (stored in particles%data)
     real(DP), pointer  :: f(:, :)           !< Forces (stored in particles%data)

     real(DP)           :: epot              !< Potential energy (per atom or total?)
     real(DP)           :: ekin              !< Kinetic energy (per atom or total?)

     real(DP)           :: fmax

     real(DP)           :: wpot(3, 3)
     real(DP)           :: wkin(3, 3)

     real(DP)           :: pressure(3, 3)

     !
     ! Output counter
     !

     integer                :: nout

     !
     ! Particles object
     !

     type(particles_t), pointer  :: p

  endtype dynamics_t


  public :: init
  interface init
     module procedure dynamics_init
  endinterface

  public :: del
  interface del
     module procedure dynamics_del
  endinterface

  public :: update
  interface update
     module procedure dynamics_update
  endinterface

  public :: print_status
  interface print_status
     module procedure dynamics_print_status
  endinterface

  public :: give_velocities
  interface give_velocities
     module procedure dynamics_give_velocities
  end interface

  public :: remove_linear_momentum
  interface remove_linear_momentum
     module procedure dynamics_remove_linear_momentum
  end interface remove_linear_momentum

contains


  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine dynamics_init(this, p, dt, mymaxtime)
    implicit none

    type(dynamics_t), intent(inout)        :: this
    type(particles_t), target, intent(in)  :: p
    real(DP), intent(in), optional         :: dt
    real(DP), intent(in), optional         :: mymaxtime

    ! ---

    this%maxtime = 100.0_DP

    this%it    = 0

    this%nout  = 0

    this%p    => p

    this%v    => NULL()
    this%f    => NULL()
    this%ti   => NULL()

    this%dt    = 0.1_DP

    this%epot  = 0.0_DP
    this%ekin  = 0.0_DP

    this%wpot  = 0.0_DP
    this%wkin  = 0.0_DP

    if (present(dt)) then
       this%dt  = dt
    endif

    if (present(mymaxtime)) then
       this%maxtime = mymaxtime
    endif

    call add_real3( &
         this%p%data, &
         V_STR, &
         F_RESTART + F_TO_TRAJ + F_COMMUNICATE + F_COMM_FORCES, &
         "angstroms/femtosecond", &
         velocity_to_Afs )

    call add_real3( &
         this%p%data, &
         F_STR, &
         F_COMMUNICATE + F_COMM_FORCES + F_TO_TRAJ, &
         "eV/angstroms", &
         energy_to_eV/length_to_A )

    call add_real_attr(this%p%data, TI_ATTR_STR)

  endsubroutine dynamics_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine dynamics_del(this)
    implicit none

    type(dynamics_t), intent(inout)  :: this

    ! ---

  endsubroutine dynamics_del


  !>
  !! Advance time and compute dynamic properties
  !!
  !! Calculates the kinetic contribution to the pressure tensor
  !<
  subroutine dynamics_update(this, it, advance_time, mpi)
    implicit none

    type(dynamics_t), intent(inout)          :: this
    integer, intent(in), optional            :: it
    logical, intent(in), optional            :: advance_time
    type(MPI_context), intent(in), optional  :: mpi

    ! ---

    if (.not. associated(this%v))   call ptr_by_name(this%p%data, V_STR, this%v)
    if (.not. associated(this%f))   call ptr_by_name(this%p%data, F_STR, this%f)
    if (.not. associated(this%ti))  call attr_by_name(this%p%data, TI_ATTR_STR, this%ti)

    if (present(mpi)) then
       call sum_in_place(mpi, this%epot)
       call sum_in_place(mpi, this%wpot)
    endif

    if (present(advance_time)) then
       if (advance_time) then
          this%it  = this%it + 1
          this%ti  = this%ti + this%dt
       endif
    else
       this%it  = this%it + 1
       this%ti  = this%ti + this%dt
    endif

    if (present(it)) then
       this%it  = it
    endif

    call compute_kinetic_energy_and_virial( &
         this%p, this%v, this%f, &
         this%wpot, this%ekin, this%fmax, this%wkin, this%pressure, &
         mpi)

  endsubroutine dynamics_update


  !>
  !! Print a status log to screen
  !!
  !! Print a status log to screen
  !<
  subroutine dynamics_print_status(this)
    implicit none

    type(dynamics_t), intent(inout)  :: this

    ! ---

    real(DP)  :: T

    ! ---

    T  = this%ekin*2/(this%p%dof*K_to_energy)

#ifdef _MP
    if (mpi_id() == ROOT) then
#endif

       if (mod(this%nout, 10) == 0) then
          write (*, '(A10,A2,A10,A2,A10,A2,A12,A2,A12,A2,A12,A2,A12,A2,A10,A2,A12)')  &
               "it", " |", &
               "t[" // trim(time_str) // "]", " |", &
               "dt[" // trim(time_str) // "]", " |", &
               "ekin[" // trim(energy_str) // "]", " |", &
               "epot[" // trim(energy_str) // "]", " |", &
               "etot[" // trim(energy_str) // "]", " |", &
               "fmax[" // trim(force_str) // "]", " |", &
               "T[K]", " |", &
               "P[" // trim(pressure_str) // "]"
       endif

       write (*, '(I10,2X,F10.1,2X,F10.6,2X,ES12.5,2X,ES12.5,2X,ES12.5,2X,ES12.5,2X,F10.3,2X,ES12.3)')  &
            this%it, this%ti, this%dt, this%ekin, this%epot, this%ekin+this%epot, this%fmax, T, tr(3, this%pressure)/3

#ifdef _MP
    endif
#endif

    if (ilog /= -1 .and. mod(this%nout, 10) == 0) then
       write (ilog, '(A10,A12,A12,A22,A22,A22,A14,A12,A14)')  &
            "it", &
            "t[" // trim(time_str) // "]", &
            "dt[" // trim(time_str) // "]", &
            "ekin[" // trim(energy_str) // "]", &
            "epot[" // trim(energy_str) // "]", &
            "etot[" // trim(energy_str) // "]", &
            "fmax[" // trim(force_str) // "]", &
            "T[K]", &
            "P[" // trim(pressure_str) // "]"
    endif

    if (ilog /= -1) then
       write (ilog, '(I10,F12.1,F12.6,ES22.13,ES22.13,ES22.13,ES14.5,F12.3,ES14.3)')  &
            this%it, this%ti, this%dt, this%ekin, this%epot, this%ekin+this%epot, this%fmax, T, tr(3, this%pressure)/3
    endif

    this%nout = this%nout+1

  endsubroutine dynamics_print_status


  !>
  !! Give initial velocities
  !!
  !! Gives atoms initial velocities. Note that because of the equipartition
  !! theorem, half of the given kinetic energy will go to potential energy
  !! and hence the routine actually doubles the given temperature value.
  !!
  !! Only mode 1 is implemented at the moment.
  !!
  !! Mode 1: Give temperature according to the Maxwell-Boltzmann distribution
  !!
  !! \f[
  !!   \rho(v_{i}^{a}) = sqrt{\frac{m_a}{2 \pi k_B T}} \exp{-\frac{0.5 m_a (v^a_i)^2}{k_B T}},
  !! \f]
  !!
  !! where \$a\$ is the atom and \$i\$ the dimension.
  !!
  !<
  subroutine dynamics_give_velocities(this, p, T, mode, ierror)
    implicit none

    type(dynamics_t), intent(inout)  :: this    !< Dynamics object
    type(particles_t), intent(inout) :: p       !< Particles (for masses)
    real(DP), intent(in)             :: T       !< Temperature to set (Kelvin)
    integer, intent(in)              :: mode    !< Mode (1 - Gaussian velocities)
    integer, intent(inout), optional :: ierror  !< Error signals

    ! ---

    integer                          :: i, j    ! Loops
    real(DP)                         :: std     ! Standard deviation for velocities

    ! ---

    ! - Verify pointers
    if (.not. associated(this%v))   call ptr_by_name(p%data, V_STR, this%v)
    if (.not. associated(p%m))   call ptr_by_name(p%data, V_STR, p%m)

    ! - Check RNG
    if(.not. rng_initialized) then
       RAISE_ERROR("dynamics_give_velocities: RNG not initialized.", ierror)
    end if

    ! - Report

    call prlog("- dynamics_give_velocities -")

    ! - Set velocities

    if(mode==1) then
       ! Maxwell-Boltzmann distribution
       call prlog("     Maxwell-Boltzmann distribution, temperature (K):" // T)
       do i = 1, p%natloc
          std = sqrt(2*Boltzmann_K*T/p%m(i))  ! Units: sqrt(eV/au) = internal
          do j = 1, 3
             VEC(this%v, i, j) = std*rng_normal1()
          end do
       end do
    else
       RAISE_ERROR("dynamics_give_velocities: Unknown mode.", ierror)
    end if

    call prlog

  end subroutine dynamics_give_velocities

  !>
  !! Remove linear momentum
  !!
  !! Remove linear momentum
  !!
  !! Note: No MPI support.
  !<
  subroutine dynamics_remove_linear_momentum(this, p, ierror)
    implicit none

    type(dynamics_t), intent(inout)  :: this    !< Dynamics object
    type(particles_t), intent(in)    :: p       !< Particles (for masses)
    integer, intent(inout), optional :: ierror  !< Error signals

    ! ---

    integer                          :: i       ! Loops
    real(DP)                         :: vavg(3) ! Average velocity
    real(DP)                         :: msum    ! Sum of masses

    ! ---

    vavg  = 0.0_DP
    msum = 0.0_DP
    do i = 1, p%nat
       vavg  = vavg + p%m(i)*VEC3(this%v, i)
       msum = msum + p%m(i)
    enddo
    vavg  = vavg/msum

    do i = 1, p%nat
       VEC3(this%v, i)  = VEC3(this%v, i) - vavg
    enddo

  endsubroutine dynamics_remove_linear_momentum

  
endmodule dynamics
