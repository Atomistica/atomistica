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
!   shared
!   classtype:berendsen_t_t classname:BerendsenT interface:callables
! @endmeta

!>
!! Berendsen temperature control
!!
!! Berendsen temperature control according to
!! H. J. C. Berendsen et al., J. Chem. Phys. 81, 3684 (1984).
!!
!! The temperature is controlled by multiplying the velocities
!! with the factor
!! \f[
!!   s = \sqrt{1 + \frac{\Delta t}{\tau} \left( \frac{T_{0}}{T} - 1 \right)}
!! \f]
!!
!! (Except in place of the temperatures the kinetic energies are used. Please
!!  see the source code.)
!<

#include "macros.inc"

module berendsen_t
  use supplib

  use particles
  use neighbors
  use dynamics

  implicit none

  private

  integer, parameter  :: n_dims       = 4
  integer, parameter  :: len_dim_str  = 15
  integer, parameter  :: ALL_DIMS     = 0

  ! This is need for xlf
  character(len_dim_str), parameter  :: STR_all          = CSTR("all")
  character(len_dim_str), parameter  :: STR_x            = CSTR("x")
  character(len_dim_str), parameter  :: STR_y            = CSTR("y")
  character(len_dim_str), parameter  :: STR_z            = CSTR("z")
  character(len_dim_str), parameter  :: dim_strs(n_dims) = &
       (/ STR_all, STR_x, STR_y, STR_z /)

  public :: berendsen_t_t
  type berendsen_t_t

     !>
     !! Rescale once at the beginning
     !<
     logical(BOOL)   :: rescale_once = .false.
    
     !>
     !! Dimensions to control, 1 - x, 2 - y, 3 - z, 0 - all
     !<
     integer   :: d  = ALL_DIMS

     !>
     !! Desired temperature
     !<
     real(DP)  :: T  = 300.0_DP

     !>
     !! Final temperature - rescaling only!
     !<
     real(DP)  :: T0  = -1.0_DP

     !>
     !! Temperature change per time unit for linear quenches
     !<
     real(DP)  :: dT = 0.0_DP

     !>
     !! Time constant
     !<
     real(DP)  :: tau  = 500.0_DP 

     !>
     !! Iteration
     !<
     integer   :: it = 0

  endtype berendsen_t_t


  public :: init
  interface init
     module procedure berendsen_t_init
  endinterface

  public :: set
  interface set
     module procedure berendsen_t_init
  endinterface

  public :: adjust_temperature
  interface adjust_temperature
     module procedure berendsen_t_adjust_temperature
  endinterface

  public :: invoke
  interface invoke
     module procedure berendsen_t_invoke
  endinterface

  public :: register
  interface register
    module procedure berendsen_t_register
  endinterface

contains


  !>
  !! Initialize the temperature control specifying parameters
  !!
  !! Initialize the temperature control specifying parameters.
  !!
  !! Unless given, the following default parameters are used:
  !!
  !! T      300 K (scale towards room temperature)
  !!
  !! d      0     (scale in all dimensions)
  !!
  !! tau    500   (in internal time units)
  !<
  subroutine berendsen_t_init(this, d, T, T0, dT, tau)
    implicit none

    type(berendsen_t_t), intent(inout)  :: this
    integer,  optional,  intent(in)     :: d
    real(DP), optional,  intent(in)     :: T
    real(DP), optional,  intent(in)     :: T0
    real(DP), optional,  intent(in)     :: dT
    real(DP), optional,  intent(in)     :: tau

    ! ---

    if (present(d)) then
       this%d = d
    endif
    if (present(T)) then
       this%T = T
    endif
    if (present(T0)) then
       this%T0 = T0
    endif
    if (present(dT)) then
       this%dT = dT
    endif
    if (present(tau)) then
       this%tau = tau
    endif

    call prlog("- berendsen_t_init -")
    if (this%rescale_once) then
       call prlog("     Rescaling velocities at the beginning of the simulation")
    endif
    call prlog("     Using Berendsen temperature control with parameters")
    call prlog("     T   = " // this%T)
    call prlog("     T0  = " // this%T0)
    call prlog("     dT  = " // this%dT)
    call prlog("     tau = " // this%tau)
    call prlog("     d   = " // this%d)
    call prlog

  endsubroutine berendsen_t_init


  !>
  !! Adjust temperature
  !!
  !! Adjust temperature. See the detailed description of the
  !! berendsen_t module for details.
  !<
  subroutine berendsen_t_adjust_temperature(this, p, v, ekin, dt, ierror)
    implicit none

    type(berendsen_t_t), intent(inout)   :: this                  !< The T control object
    type(particles_t), intent(inout)     :: p                     !< Particles object
    real(DP), intent(inout)              :: v(3, p%maxnatloc)     !< Velocities
    real(DP), intent(in)                 :: ekin                  !< Current kinetic energy (of the whole system)
    real(DP), intent(in)                 :: dt                    !< Time step
    integer, intent(inout), optional     :: ierror                !< Error passing

    ! ---

    real(DP)  :: tau
    real(DP)  :: s              ! velocity scaling factor
    real(DP)  :: desired_ekin   ! desired Ekin per dof calculated from T

    ! ---

    this%it = this%it + 1

    ! If this%T0 > 0 we do velocity rescaling, but relax the temperature
    ! exponentially
    if (this%tau > 0.0_DP .and. this%T0 > 0.0_DP) then
       this%T = this%T + dt*(this%T0-this%T)/this%tau
    endif

    ! initialize
    desired_ekin = this%T * K_to_energy / 2
    if(this%d < 0 .or. this%d > 3) then
       RAISE_ERROR("Berendsen T control: d parameter not between 0 and 3.", ierror)
    endif

    call timer_start("berendsen_t_adjust_temperature")

    if (ekin > 0.0_DP .and. ( this%it == 1 .or. .not. this%rescale_once )) then

       tau = this%tau

       ! Rescale at the beginning
       if (this%rescale_once .and. this%it == 1) then
          tau = -1.0_DP
       endif

       if (tau > 0.0_DP .and. this%T0 < 0.0_DP) then
          s = sqrt(1+dt*(desired_ekin*p%dof/ekin-1)/tau)
       else
          s = sqrt(desired_ekin*p%dof/ekin)
       endif

       if (this%d == ALL_DIMS) then
          v(:, 1:p%natloc) = v(:, 1:p%natloc) * s
       else
          v(this%d, 1:p%natloc) = v(this%d, 1:p%natloc) * s
       endif

    endif

    ! Linear temperature change
    this%T = this%T + this%dT*dt

    call timer_stop("berendsen_t_adjust_temperature")

  endsubroutine berendsen_t_adjust_temperature


  !>
  !! Adjust the temperature
  !!
  !! Adjust the temperature
  !<
  subroutine berendsen_t_invoke(this, dyn, nl, ierror)
    implicit none

    type(berendsen_t_t), intent(inout)  :: this
    type(dynamics_t), intent(inout)     :: dyn
    type(neighbors_t), intent(in)       :: nl
    integer, intent(inout), optional    :: ierror

    ! ---

    call adjust_temperature(this, dyn%p, dyn%v, dyn%ekin, dyn%dt, ierror)
    PASS_ERROR(ierror)

  endsubroutine berendsen_t_invoke


  subroutine berendsen_t_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(berendsen_t_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)             :: cfg
    type(c_ptr), intent(out)            :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("BerendsenT"), &
         CSTR("Berendsen thermostat (H.J.C. Berendsen, J.P.M. Postma, W.F. van Gunsteren, A. DiNola, J.R. Haak, J. Chem. Phys. 81, 3684 (1984)."))

    call ptrdict_register_boolean_property(m, c_loc(this%rescale_once), CSTR("rescale_once"), &
         CSTR("Rescale only once at the beginning of the simulation."))

    call ptrdict_register_real_property(m, c_loc(this%T), CSTR("T"), &
         CSTR("Initial temperature."))

    call ptrdict_register_real_property(m, c_loc(this%T0), CSTR("T0"), &
         CSTR("Target temperature (switches on velocity rescaling!)."))

    call ptrdict_register_real_property(m, c_loc(this%dT), CSTR("dT"), &
         CSTR("Linear temperature change (in temp/time)."))

    call ptrdict_register_real_property(m, c_loc(this%tau), CSTR("tau"), &
         CSTR("Temperature coupling time constant. If <= 0, instant rescaling is used."))

    call ptrdict_register_enum_property(m, c_loc(this%d), &
         n_dims, len_dim_str, dim_strs(:), &
         CSTR("d"), &
         CSTR("Dimension to thermalize: 'x', 'y', 'z' or 'all'"))

  endsubroutine berendsen_t_register

endmodule berendsen_t
