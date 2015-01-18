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
!   dependencies:verlet_support.f90
!   classtype:andersen_p_t classname:AndersenP interface:integrators
! @endmeta

!>
!! Andersen pressure control
!!
!! Andersen pressure control
!! See: H.C. Andersen, J. Chem. Phys. 72, 2384 (1980)
!<

#include "macros.inc"
#include "filter.inc"

module andersen_p
  use supplib
  use rng

  use particles
  use dynamics

  use filter
  use verlet_support

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  character(*), parameter :: ANDERSEN_P_ENERGY_STR   = "AndersenP::energy"
  character(*), parameter :: ANDERSEN_P_ETA_STR = "AndersenP::eta"

  public :: andersen_p_t
  type andersen_p_t

     !>
     !! Comma separated list of elements this integrator acts on
     !<
     character(MAX_EL_STR)  :: elements = "*"

     !>
     !! Internal filter for the list of elements this integrator acts on
     !<
     integer                :: els      = 0

     !>
     !! Target pressure
     !<
     real(DP)               :: P(3)     = 0.0_DP

     !>
     !! Fictious barostat mass
     !<
     real(DP)               :: W        = 1.0_DP

     !>
     !! Barostat energy
     !<
     real(DP), pointer      :: energy   => NULL()

     !>
     !! Internal state variables
     !<
     real(DP), pointer      :: eta(:)   => NULL()

  endtype andersen_p_t


  public :: init
  interface init
     module procedure andersen_p_init
  endinterface

  public :: step1_with_dyn
  interface step1_with_dyn
     module procedure andersen_p_step1
  endinterface

  public :: step2_with_dyn
  interface step2_with_dyn
     module procedure andersen_p_step2
  endinterface

  public :: register
  interface register
    module procedure andersen_p_register
  endinterface

contains

  !>
  !! Constructor
  !<
  subroutine andersen_p_init(this, p)
    implicit none

    type(andersen_p_t), intent(inout) :: this
    type(particles_t),  intent(inout) :: p

    ! ---

    call add_real_attr(p%data, ANDERSEN_P_ENERGY_STR, F_TO_ENER)
    call add_real3_attr(p%data, ANDERSEN_P_ETA_STR, F_RESTART)

  endsubroutine andersen_p_init


  !>
  !! Position update and velocity estimation
  !<
  subroutine andersen_p_step1(this, dyn, max_dt, max_dr, max_dr_sq)
    implicit none

    type(andersen_p_t), intent(inout) :: this
    type(dynamics_t),   intent(inout) :: dyn
    real(DP), optional, intent(in)    :: max_dt
    real(DP), optional, intent(in)    :: max_dr
    real(DP), optional, intent(inout) :: max_dr_sq

    ! ---

    type(particles_t), pointer :: p

    integer  :: i

    real(DP) :: dt, dt2, l_max_dr_sq, V, L(3), pr(3), L_dt2(3), L_dt(3)
    real(DP) :: vfac(3), rfac(3)

    ! ---

    call timer_start("andersen_p_step1")

    p => dyn%p
    call require_orthorhombic_cell(p)

    if (this%els == 0) then
       this%els = filter_from_string(this%elements, p)
    endif

    if (.not. associated(this%energy)) then
       call attr_by_name(p%data, ANDERSEN_P_ENERGY_STR, this%energy)
    endif
    if (.not. associated(this%eta)) then
       call attr_by_name(p%data, ANDERSEN_P_ETA_STR, this%eta)
    endif

    !
    ! Copy some variables
    !

    dt = dyn%dt
    dt2 = dt/2

    !
    ! Adaptive time stepping
    ! FIXME! Does this actually work?
    !

    call timestep(p, dyn%v, dyn%f, dyn%dt, max_dt, max_dr)
        
    !
    ! Integrate
    !

    V        = volume(p)
    L        = (/ p%Abox(1,1), p%Abox(2,2), p%Abox(3,3) /)
    pr       = (/ dyn%pressure(1,1), dyn%pressure(2,2), dyn%pressure(3,3) /)

    this%eta = this%eta + dt2 * V/L * ( pr - this%P )

    call verlet_v(this%els, p, dyn%v, dyn%f, dt)

    L_dt2    = L + dt*this%eta/(2*this%W)
    L_dt     = L_dt2 + dt*this%eta/(2*this%W)

    call verlet_r(this%els, p, dyn%v, dyn%f, dt, l_max_dr_sq, &
         fac = (L/L_dt2)**2)

    call set_cell(p, L_dt, scale_atoms=.false.)
    vfac     = L/L_dt
    rfac     = L_dt/L
    !$omp  parallel do default(none) &
    !$omp& shared(dyn, p) firstprivate(vfac, rfac)
    do i = 1, p%natloc
       VEC3(dyn%v, i) = vfac*VEC3(dyn%v, i)
       PNC3(p, i)     = rfac*PNC3(p, i)
    enddo

    l_max_dr_sq = max(l_max_dr_sq, maxval((L_dt-L)**2))

    !
    ! Maximum particle displacement
    !

    p%accum_max_dr  = p%accum_max_dr + sqrt(l_max_dr_sq)

    if (present(max_dr_sq)) then
       max_dr_sq  = max(max_dr_sq, l_max_dr_sq)
    endif

    call I_changed_positions(p)

    call timer_stop("andersen_p_step1")

  endsubroutine andersen_p_step1


  !>
  !! Velocity correction
  !<
  subroutine andersen_p_step2(this, dyn)
    implicit none

    type(andersen_p_t), intent(in)    :: this
    type(dynamics_t),   intent(inout) :: dyn

    ! ---

    type(particles_t), pointer :: p

    real(DP) :: pressure(3, 3), V, L(3), pr(3)

    ! ---

    call timer_start("andersen_p_step2")

    p => dyn%p

    !
    ! Communicate forces back to this processor if required
    !

#ifdef _MP
    if (mod_communicator%communicate_forces) then
       DEBUG_WRITE("- communicate_forces -")
       call communicate_forces(mod_communicator, p)
    endif
#endif

    !
    ! Update barostat momentum
    !

    call compute_kinetic_energy_and_virial(p, dyn%v, dyn%f, dyn%wpot, &
        pressure=pressure)

    V  = volume(p)
    L  = (/ p%Abox(1,1), p%Abox(2,2), p%Abox(3,3) /)
    pr = (/ pressure(1,1), pressure(2,2), pressure(3,3) /)

    !
    ! Integrate
    !

    call verlet_v(this%els, p, dyn%v, dyn%f, dyn%dt)

    this%eta = this%eta + dyn%dt/2 * V/L * ( pr - this%P )

    this%energy = 0.5_DP*this%W*dot_product(this%eta, this%eta) + dot_product(pr, L)

    call timer_stop("andersen_p_step2")
    
  endsubroutine andersen_p_step2


  subroutine andersen_p_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(andersen_p_t), target, intent(inout) :: this
    type(c_ptr),                intent(in)    :: cfg
    type(c_ptr),                intent(out)   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("AndersenP"), &
         CSTR("The Andersen barostat (NPH ensemble)."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements for which to enable this integrator."))

    call ptrdict_register_point_property(m, c_loc(this%P(1)), CSTR("P"), &
         CSTR("Target pressure"))
    call ptrdict_register_real_property(m, c_loc(this%W), CSTR("W"), &
         CSTR("Fictious barostat mass"))

  endsubroutine andersen_p_register

endmodule andersen_p
