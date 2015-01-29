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
!   classtype:peters_t_t classname:PetersT interface:callables
! @endmeta

!>
!! Peters thermostat
!!
!! Peters (DPD) thermostat
!! See: E.A.J.F. Peters, Europhys. Lett. 66, 311 (2004)
!<

#include "macros.inc"

module peters_t
  use libAtoms_module

  use io
  use logging
  use rng
  use timer
  use tls

  use data
  use particles
  use neighbors
  use dynamics

  use interpolation_kernels

  implicit none

  private

  public :: peters_t_t
  type peters_t_t

     type(particles_t), pointer     :: p
     type(neighbors_t), pointer     :: nl

     type(interpolation_kernels_t)  :: interpolation_kernel

     !
     ! Parameters
     !

     integer   :: g           = 1

     logical(BOOL)   :: energy_out  = .false.

     real(DP)  :: T           = 300._DP
     real(DP)  :: dT          = 0._DP

     real(DP)  :: gamma       = 0.001_DP

     real(DP)  :: cutoff      = -1.0_DP

     !
     ! Other stuff
     !

     integer   :: un
     real(DP)  :: sum_de

     real(DP)  :: kT
     real(DP)  :: dkT

     !
     ! Interpolation kernel
     !

     !
     ! Velocities
     !

     real(DP), pointer  :: v(:, :)

  endtype peters_t_t


  public :: init
  interface init
     module procedure peters_t_init
  endinterface

  public :: del
  interface del
     module procedure peters_t_del
  endinterface

  public :: adjust_temperature
  interface adjust_temperature
     module procedure peters_t_adjust_temperature
  endinterface

  public :: invoke
  interface invoke
     module procedure peters_t_invoke
  endinterface

  public :: register
  interface register
    module procedure peters_t_register
  endinterface

!--- Internal

  interface set_particles
     module procedure peters_t_set_particles
  endinterface

  interface set_neighbors
     module procedure peters_t_set_neighbors
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a peters_t object
  !<
  subroutine peters_t_init(this, interpolation_kernel, group, T, dT, gamma, cutoff, energy_out)
    implicit none

    type(peters_t_t), intent(inout)                      :: this
    type(interpolation_kernels_t), intent(in), optional  :: interpolation_kernel
    integer, intent(in), optional                        :: group
    real(DP), intent(in), optional                       :: T
    real(DP), intent(in), optional                       :: dT
    real(DP), intent(in), optional                       :: gamma
    real(DP), intent(in), optional                       :: cutoff
    logical, intent(in), optional                        :: energy_out

    ! ---

    if (present(group)) then
       this%g  = group
    endif

    if (present(T)) then
       this%T  = T
    endif
    if (present(dT)) then
       this%dT  = dT
    endif

    if (present(gamma)) then
       this%gamma  = gamma
    endif

    if (present(cutoff)) then
       this%cutoff  = cutoff
    endif

    if (present(energy_out)) then
       this%energy_out  = energy_out
    endif

    if (present(interpolation_kernel)) then
       this%interpolation_kernel  = interpolation_kernel
       this%cutoff  = get_cutoff(this%interpolation_kernel)
    else
       call prlog("- peters_t_init -")
       call prlog("     Using default (square/linear) kernel.")

       allocate(this%interpolation_kernel%square)
       call init(this%interpolation_kernel%square, this%cutoff)

       call prlog
    endif

  endsubroutine peters_t_init


  !>
  !! Destructor
  !!
  !! Delete a peters_t object
  !<
  subroutine peters_t_del(this)
    implicit none

    type(peters_t_t), intent(inout)   :: this

    ! ---

    if (this%energy_out) then
       call fclose(this%un)
    endif

  endsubroutine peters_t_del


  !>
  !! Initialize a peters_t object
  !!
  !! Initialize a peters_t object
  !<
  subroutine peters_t_internal_init(this)
    implicit none

    type(peters_t_t), intent(inout)   :: this

    ! ---

    write (ilog, '(A)')  "- peters_t_internal_init -"

    call rng_init

    call ptr_by_name(this%p%data, V_STR, this%v)

    this%kT      = this%T * K_to_energy
    this%dkT     = this%dT * K_to_energy

    if (this%cutoff < 0.0_DP) then
       this%cutoff  = this%nl%interaction_range
    else
       call request_interaction_range(this%nl, this%cutoff)
    endif

    write (ilog, '(5X,A,ES20.10)')  "cutoff  = ", this%cutoff

    this%sum_de  = 0.0_DP

    if (this%energy_out) then
       this%un = fopen("peters_t.out", F_WRITE)
    endif

    write (ilog, *)

  endsubroutine peters_t_internal_init


  !>
  !! Set the associated particles object
  !!
  !! Set the associated particles object
  !<
  subroutine peters_t_set_particles(this, p)
    implicit none

    type(peters_t_t), intent(inout)  :: this
    type(particles_t), target        :: p

    ! ---

    this%p  => p

    if (associated(this%nl)) then
       call peters_t_internal_init(this)
    endif

  endsubroutine peters_t_set_particles


  !>
  !! Set the associated neighbors object
  !!
  !! Set the associated neighbors object
  !<
  subroutine peters_t_set_neighbors(this, nl)
    implicit none

    type(peters_t_t), intent(inout)  :: this
    type(neighbors_t), target        :: nl

    ! ---

    this%nl  => nl

    if (associated(this%p)) then
       call peters_t_internal_init(this)
    endif

  endsubroutine peters_t_set_neighbors


  !>
  !! Adjust the temperature
  !!
  !! Carries out a single Peters step. To be called after the second
  !! integrator step.
  !<
  subroutine peters_t_adjust_temperature(this, p, nl, dt, ti, ierror)
    implicit none

    type(peters_t_t), intent(inout)   :: this
    type(particles_t), target         :: p
    type(neighbors_t), target         :: nl
    real(DP), intent(in)              :: dt
    real(DP), intent(in), optional    :: ti
    integer, intent(inout), optional  :: ierror

    ! ---

    integer   :: i, ni, j

    real(DP)  :: weight, a, b, abs_dr, dr(3), r_mi
    real(DP)  :: vij(3), muij, r_muij, dmom(3)
    real(DP)  :: de, ekin(p%nat)

    ! ---

    call timer_start("peters_t_adjust_temperature")

    if (.not. associated(this%p, p))    call set_particles(this, p)
    if (.not. associated(this%nl, nl))  call set_neighbors(this, nl)

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    !$omp  parallel default(none) &
    !$omp& shared(dt, ekin, nl, p, this) &
    !$omp& private(r_mi, ni, i, j, abs_dr, dr, r_muij, muij, weight, a, b, vij, dmom)

    call tls_init(p%nat, vec=1)

    !$omp  do
    do i = 1, p%nat
       if (p%g(i) == this%g) then
          r_mi     = 1.0_DP/p%m(i)

          ekin(i)  = dot_product(VEC3(this%v, i), VEC3(this%v, i))

          do ni = nl%seed(i), nl%last(i)
             j = GET_NEIGHBOR(nl, ni)
             if (i < j .and. p%g(j) == this%g) then

!                abs_dr  = nl%abs_dr(ni)
                DIST_SQ(p, nl, i, ni, dr, abs_dr)

                if (abs_dr < this%cutoff**2) then
                   abs_dr             = sqrt(abs_dr)

!                   dr(:)  = VEC3(nl%dr, ni)/abs_dr

                   !
                   ! Thermalize pair
                   !

                   r_muij             = r_mi + 1.0_DP/p%m(j)
                   muij               = 1.0_DP/r_muij

                   ! XXX FIXME! This changes the normalization from the previous case.
                   call value_and_derivative(this%interpolation_kernel, abs_dr, a, weight)

                   weight             = dt*r_muij*this%gamma*weight

                   a                  = muij*( 1 - exp(-weight) )
                   b                  = sqrt( this%kT*muij * ( 1 - exp(-2*weight) ) )*rng_normal1()

                   vij                = VEC3(this%v, i) - VEC3(this%v, j) + VEC(nl%dc, ni, 3)*p%shear_dv

                   dmom               = ( -a*dot_product(vij, dr) + b ) * dr

                   VEC3(tls_vec1, i)  = VEC3(tls_vec1, i) + dmom/p%m(i)
                   VEC3(tls_vec1, j)  = VEC3(tls_vec1, j) + (- dmom/p%m(j))

                endif

             endif
          enddo
       endif
!       i = p%next_particle(i)
!    enddo
    enddo

!    !$omp critical
!    VEL3(p, 1:p%nat)  = VEL3(p, 1:p%nat) + VEC3(tls_vec1, 1:p%nat)
!    !$omp end critical

    call tls_reduce(p%nat, vec1=this%v)

    !$omp end parallel

    if (this%energy_out) then
       de  = 0.0_DP

       !$omp  parallel do default(none) &
       !$omp& shared(ekin, p, this) &
       !$omp& reduction(+:de)
       do i = 1, p%nat
          if (p%g(i) == this%g) then
             de  = de + p%m(i)*( dot_product(VEC3(this%v, i), VEC3(this%v, i)) - ekin(i) )/2
          endif
       enddo

       this%sum_de  = this%sum_de + de

       if (present(ti)) then
          write (this%un, '(F12.1,4ES20.10)')  ti, de/dt, this%sum_de
       endif
    endif

    this%T   = this%T  + dt*this%dT
    this%kT  = this%kT + dt*this%dkT

    call timer_stop("peters_t_adjust_temperature")

  endsubroutine peters_t_adjust_temperature


  !>
  !! Adjust the temperature
  !!
  !! Adjust the temperature
  !<
  subroutine peters_t_invoke(this, dyn, nl, ierror)
    implicit none

    type(peters_t_t), intent(inout)   :: this
    type(dynamics_t), intent(inout)   :: dyn
    type(neighbors_t), intent(in)     :: nl
    integer, intent(inout), optional  :: ierror

    ! ---

    call adjust_temperature(this, dyn%p, nl, dyn%dt, dyn%ti, ierror)
    PASS_ERROR(ierror)

  endsubroutine peters_t_invoke


  subroutine peters_t_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(peters_t_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)             :: cfg
    type(c_ptr), intent(out)            :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("PetersT"), &
         CSTR("Peters thermostat (E.A.J.F. Peters, Europhys. Lett. 66, 311 (2004))."))

    call ptrdict_register_integer_property(m, c_loc(this%g), CSTR("group"), &
         CSTR("Group of particles to thermalize."))

    call ptrdict_register_real_property(m, c_loc(this%T), CSTR("T"), &
         CSTR("Target temperature."))
    call ptrdict_register_real_property(m, c_loc(this%dT), CSTR("dT"), &
         CSTR("Linear temperature ramp."))

    call ptrdict_register_real_property(m, c_loc(this%gamma), CSTR("gamma"), &
         CSTR("Temperature coupling dissipation constant."))

    call ptrdict_register_boolean_property(m, c_loc(this%energy_out), CSTR("energy_out"), &
         CSTR("Write energy balance to 'peter_t.out'."))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), CSTR("cutoff"), &
         CSTR("Cut-off distance. A linear weighting function is used"))

  endsubroutine peters_t_register

endmodule peters_t
