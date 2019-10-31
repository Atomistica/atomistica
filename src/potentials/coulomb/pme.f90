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
!   dependencies:pme_kernel.f90,fft3-public.f,fft_wrap.f
!   classtype:pme_t classname:PME interface:coulomb
! @endmeta

!>
!! Implementation of the smooth Particle-Mesh-Ewald method
!!
!! Implementation of the smooth Particle-Mesh-Ewald method
!!
!! See: Darden, York, Pedersen, J. Chem. Phys. 98, 10089 (1993)
!! Essmann, Perera, Berkowitz, Darden, Lee, Pedersen,
!! J. Chem. Phys. 103, 8577 (1995)
!!
!! Code adopted from ORAC under GPL license,
!! http://www.chim.unifi.it/orac/
!<

#include "macros.inc"

module pme
  use supplib

  use particles
  use neighbors

#ifdef _MP
  use mpi
#endif

  use pme_kernel

  implicit none

  private

  public :: pme_t
  type pme_t

     !
     ! Ewald parameters
     !

     real(DP)           :: alpha = 0.4_DP
     real(DP)           :: cutoff = 10.0_DP

     integer            :: grid(3)  = (/ 64, 64, 64 /)

     integer            :: order = 8

     !
     ! PME Grid
     !

     type(pme_grid_t)   :: pme_grid

     !
     ! Auxialliary stuff
     !

     real(DP)           :: sqrt_alpha
     real(DP)           :: sqrt_alpha_pi
     real(DP)           :: cutoff_sq

  endtype pme_t


  public :: init
  interface init
     module procedure pme_init
  endinterface

  public :: del
  interface del
     module procedure pme_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure pme_bind_to
  endinterface

  public :: potential
  interface potential
     module procedure pme_potential
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure pme_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure pme_register
  endinterface

  real(DP), parameter, private :: EPS = 1d-10

contains


  !>
  !! Constructor.
  !!
  !! Initialize an Pme object
  !<
  subroutine pme_init(this, cutoff, order, grid, error)
    implicit none

    type(pme_t),        intent(inout)  :: this
    real(DP), optional, intent(in)     :: cutoff
    integer,  optional, intent(in)     :: order
    integer,  optional, intent(in)     :: grid(3)
    integer,  optional, intent(out)    :: error

    ! ---

    INIT_ERROR(error)

    ASSIGN_PROPERTY(cutoff)
    ASSIGN_PROPERTY(order)
    ASSIGN_PROPERTY(grid)

  endsubroutine pme_init


  !>
  !! Destructor.
  !!
  !! Remove all dynamically allocated variables from memory.
  !<
  subroutine pme_del(this)
    implicit none

    type(pme_t), intent(inout)  :: this

    ! ---

    call del(this%pme_grid)

  endsubroutine pme_del


  !>
  !! Assign a Particles object to this Pme object
  !!
  !! Assign a Particles object to this Pme object. All subsequent operations
  !! will use this Particles object. Only a pointer to the object
  !! is copied, not the object itself.
  !<
  subroutine pme_bind_to(this, p, nl, ierror)
    implicit none

    type(pme_t), intent(inout)        :: this
    type(particles_t), intent(in)     :: p
    type(neighbors_t), intent(inout)  :: nl
    integer, optional, intent(out)    :: ierror

    !---

    ! checks
    if(.not. all(p%pbc /= 0)) then
       RAISE_ERROR("pme_bind_to: Pme summation only works for 3d periodics. You can use it for other systems as well by removing this error, but it will still be 3d.", ierror)
    end if

    write (ilog, '(A)')  "- pme_bind_to -"

    this%cutoff_sq  = this%cutoff**2
    this%alpha      = log(10.0d0)*12/this%cutoff_sq

    write (ilog, '(5X,A,F10.5)')  "cutoff   = ", this%cutoff
    write (ilog, '(5X,A,F10.5)')  "* alpha  = ", this%alpha
    write (ilog, '(5X,A,3I5)')    "grid     = ", this%grid

    this%sqrt_alpha     = sqrt(this%alpha)
    this%sqrt_alpha_pi  = sqrt(this%alpha/PI)

    call request_interaction_range(nl, this%cutoff)

    call init(this%pme_grid, this%grid, p%nat, this%order)

    write (ilog, *)

  endsubroutine pme_bind_to


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !!
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !<
  subroutine pme_potential(this, p, nl, q, phi, ierror)
    implicit none

    type(pme_t),        intent(inout) :: this
    type(particles_t),  intent(in)    :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(in)    :: q(p%maxnatloc)
    real(DP),           intent(inout) :: phi(p%maxnatloc)
    integer,  optional, intent(inout) :: ierror

    !---

    real(DP)  :: Abox(3, 3), Bbox(3, 3)

    !---

    real(DP)  :: dq_i

    real(DP)  :: hlp1, hlp(3)

    integer   :: i, ni, j
    real(DP)  :: rij(3), abs_rij, abs_rij_sq, virial(3, 3)

    real(DP)  :: epot_rec

    ! FIXME! This is on stack. Move somewhere else?
    real(DP)  :: x(p%nat), y(p%nat), z(p%nat)

    !---

    call timer_start('pme_potential')

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    !
    ! Direct sum
    !

    ! This loop ca. halved execution time with 2 threads, 500 atoms
    !$omp  parallel default(none) &
    !$omp& shared(p, this, nl, q, phi) &
    !$omp& private(dq_i, x) &
    !$omp& private(hlp, j, i, abs_rij) &
    !$omp& private(hlp1, rij, abs_rij_sq)

    call tls_init(p%maxnatloc, sca=1)  ! is called tls_sca1

    !$omp do
    do i = 1, p%nat

       dq_i = q(i)
       
       if (abs(dq_i) > EPS) then

          tls_sca1(i)        = tls_sca1(i) - 2*q(i)*this%sqrt_alpha_pi

          ni_loop: do ni = nl%seed(i), nl%last(i)
             j = nl%neighbors(ni)
             if (abs(q(j)) > EPS) then
                if (i < j) then
                   DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)
                   if (abs_rij_sq < this%cutoff_sq) then
                      abs_rij  = sqrt(abs_rij_sq)

                      hlp1 = erfc(this%sqrt_alpha*abs_rij)/abs_rij
                      tls_sca1(i) = tls_sca1(i) + q(j) * hlp1
                      tls_sca1(j) = tls_sca1(j) + dq_i * hlp1
                   endif
                else if (i == j) then
                   abs_rij  = GET_ABS_DRJ(p, nl, i, j, ni)
                   if (abs_rij < this%cutoff) then
                      hlp1 = erfc(this%sqrt_alpha*abs_rij)/abs_rij
                      tls_sca1(i) = tls_sca1(i) + q(j) * hlp1
                   endif
                endif
             endif

          enddo ni_loop

       endif
    enddo

    call tls_reduce(p%nat, sca1=phi)
    !$omp end parallel

    !
    ! Reciprocal sum (call stand-alone Darden routine)
    !

    x = POS(p, 1:p%nat, 1)
    y = POS(p, 1:p%nat, 2)
    z = POS(p, 1:p%nat, 3)

    epot_rec = 0.0_DP
    virial = 0.0_DP

    call get_true_cell(p, Abox, Bbox, error=ierror)
    PASS_ERROR(ierror)
  
    call potential_and_field( &
         this%pme_grid, x, y, z, q, Bbox, volume(p), this%sqrt_alpha, &
         epot_rec, virial, phi)

    call timer_stop('pme_potential')

  endsubroutine pme_potential


  !>
  !! Calculate the electrostatic potential and the electric field
  !!
  !! Return the electrostatic potential and the electric field (on each atom)
  !! alongside the total Coulomb energy. It uses the position from the
  !! associated Particles object and the charges from the respective charge
  !! array.
  !<
  subroutine pme_energy_and_forces(this, p, nl, q, epot, f, wpot, error)
    implicit none

    type(pme_t),       intent(inout) :: this
    type(particles_t), intent(in)    :: p
    type(neighbors_t), intent(inout) :: nl
    real(DP),          intent(in)    :: q(p%nat)
    real(DP),          intent(inout) :: epot
    real(DP),          intent(inout) :: f(3, p%nat)
    real(DP),          intent(inout) :: wpot(3, 3)
    integer, optional, intent(inout) :: error

    !---

    real(DP)  :: Abox(3, 3), Bbox(3, 3)

    real(DP)  :: dq_i

    real(DP)  :: hlp1, hlp(3)

    integer   :: i, ni, j
    real(DP)  :: rij(3), abs_rij, abs_rij_sq, wpot_dir(3, 3), wpot_rec(3, 3)

    real(DP)  :: epot_dir, epot_rec, epot_self

    ! FIXME! This is on stack. Move somewhere else?
    real(DP)  :: phi(p%nat), x(p%nat), y(p%nat), z(p%nat)
    real(DP)  :: Ex(p%nat), Ey(p%nat), Ez(p%nat)

    !---

    call timer_start('pme_energy_and_forces')

    call update(nl, p, error)
    PASS_ERROR(error)

    !
    ! Direct sum
    !

    epot_dir = 0.0_DP
    wpot_dir = 0.0_DP

    ! This loop ca. halved execution time with 2 threads, 500 atoms
    !$omp  parallel default(none) &
    !$omp& shared(p, this, nl, q, f) &
    !$omp& private(dq_i, x) &
    !$omp& private(hlp, j, i, abs_rij) &
    !$omp& private(hlp1, rij, abs_rij_sq) &
    !$omp& reduction(+:epot_dir) reduction(+:wpot_dir)

    call tls_init(p%maxnatloc, vec=1)

    !$omp do
    do i = 1, p%nat

       dq_i = q(i)
       
       if (abs(dq_i) > EPS) then

          ni_loop: do ni = nl%seed(i), nl%last(i)
             j = nl%neighbors(ni)
             if (abs(q(j)) > EPS) then
                if (i < j) then
                   DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)
                   if (abs_rij_sq < this%cutoff_sq) then
                      abs_rij = sqrt(abs_rij_sq)

                      hlp1     = erfc(this%sqrt_alpha*abs_rij)/abs_rij
                      epot_dir = epot_dir + dq_i*q(j) * hlp1

                      hlp = rij/(abs_rij**3) * &
                           ( erfc(this%sqrt_alpha*abs_rij) + &
                           2*this%sqrt_alpha_pi*abs_rij* &
                           exp(-this%alpha*abs_rij_sq) )

                      hlp = dq_i*q(j)*hlp

                      VEC3(tls_vec1, i) = VEC3(tls_vec1, i) + hlp
                      VEC3(tls_vec1, j) = VEC3(tls_vec1, j) - hlp

                      wpot_dir = wpot_dir - outer_product(rij, hlp)
                   endif
                else if (i == j) then
                   DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)
                   if (abs_rij_sq < this%cutoff_sq) then
                      abs_rij  = sqrt(abs_rij_sq)

                      hlp1 = erfc(this%sqrt_alpha*abs_rij)/abs_rij
                      epot_dir    = epot_dir + 0.5_DP * q(j)*q(j) * hlp1

                      hlp = rij/(abs_rij**3) * &
                           ( erfc(this%sqrt_alpha*abs_rij) + &
                           2*this%sqrt_alpha_pi*abs_rij* &
                           exp(-this%alpha*abs_rij_sq) )

                      hlp       = 0.5*dq_i*q(j)*hlp
                      wpot_dir  = wpot_dir - outer_product(rij, hlp)
                   endif
                endif
             endif

          enddo ni_loop

       endif
    enddo

    call tls_reduce(p%nat, vec1=f)
    !$omp end parallel

    !
    ! Reciprocal sum (call stand-alone Darden routine)
    !

    x = POS(p, 1:p%nat, 1)
    y = POS(p, 1:p%nat, 2)
    z = POS(p, 1:p%nat, 3)

    Ex = 0.0_DP
    Ey = 0.0_DP
    Ez = 0.0_DP

    epot_rec = 0.0_DP
    wpot_rec = 0.0_DP

    call get_true_cell(p, Abox, Bbox)

    call potential_and_field( &
         this%pme_grid, x, y, z, q, Bbox, volume(p), this%sqrt_alpha, &
         epot_rec, wpot_rec, phi, Ex, Ey, Ez)

    !
    ! Self energy
    !

    epot_self = 0
    
    do i = 1, p%nat
       epot_self = epot_self + q(i)**2
    enddo
    epot_self = epot_self*this%sqrt_alpha_pi

    VEC(f, 1:p%nat, 1) = VEC(f, 1:p%nat, 1) + q*Ex
    VEC(f, 1:p%nat, 2) = VEC(f, 1:p%nat, 2) + q*Ey
    VEC(f, 1:p%nat, 3) = VEC(f, 1:p%nat, 3) + q*Ez

    epot = epot + epot_dir + epot_rec - epot_self

    wpot = wpot + wpot_dir + wpot_rec

    call timer_stop('pme_energy_and_forces')

  endsubroutine pme_energy_and_forces


  !>
  !! Expose potential parameters
  !<
  subroutine pme_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(pme_t), target, intent(inout) :: this
    type(c_ptr),         intent(in)    :: cfg
    type(c_ptr),         intent(out)   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("PME"), &
         CSTR("Particle-mesh Ewald summation."))

    call ptrdict_register_real_property(m, c_loc(this%alpha), CSTR("alpha"), &
         CSTR("Gaussian broadening parameter."))
    call ptrdict_register_real_property(m, c_loc(this%cutoff), CSTR("cutoff"), &
         CSTR("Cut-off for the real-space sum."))
    call ptrdict_register_intpoint_property(m, c_loc(this%grid(1)), &
         CSTR("grid"), &
         CSTR("Dimension of the reciprocal space grid."))
    call ptrdict_register_integer_property(m, c_loc(this%order), CSTR("order"), &
         CSTR("Order of the polynomial interpolation."))

  endsubroutine pme_register

endmodule pme
