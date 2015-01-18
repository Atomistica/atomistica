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
!   classtype:direct_coulomb_t classname:DirectCoulomb interface:coulomb
! @endmeta

!>
!! Evaluates the Coulomb interaction by direct summation over all pairs.
!!
!! Evaluate the Coulomb interaction by direct summation over all pairs.
!! This scales N^2 and only works for non-periodic systems.
!<

#include "macros.inc"

module direct_coulomb
  use supplib

  use particles
  use neighbors

  implicit none

  private

  public :: direct_coulomb_t
  type direct_coulomb_t

     real(DP)   :: epsilon_r  = 1.0_DP    !> Relative dielectric constant

  endtype direct_coulomb_t


  public :: init
  interface init
     module procedure direct_coulomb_init
  endinterface

  public :: del
  interface del
     module procedure direct_coulomb_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure direct_coulomb_bind_to
  endinterface

  public :: potential
  interface potential
     module procedure direct_coulomb_potential
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure direct_coulomb_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure direct_coulomb_register
  endinterface

contains

  !>
  !! Constructor.
  !!
  !! Initialize a DirectCoulomb object
  !<
  subroutine direct_coulomb_init(this, epsilon_r, error)
    implicit none

    type(direct_coulomb_t), intent(inout)  :: this       !> DirectCoulomb object
    real(DP), intent(in), optional         :: epsilon_r  !> Relative constant of permittivity
    integer, intent(out), optional         :: error

    !---

    INIT_ERROR(error)

#ifdef _MP
    RAISE_ERROR("The DirectCoulomb module does not work with MPI.", error)
#endif

    if (present(epsilon_r)) then
       this%epsilon_r  = epsilon_r
    endif

  endsubroutine direct_coulomb_init


  !>
  !! Destructor.
  !!
  !! Uninitialize a DirectCoulomb object
  !<
  subroutine direct_coulomb_del(this)
    implicit none

    type(direct_coulomb_t), intent(inout)  :: this !> DirectCoulomb object

    !---
    
  endsubroutine direct_coulomb_del


  !>
  !! Assign a Neighbors object to this DirectCoulomb object
  !!
  !! Does nothing.
  !<
  subroutine direct_coulomb_bind_to(this, p, nl, ierror)
    implicit none

    type(direct_coulomb_t), intent(inout)  :: this    !< DirectCoulomb object
    type(particles_t), intent(in)          :: p       !< Particles object
    type(neighbors_t), intent(inout)       :: nl      !< Neighbors object
    integer, intent(inout), optional       :: ierror  !< Error passing

    ! ---

  endsubroutine direct_coulomb_bind_to


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !!
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !<
  subroutine direct_coulomb_potential(this, p, nl, q, phi, ierror)
    implicit none

    type(direct_coulomb_t), intent(inout)  :: this
    type(particles_t), intent(in)          :: p
    type(neighbors_t), intent(in)          :: nl
    real(DP), intent(in)                   :: q(p%maxnatloc)
    real(DP), intent(inout)                :: phi(p%maxnatloc)
    integer, intent(inout), optional       :: ierror

    ! ---

    integer   :: i, j
    real(DP)  :: abs_dr, dr(3)

    ! ---

    do i = 1, p%natloc-1
       do j = i+1, p%natloc
          dr      = POS3(p, i)-POS3(p, j)
          abs_dr  = 1.0_DP/(this%epsilon_r*sqrt(dot_product(dr, dr)))
          phi(i)  = phi(i) + q(j)*abs_dr
          phi(j)  = phi(j) + q(i)*abs_dr
       enddo
    enddo

  endsubroutine direct_coulomb_potential


  !!>
  !! Calculate the electrostatic potential and the electric field
  !!
  !! Return the electrostatic potential and the electric field (on each atom) alongside
  !! the total Coulomb energy. It uses the position from the associated Particles object
  !! and the charges from the respective charge array.
  !!<
  subroutine direct_coulomb_energy_and_forces(this, p, nl, q, epot, f, wpot, error)
    implicit none

    type(direct_coulomb_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    real(DP),               intent(in)    :: q(p%maxnatloc)
    real(DP),               intent(inout) :: epot
    real(DP),               intent(inout) :: f(3, p%maxnatloc)
    real(DP),               intent(inout) :: wpot(3, 3)
    integer,      optional, intent(inout) :: error

    ! --

    integer  :: i, j
    real(DP) :: abs_dr, dr(3), df(3), fac, q_i
    real(DP) :: energy, virial(3, 3)

    ! ---

    energy  = 0.0_DP
    virial  = 0.0_DP
    fac     = 1.0_DP/this%epsilon_r

    !$omp  parallel default(none) &
    !$omp& firstprivate(fac) &
    !$omp& shared(f, p, q) &
    !$omp& private(abs_dr, df, dr, i, j, q_i) &
    !$omp& reduction(+:energy) reduction(+:virial)

    call tls_init(p%natloc, vec=1)

    !$omp do
    do i = 1, p%natloc-1
       q_i  = q(i)

       do j = i+1, p%natloc
          dr                = POS3(p, i)-POS3(p, j)
          abs_dr            = sqrt(dot_product(dr, dr))
          df                = fac*q_i*q(j)*dr/(abs_dr**3)
          VEC3(tls_vec1, i) = VEC3(tls_vec1, i) + df
          VEC3(tls_vec1, j) = VEC3(tls_vec1, j) - df
          energy            = energy + fac*q_i*q(j)/abs_dr
          virial            = virial - outer_product(dr, df)
       enddo
    enddo

    call tls_reduce(p%natloc, vec1=f)

    !$omp end parallel

    epot  = epot + energy
    wpot  = wpot + virial

  endsubroutine direct_coulomb_energy_and_forces


  !!>
  !! Registry
  !!
  !! Expose parameters to the user
  !!<
  subroutine direct_coulomb_register(this, cfg, m)
    implicit none

    type(direct_coulomb_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("DirectCoulomb"), &
         CSTR("Evaluation of the Coulomb potential by direct summation."))

    call ptrdict_register_real_property(m, c_loc(this%epsilon_r), &
         CSTR("epsilon_r"), &
         CSTR("Relative constant of permittivity."))

  endsubroutine direct_coulomb_register

endmodule direct_coulomb
