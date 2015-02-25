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
!   classtype:cutoff_coulomb_t classname:CutoffCoulomb interface:coulomb
! @endmeta

!>
!! Coulomb evaluation by cutting-off the interaction at a certain distance.
!!
!! Coulomb evaluation by cutting-off the interaction at a certain distance.
!!
!! Direct evaluation of the Coulomb potential with a cut-off radius. Beyond
!! that cut-off, the interaction energies are set to zero instantly. No smooth
!! cut-off function is used here.
!!
!! USE WITH CARE!
!<

#include "macros.inc"

module cutoff_coulomb
  use supplib

  use particles
  use neighbors

  implicit none

  private

  public :: cutoff_coulomb_t
  type cutoff_coulomb_t

     real(DP)   :: epsilon_r  = 1.0_DP    !> Relative dielectric constant
     real(DP)   :: cutoff     = 10.0_DP   !> Cut-off distance

  endtype cutoff_coulomb_t


  public :: init
  interface init
     module procedure cutoff_coulomb_init
  endinterface

  public :: del
  interface del
     module procedure cutoff_coulomb_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure cutoff_coulomb_bind_to
  endinterface

  public :: potential
  interface potential
     module procedure cutoff_coulomb_potential
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure cutoff_coulomb_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure cutoff_coulomb_register
  endinterface

contains

  !>
  !! Constructor.
  !!
  !! Initialize a CutoffCoulomb object
  !<
  subroutine cutoff_coulomb_init(this, cutoff, epsilon_r, error)
    implicit none

    type(cutoff_coulomb_t), intent(inout)  :: this
    real(DP), intent(in), optional         :: cutoff
    real(DP), intent(in), optional         :: epsilon_r
    integer, intent(out), optional         :: error

    !---

    INIT_ERROR(error)

    call prscrlog("Warning: CutoffCoulomb does not interpolate to zero and should not really be used unless you add interpolation.")

#ifdef _MP
    RAISE_ERROR("The DirectCoulomb module does not (yet) work with MPI.", error)
#endif

    if (present(cutoff)) then
       this%cutoff  = cutoff
    endif
       
    if (present(epsilon_r)) then
       this%epsilon_r  = epsilon_r
    endif

  endsubroutine cutoff_coulomb_init


  !>
  !! Destructor.
  !!
  !! Uninitialize a CutoffCoulomb object
  !<
  subroutine cutoff_coulomb_del(this)
    implicit none

    type(cutoff_coulomb_t), intent(inout)  :: this !> CutoffCoulomb object

    !---

  endsubroutine cutoff_coulomb_del


  !>
  !! Assign a Neighbors object to this CutoffCoulomb object
  !!
  !! Assign a Neighbors object to this CutoffCoulomb object. All subsequent operations
  !! will use this Neighbors object. Only a pointer to the object
  !! is copied, not the object itself.
  !<
  subroutine cutoff_coulomb_bind_to(this, p, nl, ierror)
    implicit none

    type(cutoff_coulomb_t), intent(inout)  :: this    !< CutoffCoulomb object
    type(particles_t), intent(in)          :: p       !< Particles object
    type(neighbors_t), intent(inout)       :: nl      !< Neighbors object
    integer, intent(inout), optional       :: ierror  !< Error passing

    ! ---

    call request_interaction_range(nl, this%cutoff)

  endsubroutine cutoff_coulomb_bind_to


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !!
  !! Calculate the electrostatic potential of every atom (for variable charge models)
  !<
  subroutine cutoff_coulomb_potential(this, p, nl, q, phi, ierror)
    implicit none

    type(cutoff_coulomb_t), intent(inout)  :: this
    type(particles_t), intent(in)          :: p
    type(neighbors_t), intent(inout)       :: nl
    real(DP), intent(in)                   :: q(p%maxnatloc)
    real(DP), intent(inout)                :: phi(p%maxnatloc)
    integer, intent(inout), optional       :: ierror

    ! ---

    integer             :: i, j
    integer(NEIGHPTR_T) :: ni

    real(DP) :: abs_dr, dr(3), cutoff_sq

    ! ---

    call timer_start('cutoff_coulomb_potential')

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    cutoff_sq  = this%cutoff**2

    !$omp  parallel default(none) &
    !$omp& shared(this, nl, q, p, cutoff_sq, phi) &
    !$omp& private(i, ni, j, dr, abs_dr)

    call tls_init(p%maxnatloc, sca=1)  ! is called tls_sca1

    !$omp do
    do i = 1, p%natloc
       do ni = nl%seed(i), nl%last(i)
          j = GET_NEIGHBOR(nl, ni)
          if(i < j) then
             DISTJ_SQ(p, nl, i, ni, j, dr, abs_dr)
             if (abs_dr < cutoff_sq) then
                abs_dr  = 1.0_DP/(this%epsilon_r*sqrt(abs_dr))
                tls_sca1(i)  = tls_sca1(i) + q(j)*abs_dr
                tls_sca1(j)  = tls_sca1(j) + q(i)*abs_dr
             endif
          end if
       enddo
    enddo

    call tls_reduce(p%nat, sca1=phi)
    !$omp end parallel

    call timer_stop('cutoff_coulomb_potential')

  endsubroutine cutoff_coulomb_potential


  !!>
  !! Calculate the electrostatic potential and the electric field
  !!
  !! Return the electrostatic potential and the electric field (on each atom) alongside
  !! the total Coulomb energy. It uses the position from the associated Particles object
  !! and the charges from the respective charge array.
  !!<
  subroutine cutoff_coulomb_energy_and_forces(this, p, nl, q, epot, f, wpot, error)
    implicit none

    type(cutoff_coulomb_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    real(DP),               intent(in)    :: q(p%nat)
    real(DP),               intent(inout) :: epot
    real(DP),               intent(inout) :: f(3, p%nat)
    real(DP),               intent(inout) :: wpot(3, 3)
    integer,      optional, intent(out)   :: error

    ! --

    integer             :: i, j
    integer(NEIGHPTR_T) :: ni

    real(DP) :: abs_dr, dr(3), df(3), fac, q_i, cutoff_sq

    ! ---

    INIT_ERROR(error)

    call timer_start('cutoff_coulomb_energy_and_forces')

    call update(nl, p, error)
    PASS_ERROR(error)

    cutoff_sq  = this%cutoff**2
    fac        = 1.0_DP/this%epsilon_r
    do i = 1, p%natloc
       q_i  = q(i)

       do ni = nl%seed(i), nl%last(i)
          DISTJ_SQ(p, nl, i, ni, j, dr, abs_dr)

          if (i <= j .and. abs_dr < cutoff_sq) then
             abs_dr        = sqrt(abs_dr)
             df            = q_i*q(j)*fac*dr/(abs_dr**3)
             if (i == j) then
                epot        = epot + 0.5_DP*fac*q_i*q(j)/abs_dr
                wpot        = wpot - outer_product(dr, 0.5_DP*df)
             else
                VEC3(f, i)  = VEC3(f, i) + df
                VEC3(f, j)  = VEC3(f, j) - df

                epot        = epot + fac*q_i*q(j)/abs_dr
                wpot        = wpot - outer_product(dr, df)
             endif
          endif
       enddo
    enddo

    call timer_stop('cutoff_coulomb_energy_and_forces')

  endsubroutine cutoff_coulomb_energy_and_forces


  !!>
  !! Registry
  !!
  !! Expose parameters to the user
  !!<
  subroutine cutoff_coulomb_register(this, cfg, m)
    implicit none

    type(cutoff_coulomb_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("CutoffCoulomb"), &
         CSTR("Evaluation of the Coulomb potential by direct summation with a cutoff."))

    call ptrdict_register_real_property(m, c_loc(this%epsilon_r), &
         CSTR("epsilon_r"), &
         CSTR("Relative constant of permittivity."))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), CSTR("cutoff"), &
         CSTR("Cutoff radius."))

  endsubroutine cutoff_coulomb_register

endmodule cutoff_coulomb
