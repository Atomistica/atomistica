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
!   classtype:gaussian_charges_t classname:GaussianCharges interface:coulomb
! @endmeta

!>
!! Gaussian broadening of charges.
!!
!! Gaussian broadening of charges.
!!
!! Assigns a shape to the charges on the atoms. Charge distributions is given
!!
!! \f[
!!   \rho(\vec{r}) = \sum_i q_i f_i(\vec{r}-\vec{r_i}) \right)
!! \f]
!!
!! where \f$f(\vec{r})\f$ the shape. This module implements Gaussian shapes.
!!
!! Note that this module does not compute the contribution of the singular,
!! long-ranged \f$1/r\f$ term.
!!
!! This module is required for both the TightBinding and VariableCharge
!! modules. For tight-binding calculations, \f$Z_i=0\f$ for all \f$i\f$.
!<

#include "macros.inc"
#include "filter.inc"

module gaussian_charges
  use supplib

  use particles
  use neighbors
  use filter

  implicit none

  private

  integer, parameter  :: GAUSSIAN_CHARGES_MAX_EL  = 3

  !
  ! The module for the computation of energies/potentials
  !

  type gaussian_charges_db_t

     integer   :: nel = -1
     integer   :: nU

     character :: el(2, GAUSSIAN_CHARGES_MAX_EL)   !< Atom type

     real(DP)  :: U(GAUSSIAN_CHARGES_MAX_EL)       !< Hubbard U

  endtype gaussian_charges_db_t

  public :: gaussian_charges_t
  type gaussian_charges_t

     !
     ! Element on which to apply the force
     !

     character(MAX_EL_STR) :: elements = "*"
     integer               :: els

     !
     ! real space sampling
     !

     real(DP) :: cutoff = 5.0_DP
     real(DP) :: cutoff_sq

     !
     ! Hubbard-U
     !

     real(DP), allocatable :: U(:)

     !
     ! Database
     !
     
     type(gaussian_charges_db_t) :: db

  endtype gaussian_charges_t


  public :: init
  interface init
     module procedure gaussian_charges_init
  endinterface

  public :: del
  interface del
     module procedure gaussian_charges_del
  endinterface

  public :: set_Hubbard_U
  interface set_Hubbard_U
     module procedure gaussian_charges_set_Hubbard_U
  endinterface

  public :: bind_to
  interface bind_to
     module procedure gaussian_charges_bind_to
  endinterface

  public :: potential
  interface potential
     module procedure gaussian_charges_potential
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure gaussian_charges_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure gaussian_charges_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a GaussianCharges object
  !<
  subroutine gaussian_charges_init(this, p, U, elements, cutoff, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(gaussian_charges_t),    intent(inout) :: this
    type(particles_t), optional, intent(in)    :: p
    real(DP),          optional, intent(in)    :: U(*)
    character(*),      optional, intent(in)    :: elements
    real(DP),          optional, intent(in)    :: cutoff
    integer,           optional, intent(inout) :: error

    ! ---

    call prlog("- gaussian_charges_init -")

    ASSIGN_PROPERTY(elements)
    ASSIGN_PROPERTY(cutoff)

    if (present(U)) then
       call set_Hubbard_U(this, p, U, error=error)
       PASS_ERROR(error)
    endif

    call prlog

  endsubroutine gaussian_charges_init


  !>
  !! Destructor
  !!
  !! Free all internal data buffers
  !<
  subroutine gaussian_charges_del(this)
    implicit none

    type(gaussian_charges_t), intent(inout) :: this

    ! ---

    if (allocated(this%U)) then
       deallocate(this%U)
    endif

  endsubroutine gaussian_charges_del


  !>
  !! Set the Hubbard U values and nuclear charges
  !!
  !! Set the Hubbard U values and nuclear charges Z. U and Z values are passed
  !! and stored per element, not per atom.
  !!
  !! Note that this needs to be called before *bind_to*!
  !<
  subroutine gaussian_charges_set_Hubbard_U(this, p, U, error)
    implicit none

    type(gaussian_charges_t), intent(inout) :: this  !> GaussianCharges object
    type(particles_t),        intent(in)    :: p
    real(DP),                 intent(in)    :: U(p%nel)
    integer,        optional, intent(out)   :: error

    ! ---

    INIT_ERROR(error)

    call prlog("- gaussian_charges_set_Hubbard_U -")
    call prlog("U = "//U)

    !this%db%nel = p%nel
    this%db%nU = p%nel

    !call resize(this%U, p%nel)
    this%db%U(1:p%nel) = U

    call prlog
    
  endsubroutine gaussian_charges_set_Hubbard_U


  !>
  !! Assign a Particles and a Neighbors object to this GaussianCharges object
  !!
  !! Assign a Particles and a Neighbors object to this GaussianCharges object. All subsequent operations
  !! will use this Particles object. Only a pointer to the object
  !! is copied, not the object itself.
  !!
  !! Note that this needs to be called *after* set_Hubbard_U!
  !<
  subroutine gaussian_charges_bind_to(this, p, nl, ierror)
    implicit none

    type(gaussian_charges_t), intent(inout) :: this
    type(particles_t),        intent(in)    :: p
    type(neighbors_t),        intent(inout) :: nl
    integer,        optional, intent(inout) :: ierror

    ! ---

    integer :: i, j, Z

    ! ---

    call prlog("- gaussian_charges_bind_to -")
    call prlog("elements = " // trim(this%elements))

    this%els  = filter_from_string(this%elements, p, ierror=ierror)
    PASS_ERROR(ierror)

    !
    ! Copy parameters to per element array
    !

    if (this%db%nel > 0 .and. this%db%nel /= this%db%nU ) then

       write (*, '(A,I2)')  "nel = ", this%db%nel
       write (*, '(A,I2)')  "nU  = ", this%db%nU

       RAISE_ERROR("The number of entries must be identical for all parameters.", ierror)
    endif

    !
    ! Convert units of Hubbard U's
    !

    call resize(this%U, p%nel)

    this%U = 0.0_DP
    if (this%db%nel > 0) then
       do j = 1, p%nel
          do i = 1, this%db%nel
             Z = atomic_number(a2s(this%db%el(:,i)))
             if (Z <= 0 .or. Z > MAX_Z) then
                RAISE_ERROR("Unknown element '" // trim(a2s(this%db%el(:,i))) // "'.", ierror)
             endif

             if (Z == p%el2Z(j)) then
                call prlog("     " // ElementName(Z) // " - " // j)
                this%U(j)  = this%db%U(i) / (Hartree*Bohr)
                call prlog("     - U = " // this%db%U(i) // " (" // this%U(j) // ")")
             endif
          enddo
       enddo
    else
       if (this%db%nU > 0) then
          this%U(1:this%db%nU) = this%db%U(1:this%db%nU) / (Hartree*Bohr)
       endif
    endif

    this%cutoff_sq = this%cutoff**2

    call request_interaction_range(nl, this%cutoff)
    call prlog("cutoff   = " // this%cutoff)

    call prlog

  endsubroutine gaussian_charges_bind_to


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge
  !! models)
  !!
  !! Difference between the Ewald potential and the real one due to
  !! the Gaussian charge distribution. You always need an additional
  !! DirectCoulomb, Ewald, etc.
  !<
  subroutine gaussian_charges_potential(this, p, nl, q, phi, ierror)
    implicit none

    type(gaussian_charges_t), intent(inout) :: this
    type(particles_t),        intent(in)    :: p
    type(neighbors_t),        intent(inout) :: nl
    real(DP),                 intent(in)    :: q(p%maxnatloc)
    real(DP),                 intent(inout) :: phi(p%maxnatloc)
    integer,        optional, intent(inout) :: ierror

    !---

    real(DP) :: U_i_sq, abs_rij, hlp, expi, expj, src, fac, fac2, efac
    real(DP) :: avg, fi1, fj1, fi2, fj2, U_i, U_j, q_i, q_j, Z_i, Z_j

    integer             :: i, j
    integer(NEIGHPTR_T) :: ni

    integer, parameter :: sqrt_pi_2 = sqrt(PI/2)

    !---

    call timer_start('gaussian_charges_potential')

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    !$omp  parallel default(none) &
    !$omp& shared(nl, this, phi, p, q) &
    !$omp& private(i, U_i_sq, j, ni, abs_rij, hlp)
    
    call tls_init(size(phi), sca=1)  ! is called tls_sca1 (=phi)
    
    !$omp do 
    do i = 1, p%natloc

       if (IS_EL(this%els, p, i)) then

          if (this%U(p%el(i)) > 0.0_DP) then

             !
             ! Atom i has a Gaussian charge cloud
             !

             U_i_sq = this%U(p%el(i))**2

             Gaussian_ni_loop: do ni = nl%seed(i), nl%last(i)
                j = GET_NEIGHBOR(nl, ni)

                if (i <= j .and. IS_EL(this%els, p, j)) then
                   abs_rij = GET_ABS_DRJ(p, nl, i, j, ni)
                   if (abs_rij < this%cutoff) then

                      if (this%U(p%el(j)) > 0.0_DP) then
                         hlp = -erfc(sqrt(PI/2*U_i_sq * this%U(p%el(j))**2/(U_i_sq + this%U(p%el(j))**2)) * abs_rij)/(abs_rij)

                         tls_sca1(i) = tls_sca1(i) + q(j)*hlp
                         if (i /= j)  tls_sca1(j) = tls_sca1(j) + q(i)*hlp
                      else
                         hlp = -erfc(sqrt_pi_2*this%U(p%el(i))*abs_rij)/(abs_rij)

                         tls_sca1(i) = tls_sca1(i) + q(j)*hlp
                         tls_sca1(j) = tls_sca1(j) + q(i)*hlp
                      endif

                   endif
                endif

             enddo Gaussian_ni_loop

          else

             !
             ! Atom i is a point charge
             !
             
             Gaussian_ni_loop2: do ni = nl%seed(i), nl%last(i)
                j = GET_NEIGHBOR(nl, ni)

                if (i < j .and. IS_EL(this%els, p, j)) then
                   if (this%U(p%el(j)) > 0.0_DP) then
                      abs_rij = GET_ABS_DRJ(p, nl, i, j, ni)

                      if (abs_rij < this%cutoff) then
                         hlp = -erfc(sqrt_pi_2*this%U(p%el(j))*abs_rij) &
                              /(abs_rij)

                         tls_sca1(i) = tls_sca1(i) + q(j)*hlp
                         tls_sca1(j) = tls_sca1(j) + q(i)*hlp
                      endif
                   endif
                endif

             enddo Gaussian_ni_loop2

          endif

          tls_sca1(i) = tls_sca1(i) + q(i)*this%U(p%el(i))! - this%jellium_potential
             
       endif

    enddo

    call tls_reduce(p%natloc, sca1=phi)
    !$omp end parallel

    call timer_stop('gaussian_charges_potential')

  endsubroutine gaussian_charges_potential


  !>
  !! Calculate the electrostatic potential and the electric field
  !!
  !! Difference between the Ewald (point charge) force and the force
  !! due to the Gaussian charge distribution
  !<
  subroutine gaussian_charges_energy_and_forces(this, p, nl, q, epot, f, wpot, &
       error)
    implicit none

    type(gaussian_charges_t), intent(inout) :: this
    type(particles_t),        intent(in)    :: p
    type(neighbors_t),        intent(inout) :: nl
    real(DP),                 intent(in)    :: q(p%maxnatloc)
    real(DP),                 intent(inout) :: epot
    real(DP),                 intent(inout) :: f(3, p%maxnatloc)
    real(DP),                 intent(inout) :: wpot(3, 3)
    integer,        optional, intent(inout) :: error

    !---

    real(DP) :: U_i_sq, q_i, q_j, rij(3), abs_rij, abs_rij_sq
    real(DP) :: c, df(3), hlp, sqrt_pi_2, src, fac, fac2, efac, expi, expj
    real(DP) :: avg, e, ffac, fi1, fj1, fi2, fj2, U_i, U_j, Z_i, Z_j

    integer             :: i, j
    integer(NEIGHPTR_T) :: ni

    !---

    call timer_start('gaussian_charges_energy_and_forces')

    call update(nl, p, error)
    PASS_ERROR(error)

    sqrt_pi_2  = sqrt(PI/2)
    
    do i = 1, p%natloc
       
       if (IS_EL(this%els, p, i)) then
             
          q_i = q(i)

          if (this%U(p%el(i)) > 0.0_DP) then

             !
             ! Atom i has a Gaussian charge cloud
             !

             U_i_sq = this%U(p%el(i))**2

             ni_loop: do ni = nl%seed(i), nl%last(i)
                j = GET_NEIGHBOR(nl, ni)
                   
                if (i <= j .and. IS_EL(this%els, p, j)) then
                   DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)

                   q_j = q(j)

                   if (abs_rij_sq < this%cutoff_sq) then
                      abs_rij  = sqrt(abs_rij_sq)

                      if (this%U(p%el(j)) > 0.0_DP) then
                         c    = PI/2*U_i_sq * this%U(p%el(j))**2 &
                              /(U_i_sq + this%U(p%el(j))**2)
                         hlp  = -erfc(sqrt(PI/2*U_i_sq * this%U(p%el(j))**2 &
                              /(U_i_sq + this%U(p%el(j))**2)) * abs_rij)&
                              /(abs_rij)
                      else
                         c    = PI/2*U_i_sq
                         hlp  = -erfc(sqrt_pi_2*this%U(p%el(i))*abs_rij) &
                              /(abs_rij)
                      endif

                      df = -rij/(abs_rij**3) *  &
                           ( erfc(sqrt(c) * abs_rij) &
                           + 2*sqrt(c/PI)*exp(-c*abs_rij_sq)*abs_rij &
                           )

                      if (i == j) then

                         epot     = epot + 0.5_DP*q_i*q_i*hlp
                         wpot     = wpot - outer_product(rij, 0.5_DP*q_i*q_j*df)
                      else

                         VEC3(f, i)  = VEC3(f, i) + q_i*q_j*df
                         VEC3(f, j)  = VEC3(f, j) - q_i*q_j*df

                         if (j > p%natloc) then
                            epot  = epot + 0.5_DP*q_i*q_j*hlp
                            wpot  = wpot - outer_product(rij, 0.5_DP*q_i*q_j*df)
                         else
                            epot  = epot + q_i*q_j*hlp
                            wpot  = wpot - outer_product(rij, q_i*q_j*df)
                         endif
                      endif
                   endif
                endif
             enddo ni_loop

          else

             !
             ! Atom i is a point charge
             !

             ni_loop2: do ni = nl%seed(i), nl%last(i)
                j = GET_NEIGHBOR(nl, ni)

                if (i < j .and. IS_EL(this%els, p, j)) then

                   if (this%U(p%el(j)) > 0.0_DP) then
                      DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)

                      q_j = q(j)

                      if (abs_rij_sq < this%cutoff_sq) then
                         abs_rij  = sqrt(abs_rij_sq)

                         c   = PI/2*this%U(p%el(j))**2
                         hlp = -erfc(sqrt_pi_2*this%U(p%el(j))*abs_rij) &
                              /(abs_rij)

                         df = -rij/(abs_rij**3) *  &
                              ( erfc(sqrt(c) * abs_rij) &
                              + 2*sqrt(c/PI)*exp(-c*abs_rij_sq)*abs_rij &
                              )

                         VEC3(f, i)  = VEC3(f, i) + q_i*q_j*df
                         VEC3(f, j)  = VEC3(f, j) - q_i*q_j*df

                         if (j > p%natloc) then
                            epot  = epot + 0.5_DP*q_i*q_j*hlp
                            wpot  = wpot - outer_product(rij, 0.5_DP*q_i*q_j*df)
                         else
                            epot  = epot + q_i*q_j*hlp
                            wpot  = wpot - outer_product(rij, q_i*q_j*df)
                         endif
                      endif

                   endif

                endif
             enddo ni_loop2

          endif

          epot   = epot + 0.5_DP*q(i)*q(i)*this%U(p%el(i)) 
                      
       endif

    enddo

    call timer_stop('gaussian_charges_energy_and_forces')

  endsubroutine gaussian_charges_energy_and_forces


  !!>
  !! Registry
  !!
  !! Expose parameters to the user
  !!<
  subroutine gaussian_charges_register(this, cfg, m)
    implicit none

    type(gaussian_charges_t), target, intent(inout) :: this
    type(c_ptr),                      intent(in)    :: cfg
    type(c_ptr),                      intent(out)   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("GaussianCharges"), &
         CSTR("This module assign a shape to each charge."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements."))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), &
         CSTR("cutoff"), &
         CSTR("Cutoff of the correction to the Coulomb potential."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, GAUSSIAN_CHARGES_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%U), GAUSSIAN_CHARGES_MAX_EL, c_loc(this%db%nU), &
         CSTR("U"), CSTR("Hubbard U."))

  endsubroutine gaussian_charges_register

endmodule gaussian_charges
