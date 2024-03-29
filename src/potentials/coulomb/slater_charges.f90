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
!   dependencies:damp_short_gamma.f90,coulomb_short_gamma.f90
!   classtype:slater_charges_t classname:SlaterCharges interface:coulomb
! @endmeta

!>
!! Slater-type broadening of charges.
!!
!! Slater-type broadening of charges.
!!
!! Assigns a shape to the charges on the atoms. Charge distributions is given
!!
!! \f[
!!   \rho(\vec{r}) = \sum_i \left( Z_i \delta^3(\vec{r}-\vec{r_i}) + (q_i-Z_i) f_i(\vec{r}-\vec{r_i}) \right)
!! \f]
!!
!! where \f$Z_i\f$ is the effective nuclear charge and \f$f(\vec{r})\f$ the
!! shape.
!! This module implements Slater-type (exponential) shapes.
!!
!! Note that this module does not compute the contribution of the singular,
!! long-ranged \f$1/r\f$ term.
!!
!! This module is required for both the TightBinding and VariableCharge
!! modules. For tight-binding calculations, \f$Z_i=0\f$ for all \f$i\f$.
!<

#include "macros.inc"
#include "filter.inc"

module slater_charges
  use supplib

  use particles
  use neighbors
  use filter

  !
  ! DFTB3
  !
  use coulomb_short_gamma
  use damp_short_gamma
 
  implicit none

  private

  integer, parameter  :: SLATER_CHARGES_MAX_EL = 16

  !
  ! The module for the computation of energies/potentials
  !

  type slater_charges_db_t

     integer   :: nel = -1
     integer   :: nU = -1, nZ = -1

     character :: el(2, SLATER_CHARGES_MAX_EL)   !< Atom type

     real(DP)  :: U(SLATER_CHARGES_MAX_EL)       !< Hubbard U
     real(DP)  :: Z(SLATER_CHARGES_MAX_EL)       !< Nuclear charge

     ! for DFTB3
     real(DP)  :: dU(0:116) = 0.0_DP             !< Hubbard derivative with respect to atomic charge
 
  endtype slater_charges_db_t

  public :: slater_charges_t
  type slater_charges_t

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
     ! DFTB3
     !

     logical   :: dftb3 = .false.

     logical   :: damp_gamma = .false.
     real(DP)  :: zeta = 4.00000_DP

     !
     ! Hubbard-U
     !

     real(DP), allocatable :: U(:)

     !
     ! Effective nuclear charge
     !

     real(DP), allocatable :: Z(:)

     !
     ! for DFTB3
     ! Hubbard derivative with respect to atomic charge
     !

     real(DP), allocatable :: dU(:)

     !
     ! Database
     !
     
     type(slater_charges_db_t) :: db

  endtype slater_charges_t


  public :: init
  interface init
     module procedure slater_charges_init
  endinterface

  public :: del
  interface del
     module procedure slater_charges_del
  endinterface

  public :: set_Hubbard_U
  interface set_Hubbard_U
     module procedure slater_charges_set_Hubbard_U
  endinterface

  public :: bind_to
  interface bind_to
     module procedure slater_charges_bind_to
  endinterface

  public :: potential
  interface potential
     module procedure slater_charges_potential
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure slater_charges_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure slater_charges_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a SlaterCharges object
  !<
  subroutine slater_charges_init(this, p, U, elements, cutoff, error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(slater_charges_t),      intent(inout) :: this
    type(particles_t), optional, intent(in)    :: p
    real(DP),          optional, intent(in)    :: U(*)
    character(*),      optional, intent(in)    :: elements
    real(DP),          optional, intent(in)    :: cutoff
    integer,           optional, intent(inout) :: error

    ! ---

    call prlog("- slater_charges_init -")

    ASSIGN_PROPERTY(elements)
    ASSIGN_PROPERTY(cutoff)

    if (present(U)) then
       call set_Hubbard_U(this, p, U, error=error)
       PASS_ERROR(error)
    endif

    call prlog

  endsubroutine slater_charges_init


  !>
  !! Destructor
  !!
  !! Free all internal data buffers
  !<
  subroutine slater_charges_del(this)
    implicit none

    type(slater_charges_t), intent(inout) :: this

    ! ---

    if (allocated(this%U)) then
       deallocate(this%U)
    endif
    if (allocated(this%Z)) then
       deallocate(this%Z)
    endif

  endsubroutine slater_charges_del


  !>
  !! Set the Hubbard U values and nuclear charges
  !!
  !! Set the Hubbard U values and nuclear charges Z. U and Z values are passed
  !! and stored per element, not per atom.
  !!
  !! Note that this needs to be called before *bind_to*!
  !<
  subroutine slater_charges_set_Hubbard_U(this, p, U, Z, error)
    implicit none

    type(slater_charges_t), intent(inout) :: this  !> SlaterCharges object
    type(particles_t),      intent(in)    :: p
    real(DP),               intent(in)    :: U(p%nel)
    real(DP),     optional, intent(in)    :: Z(p%nel)
    integer,      optional, intent(out)   :: error

    ! ---

    INIT_ERROR(error)

    call prlog("- slater_charges_set_Hubbard_U -")
    call prlog("U = "//U)
    if (present(Z)) then
       call prlog("Z = "//Z)
    endif

    if (p%nel > SLATER_CHARGES_MAX_EL) then
       RAISE_ERROR("Number of elements > SLATER_CHARGES_MAX_EL", error)
    endif

    !this%db%nel = p%nel
    this%db%nU = p%nel
    this%db%nZ = p%nel

    !call resize(this%U, p%nel)
    this%db%U(1:p%nel) = U
    
    !call resize(this%Z, p%nel)
    if (present(Z)) then
       this%db%Z(1:p%nel) = Z
    else
       this%db%Z = 0.0_DP
    endif

    call prlog

  endsubroutine slater_charges_set_Hubbard_U


  !>
  !! Assign a Particles and a Neighbors object to this SlaterCharges object
  !!
  !! Assign a Particles and a Neighbors object to this SlaterCharges object. All subsequent operations
  !! will use this Particles object. Only a pointer to the object
  !! is copied, not the object itself.
  !!
  !! Note that this needs to be called *after* set_Hubbard_U!
  !<
  subroutine slater_charges_bind_to(this, p, nl, ierror)
    implicit none

    type(slater_charges_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    integer,      optional, intent(inout) :: ierror

    ! ---

    integer :: i, j, Z

    ! ---

    call prlog("- slater_charges_bind_to -")
    call prlog("elements   = " // trim(this%elements))
    call prlog("dftb3      = " // this%dftb3)
    call prlog("damp_gamma = " // this%damp_gamma)
    call prlog("zeta       = " // this%zeta)

    this%els  = filter_from_string(this%elements, p, ierror=ierror)
    PASS_ERROR(ierror)

    !
    ! Copy parameters to per element array
    !

    if (this%db%nel > 0 .and. &
         ( this%db%nel /= this%db%nU .or. &
         this%db%nel /= this%db%nZ )) then

       write (*, '(A,I2)')  "nel = ", this%db%nel
       write (*, '(A,I2)')  "nU  = ", this%db%nU
       write (*, '(A,I2)')  "nZ  = ", this%db%nZ

       RAISE_ERROR("The number of entries must be identical for all parameters.", ierror)
    endif

    !
    ! Convert units of Hubbard U's
    !

    call resize(this%U, p%nel)
    call resize(this%dU, p%nel)
    call resize(this%Z, p%nel)

    if (this%dftb3) then
       call prlog("element Hubbard-U (in Hartree) (in eV) derivative (in Hartree) (in eV) nuclear charge")
       call prlog("=======           ============ =======            ============ ======= ==============")
    else
       call prlog("element Hubbard-U (in Hartree) (in eV) nuclear charge")
       call prlog("=======           ============ ======= ==============")
    endif

    this%U = 0.0_DP
    this%Z = 0.0_DP
    if (this%db%nel > 0) then
       do j = 1, p%nel
          do i = 1, this%db%nel
             Z = atomic_number(a2s(this%db%el(:,i)))
             if (Z <= 0 .or. Z > MAX_Z) then
                RAISE_ERROR("Unknown element '" // trim(a2s(this%db%el(:,i))) // "'.", ierror)
             endif

             if (Z == p%el2Z(j)) then
                this%U(j)  = this%db%U(i) / (Hartree*Bohr)
                this%dU(j) = this%db%dU(Z) / (Hartree*Bohr)
                this%Z(j)  = this%db%Z(i)
                if (this%dftb3) then
                   write (ilog, '(5X,A7,11X,F12.3,F8.3,12X,F12.3,F8.3,F15.3)') ElementName(Z), this%db%U(i), this%U(j), this%db%dU(Z), this%dU(j), this%db%Z(i)
                else
                   write (ilog, '(5X,A7,11X,F12.3,F8.3,F15.3)') ElementName(Z), this%db%U(i), this%U(j), this%db%Z(i)
                endif
             endif
          enddo
       enddo
    else
       if (this%db%nU /= p%nel .or. this%db%nZ /= p%nel) then
          RAISE_ERROR("this%db%nU /= p%nel .or. this%db%nZ /= p%nel", ierror)
       endif
       do i = 1, p%nel
          Z = p%el2Z(i)
          this%U(i)  = this%db%U(i) / (Hartree*Bohr)
          this%dU(i) = this%db%dU(Z) / (Hartree*Bohr)
          this%Z(i)  = this%db%Z(i)
          if (this%dftb3) then
             write (ilog, '(5X,A7,11X,F12.3,F8.3,12X,F12.3,F8.3,F15.3)') ElementName(Z), this%db%U(i), this%U(i), this%db%dU(Z), this%dU(i), this%db%Z(i)
          else
             write (ilog, '(5X,A7,11X,F12.3,F8.3,F15.3)') ElementName(Z), this%db%U(i), this%U(i), this%db%Z(i)
          endif
       enddo
    endif

    ! U_i is converted such that the charge is f_i(r) ~ exp(-U_i r), i.e.
    ! U_i = tau_i from Elstner's paper or U_i = 2 zeta_i from
    ! Streitz-Mintmire's paper.

    this%U = 16*this%U/5

    this%cutoff_sq = this%cutoff**2

    call request_interaction_range(nl, this%cutoff)
    call prlog("cutoff   = " // this%cutoff)

    call prlog

  endsubroutine slater_charges_bind_to


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge
  !! models)
  !!
  !! Difference between the Ewald potential and the real one due to
  !! the Gaussian charge distribution. You always need an additional
  !! DirectCoulomb, Ewald, etc.
  !<
  subroutine slater_charges_potential(this, p, nl, q, phi, ierror)
    implicit none

    type(slater_charges_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    real(DP),               intent(in)    :: q(p%maxnatloc)
    real(DP),               intent(inout) :: phi(p%maxnatloc)
    integer,      optional, intent(inout) :: ierror

    !---

    real(DP) :: abs_rij, hlp, expi, expj, src, fac, fac2, efac
    real(DP) :: avg, fi1, fj1, fi2, fj2, U_i, U_j, q_i, q_j, Z_i, Z_j

    !
    ! DFTB3 
    !

    real(DP)     :: dU_i, dU_j, dq_i, dq_j
    real(DP)     :: cgma_ac, cgma_bc, cgma_ca, cgma_cb

    integer             :: i, j, atomic_number_i, atomic_number_j
    integer(NEIGHPTR_T) :: ni

    !---

    call timer_start('slater_charges_potential')

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    ! U_i has been converted such that the charge is f_i(r) ~ exp(-U_i r),
    ! i.e. U_i = tau_i from Elstner's paper or U_i = 2 zeta_i from
    ! Streitz-Mintmire's paper.

    !$omp  parallel default(none) &
    !$omp& shared(nl, this, phi, p, q) &
    !$omp& private(i, j, q_i, q_j, U_i, U_j, Z_i, Z_j) &
    !$omp& private(ni, abs_rij, hlp, src, fac, avg) &
    !$omp& private(fac2, efac, fi1, fi2, fj1, fj2, expi, expj) &
    !$omp& private(atomic_number_i, atomic_number_j, dU_i, dq_i, dU_j, dq_j) &
    !$omp& private(cgma_ac, cgma_bc, cgma_ca, cgma_cb)  
     
    call tls_init(size(phi), sca=1)  ! is called tls_sca1 (=phi)

    !$omp do 
    do i = 1, p%natloc

       if (IS_EL(this%els, p, i)) then

          q_i = q(i)
          U_i = this%U(p%el(i))
          Z_i = this%Z(p%el(i))

          atomic_number_i = p%Z(i)
          dU_i = this%dU(p%el(i))
          dq_i = q_i - Z_i 

          !
          ! Atom i has a Gaussian charge cloud
          !

          Slater_ni_loop: do ni = nl%seed(i), nl%last(i)
             j = GET_NEIGHBOR(nl, ni)
            
             if (i <= j .and. IS_EL(this%els, p, j)) then
                abs_rij = GET_ABS_DRJ(p, nl, i, j, ni)
                if (abs_rij < this%cutoff) then

                   q_j = q(j)
                   U_j = this%U(p%el(j))
                   Z_j = this%Z(p%el(j))

                   atomic_number_j = p%Z(j)
                   dU_j = this%dU(p%el(j))
                   dq_j = q_j - Z_j 

                   !
                   ! Nuclear repulsion integrals
                   !

                   hlp = -(0.5_DP*U_i+1.0_DP/abs_rij)*exp(-U_i*abs_rij)
                   tls_sca1(i) = tls_sca1(i) + Z_j*hlp
                   hlp = -(0.5_DP*U_j+1.0_DP/abs_rij)*exp(-U_j*abs_rij)
                   tls_sca1(j) = tls_sca1(j) + Z_i*hlp

                   if (abs(U_i - U_j) < 1d-6) then

                      !
                      ! Coulomb integrals
                      !

                      src = 1.0_DP/(U_i+U_j)
                      fac = U_i*U_j*src
                      avg = 1.6_DP*(fac+fac*fac*src)
                      fac = avg*abs_rij
                      fac2 = fac*fac
                      efac = exp(-fac)/(48*abs_rij)

                      hlp = -(48 + 33*fac + fac2*(9+fac))*efac

                      !
                      ! XH correction for DFTB2
                      !

                      if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
                         hlp = hlp*hij(abs_rij, U_i, U_j, this%zeta)
                      endif

                      tls_sca1(i) = tls_sca1(i) + dq_j*hlp
                      tls_sca1(j) = tls_sca1(j) + dq_i*hlp

                      !
                      ! DFTB3 Hamiltonian
                      !

                      if (this%dftb3) then

                         if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then

                            cgma_ac = capital_short_gamma(abs_rij, dU_i, U_i, U_j, this%zeta)
                            cgma_ca = capital_short_gamma(abs_rij, dU_j, U_j, U_i, this%zeta)    
                 
                         else
 
                            cgma_ac = capital_short_gamma(abs_rij, dU_i, U_i, U_j)
                            cgma_ca = capital_short_gamma(abs_rij, dU_j, U_j, U_i)    
                     
                         endif

                         cgma_bc = cgma_ca
                         cgma_cb = cgma_ac 
 
                         tls_sca1(i) = tls_sca1(i) &
                                     - dq_i*dq_j*cgma_ac*2.0_DP/3.0_DP &
                                     - dq_j**2*cgma_ca/3.0_DP
                         tls_sca1(j) = tls_sca1(j) &
                                     - dq_j*dq_i*cgma_bc*2.0_DP/3.0_DP &
                                     - dq_i**2*cgma_cb/3.0_DP

                      endif

                   else

                      !
                      ! Coulomb integrals
                      !

                      fi1 = 1.0_DP/(2*(U_i**2-U_j**2)**2)
                      fj1 = -U_i**4*U_j*fi1
                      fi1 = -U_j**4*U_i*fi1
                         
                      fi2 = 1.0_DP/((U_i**2-U_j**2)**3)
                      fj2 = -(U_i**6-3*U_i**4*U_j**2)*fi2
                      fi2 =  (U_j**6-3*U_j**4*U_i**2)*fi2

                      expi = exp(-U_i*abs_rij)
                      expj = exp(-U_j*abs_rij)

                      hlp = expi*(fi1+fi2/abs_rij) + expj*(fj1+fj2/abs_rij)

                      !
                      ! XH correction for DFTB2 
                      !

                      if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
                         hlp = hlp*hij(abs_rij, U_i, U_j, this%zeta)
                      endif

                      tls_sca1(i) = tls_sca1(i) + (q_j-Z_j)*hlp
                      tls_sca1(j) = tls_sca1(j) + (q_i-Z_i)*hlp

                      !
                      ! DFTB3 Hamiltonian
                      !

                      if (this%dftb3) then

                         if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then

                            cgma_ac = capital_short_gamma(abs_rij, dU_i, U_i, U_j, this%zeta)
                            cgma_ca = capital_short_gamma(abs_rij, dU_j, U_j, U_i, this%zeta)    
  
                         else

                            cgma_ac = capital_short_gamma(abs_rij, dU_i, U_i, U_j)
                            cgma_ca = capital_short_gamma(abs_rij, dU_j, U_j, U_i)    
                     
                         endif

                         cgma_bc = cgma_ca
                         cgma_cb = cgma_ac 
                      
                         tls_sca1(i) = tls_sca1(i) &
                                     - dq_i*dq_j*cgma_ac*2.0_DP/3.0_DP &
                                     - dq_j**2*cgma_ca/3.0_DP
                         tls_sca1(j) = tls_sca1(j) &
                                     - dq_j*dq_i*cgma_bc*2.0_DP/3.0_DP &
                                     - dq_i**2*cgma_cb/3.0_DP

                      endif

                   endif
                      
                endif
             endif

          enddo Slater_ni_loop

          tls_sca1(i) = tls_sca1(i) + 5*q_i*U_i/16
 
          !
          ! DFTB3 on-site correction
          ! 

          if (this%dftb3) tls_sca1(i) = tls_sca1(i) - 0.50_DP*q_i**2*dU_i
                  
       endif

    enddo

    call tls_reduce(p%natloc, sca1=phi)
    !$omp end parallel

    call timer_stop('slater_charges_potential')

  endsubroutine slater_charges_potential


  !>
  !! Calculate the electrostatic potential and the electric field
  !!
  !! Difference between the Ewald (point charge) force and the force
  !! due to the Gaussian charge distribution
  !<
  subroutine slater_charges_energy_and_forces(this, p, nl, q, epot, f, wpot, error)
    implicit none

    type(slater_charges_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    real(DP),               intent(in)    :: q(p%maxnatloc)
    real(DP),               intent(inout) :: epot
    real(DP),               intent(inout) :: f(3, p%maxnatloc)
    real(DP),               intent(inout) :: wpot(3, 3)
    integer,      optional, intent(inout) :: error

    !---

    real(DP) :: q_i, q_j, rij(3), abs_rij, abs_rij_sq
    real(DP) :: df(3), hlp, src, fac, fac2, efac, expi, expj
    real(DP) :: avg, e, ffac, fi1, fj1, fi2, fj2, U_i, U_j, Z_i, Z_j

    !
    ! DFTB3
    !

    real(DP)     :: dU_i, dU_j, dq_i, dq_j
    real(DP)     :: edftb3, fdftb3

    integer             :: i, j, atomic_number_i, atomic_number_j
    integer(NEIGHPTR_T) :: ni

    !---

    call timer_start('slater_charges_energy_and_forces')

    call update(nl, p, error)
    PASS_ERROR(error)

    ! U_i has been converted such that the charge is f_i(r) ~ exp(-U_i r),
    ! i.e. U_i = tau_i from Elstner's paper or U_i = 2 zeta_i from
    ! Streitz-Mintmire's paper.

    do i = 1, p%natloc

       if (IS_EL(this%els, p, i)) then

          q_i = q(i)
          U_i = this%U(p%el(i))
          Z_i = this%Z(p%el(i))

          atomic_number_i = p%Z(i)
          dU_i = this%dU(p%el(i))
          dq_i = q_i - Z_i

          !
          ! Atom i has a Gaussian charge cloud
          !

          Slater_ni_loop: do ni = nl%seed(i), nl%last(i)
             j   = GET_NEIGHBOR(nl, ni)
                
             if (i <= j .and. IS_EL(this%els, p, j)) then
                DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)
                   
                if (abs_rij_sq < this%cutoff_sq) then
                   abs_rij  = sqrt(abs_rij_sq)

                   q_j = q(j)
                   U_j = this%U(p%el(j))
                   Z_j = this%Z(p%el(j))

                   atomic_number_j = p%Z(j)
                   dU_j = this%dU(p%el(j))
                   dq_j = q_j - Z_j

                   !
                   ! Nuclear repulsion integrals
                   !

                   hlp = -(0.5_DP*U_i+1.0_DP/abs_rij)*exp(-U_i*abs_rij)
                   fac2 = -(0.5_DP*U_i**2+U_i/abs_rij+1.0_DP/abs_rij**2) * &
                        exp(-U_i*abs_rij)
                   fac = q_i*Z_j
                   ffac = fac*fac2
                   e = fac*hlp

                   hlp = -(0.5_DP*U_j+1.0_DP/abs_rij)*exp(-U_j*abs_rij)
                   fac2 = -(0.5_DP*U_j**2+U_j/abs_rij+1.0_DP/abs_rij**2) * &
                        exp(-U_j*abs_rij)
                   fac = Z_i*q_j
                   ffac = ffac + fac*fac2
                   e = e + fac*hlp

                   !
                   ! Coulomb integrals
                   !

                   if (abs(U_i - U_j) < 1d-6) then

                      src = 1.0_DP/(U_i+U_j)
                      fac = U_i*U_j*src
                      avg = 1.6_DP*(fac+fac*fac*src)
                      fac = avg*abs_rij
                      fac2 = fac*fac
                      efac = exp(-fac)/(48*abs_rij)
                      hlp = -(48 + 33*fac + fac2*(9+fac))*efac

                      fac2 = &
                           ( hlp/abs_rij + avg*hlp &
                           + (33*avg + 18*fac*avg + 3*fac2*avg)*efac &
                           )

                      !
                      ! XH correction for DFTB2  
                      !

                      if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
                         fac2 = fac2*hij(abs_rij, U_i, U_j, this%zeta) &
                              + Sgij(abs_rij, U_i)*part_deriv_hij_wrt_r(abs_rij, U_i, U_j, this%zeta)
                      endif

                      ! 
                      ! DFTB3 forces
                      !

                      fdftb3 = 0.0_DP                      
  
                      if (this%dftb3) then

                         if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
 
                            fdftb3 = 1.0_DP/3.0_DP*(dq_i**2*dq_j*derivative_capital_short_gamma(abs_rij, dU_i, U_i, U_j, this%zeta) &
                                                  + dq_i*dq_j**2*derivative_capital_short_gamma(abs_rij, dU_j, U_j, U_i, this%zeta))
                           
                         else
           
                            fdftb3 = 1.0_DP/3.0_DP*(dq_i**2*dq_j*derivative_capital_short_gamma(abs_rij, dU_i, U_i, U_j) &
                                                  + dq_i*dq_j**2*derivative_capital_short_gamma(abs_rij, dU_j, U_j, U_i))
                        
                         endif
                           
                      endif
                          
                   else

                      fi1 = 1.0_DP/(2*(U_i**2-U_j**2)**2)
                      fj1 = -U_i**4*U_j*fi1
                      fi1 = -U_j**4*U_i*fi1
                         
                      fi2 = 1.0_DP/((U_i**2-U_j**2)**3)
                      fj2 = -(U_i**6-3*U_i**4*U_j**2)*fi2
                      fi2 =  (U_j**6-3*U_j**4*U_i**2)*fi2

                      expi = exp(-U_i*abs_rij)
                      expj = exp(-U_j*abs_rij)

                      hlp = expi*(fi1+fi2/abs_rij) + expj*(fj1+fj2/abs_rij)

                      fac2 = &
                           ( expi*( &
                           U_i*(fi1+fi2/abs_rij) + fi2/(abs_rij_sq) &
                           ) &
                           + expj*( &
                           U_j*(fj1+fj2/abs_rij) + fj2/(abs_rij_sq) &
                           ) &
                           )

                      !
                      ! XH correction for DFTB2
                      !

                      if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
                         fac2 = fac2*hij(abs_rij, U_i, U_j, this%zeta) &
                              + Sfij(abs_rij, U_i, U_j)*part_deriv_hij_wrt_r(abs_rij, U_i, U_j, this%zeta)
                      endif

                      ! 
                      ! DFTB3 forces
                      !

                      fdftb3 = 0.0_DP                      
                        
                      if (this%dftb3) then

                         if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
 
                            fdftb3 = 1.0_DP/3.0_DP*(dq_i**2*dq_j*derivative_capital_short_gamma(abs_rij, dU_i, U_i, U_j, this%zeta) &
                                                  + dq_i*dq_j**2*derivative_capital_short_gamma(abs_rij, dU_j, U_j, U_i, this%zeta))
                           
                         else
           
                            fdftb3 = 1.0_DP/3.0_DP*(dq_i**2*dq_j*derivative_capital_short_gamma(abs_rij, dU_i, U_i, U_j) &
                                                  + dq_i*dq_j**2*derivative_capital_short_gamma(abs_rij, dU_j, U_j, U_i))
                        
                         endif
                           
                      endif
                          
                   endif

                   !
                   ! DFTB3 energy 
                   !
                
                   edftb3 = 0.0_DP
                   
                   if (this%dftb3) then
                      
                      if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
                      
                         edftb3 = - 1.0_DP/3.0_DP*(dq_i**2*dq_j*capital_short_gamma(abs_rij, dU_i, U_i, U_j, this%zeta) &
                                                 + dq_i*dq_j**2*capital_short_gamma(abs_rij, dU_j, U_j, U_i, this%zeta))
                      
                      else
                      
                         edftb3 = - 1.0_DP/3.0_DP*(dq_i**2*dq_j*capital_short_gamma(abs_rij, dU_i, U_i, U_j) &
                                                 + dq_i*dq_j**2*capital_short_gamma(abs_rij, dU_j, U_j, U_i))
                                         
                      endif

                   endif

                   !
                   ! XH correction for DFTB2
                   !

                   if (this%damp_gamma .and. (atomic_number_i == 1 .or. atomic_number_j == 1)) then
                      hlp = hlp*hij(abs_rij, U_i, U_j, this%zeta)
                   endif

                   fac = q_i*q_j-q_i*Z_j-Z_i*q_j
                   ffac = ffac + fac*fac2 + fdftb3
                   df = ffac * rij/abs_rij
                   e = e + fac*hlp + edftb3

                   VEC3(f, i) = VEC3(f, i) + df
                   VEC3(f, j) = VEC3(f, j) - df

                   if (j > p%natloc) then
                      epot  = epot + 0.5_DP*e
                      wpot  = wpot - outer_product(rij, 0.5_DP*df)
                   else
                      epot  = epot + e
                      wpot  = wpot - outer_product(rij, df)
                   endif
                      
                endif
             endif

          enddo Slater_ni_loop

          epot   = epot + 5*q_i*q_i*U_i/32

          !
          ! DFTB3 on-site energy 
          !          
 
          if (this%dftb3) epot = epot - q_i**3*dU_i/6.0_DP
                  
       endif

    enddo

    call timer_stop('slater_charges_energy_and_forces')

  endsubroutine slater_charges_energy_and_forces


  !!>
  !! Registry
  !!
  !! Expose parameters to the user
  !!<
  subroutine slater_charges_register(this, cfg, m)
    implicit none

    type(slater_charges_t), target, intent(inout) :: this
    type(c_ptr),                    intent(in)    :: cfg
    type(c_ptr),                    intent(out)   :: m

    ! ---

    type(c_ptr) :: subm

    integer :: i

    ! ---

    m = ptrdict_register_section(cfg, CSTR("SlaterCharges"), &
         CSTR("This module assigns a Slater (exponential) shape to each charge and compute coulomb and nuclear repulsion integrals."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements."))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), &
         CSTR("cutoff"), &
         CSTR("Cutoff of the correction to the Coulomb potential."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, SLATER_CHARGES_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_boolean_property(m, c_loc(this%dftb3), &
         CSTR("dftb3"), &
         CSTR("Enable DFTB3's nonlinear Hubbard U."))

    call ptrdict_register_boolean_property(m, c_loc(this%damp_gamma), &
         CSTR("damp_gamma"), &
         CSTR("Enable DFTB3's damping of interaction with Hydrogen atoms."))
    call ptrdict_register_real_property(m, c_loc(this%zeta), &
         CSTR("zeta"), &
         CSTR("Exponent of damping function."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%U), SLATER_CHARGES_MAX_EL, c_loc(this%db%nU), &
         CSTR("U"), CSTR("Hubbard U."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Z), SLATER_CHARGES_MAX_EL, c_loc(this%db%nZ), &
         CSTR("Z"), CSTR("Nuclear charge."))

    subm = ptrdict_register_section(m, CSTR("HubbardDerivatives"), &
         CSTR("Derivative of Hubbard-U for DFTB3."))

    do i = 1, 116
       call ptrdict_register_real_property(subm, c_loc(this%db%dU(i)), &
            CSTR(trim(ElementName(i))), &
            CSTR("Derivative of Hubbard-U for element "//trim(ElementName(i))//"."))
    enddo

  endsubroutine slater_charges_register

endmodule slater_charges
