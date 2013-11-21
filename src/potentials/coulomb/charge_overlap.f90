!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
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
!   classtype:charge_overlap_t classname:ChargeOverlap interface:coulomb
! @endmeta

!>
!! Gaussian or Slater-type broadening of charges.
!!
!! Gaussian or Slater-type broadening of charges.
!!
!! Assigns a shape to the charges on the atoms. Charge distributions is given
!!
!! \f[
!!   \rho(\vec{r}) = \sum_i \left( Z_i \delta^3(\vec{r}-\vec{r_i}) + (q_i-Z_i) f_i(\vec{r}-\vec{r_i}) \right)
!! \f]
!!
!! where \f$Z_i\f$ is the effective nuclear charge and \f$f(\vec{r})\f$ the shape.
!! Currently Gaussians and Slater-type (exponential) shapes are supported.
!!
!! Note that this module does not compute the electrostatic potential and field of the
!! singular, long-ranged \f$1/r\f$ term.
!!
!! This module is required for both the TightBinding and VariableCharge modules. For
!! tight-binding calculations, \f$Z_i=0\f$ for all \f$i\f$.
!<

#include "macros.inc"
#include "filter.inc"

module charge_overlap
  use supplib

  use particles
  use neighbors
  use filter

  implicit none

  private

  integer, parameter                   :: n_shapes = 4
  integer, parameter                   :: NONE = 0
  integer, parameter                   :: SELF_ENERGY_ONLY = 1
  integer, parameter                   :: GAUSSIAN = 2
  integer, parameter                   :: SLATER = 3

  integer, parameter                   :: len_shape_str = 20
  character(len_shape_str), parameter  :: str0 = CSTR("NONE")
  character(len_shape_str), parameter  :: str1 = CSTR("self-energy-only")
  character(len_shape_str), parameter  :: str2 = CSTR("Gaussian")
  character(len_shape_str), parameter  :: str3 = CSTR("Slater")
  character(len_shape_str), parameter  :: shape_strs(n_shapes) = &
       (/ str0, str1, str2, str3 /)

  public :: n_shapes, SELF_ENERGY_ONLY, GAUSSIAN, SLATER, len_shape_str
  public :: shape_strs

  integer, parameter  :: CHARGE_TRANSFER_MAX_EL  = 3

  !
  ! The module for the computation of energies/potentials
  !

  type charge_overlap_db_t

     integer   :: nel = -1
     integer   :: nU, nZ

     character :: el(2, CHARGE_TRANSFER_MAX_EL)   !< Atom type

     real(DP)  :: U(CHARGE_TRANSFER_MAX_EL)       !< Hubbard U
     real(DP)  :: Z(CHARGE_TRANSFER_MAX_EL)       !< Nuclear charge

  endtype charge_overlap_db_t

  public :: charge_overlap_t
  type charge_overlap_t

     integer :: shape = NONE

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
     ! Effective nuclear charge
     !

     real(DP), allocatable :: Z(:)

     !
     ! Database
     !
     
     type(charge_overlap_db_t) :: db

  endtype charge_overlap_t


  public :: init
  interface init
     module procedure charge_overlap_init
  endinterface

  public :: del
  interface del
     module procedure charge_overlap_del
  endinterface

  public :: set_Hubbard_U
  interface set_Hubbard_U
     module procedure charge_overlap_set_Hubbard_U
  endinterface

  public :: bind_to
  interface bind_to
     module procedure charge_overlap_bind_to
  endinterface

  public :: potential
  interface potential
     module procedure charge_overlap_potential
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure charge_overlap_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure charge_overlap_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a ChargeOverlap object
  !<
  subroutine charge_overlap_init(this, p, U, elements, cutoff, shape, &
       error)
    use, intrinsic :: iso_c_binding

    implicit none

    type(charge_overlap_t),      intent(inout) :: this
    type(particles_t), optional, intent(in)    :: p
    real(DP),          optional, intent(in)    :: U(*)
    character(*),      optional, intent(in)    :: elements
    real(DP),          optional, intent(in)    :: cutoff
    integer,           optional, intent(in)    :: shape
    integer,           optional, intent(inout) :: error

    ! ---

    call prlog("- charge_overlap_init -")

    ASSIGN_PROPERTY(elements)
    ASSIGN_PROPERTY(cutoff)
    ASSIGN_PROPERTY(shape)

    if (this%shape == NONE) then
       RAISE_ERROR("Please specify a shape for the ChargeOverlap module.", error)
    endif

    call prlog("     shape = " // shape_strs(this%shape+1)(1:len_trim(shape_strs(this%shape+1))-1))

    if (present(U)) then
       call set_Hubbard_U(this, p, U, error=error)
       PASS_ERROR(error)
    endif

    call prlog

  endsubroutine charge_overlap_init


  !>
  !! Destructor
  !!
  !! Free all internal data buffers
  !<
  subroutine charge_overlap_del(this)
    implicit none

    type(charge_overlap_t), intent(inout) :: this

    ! ---

    if (allocated(this%U)) then
       deallocate(this%U)
    endif
    if (allocated(this%Z)) then
       deallocate(this%Z)
    endif

  endsubroutine charge_overlap_del


  !>
  !! Set the Hubbard U values and nuclear charges
  !!
  !! Set the Hubbard U values and nuclear charges Z. U and Z values are passed
  !! and stored per element, not per atom.
  !!
  !! Note that this needs to be called before *bind_to*!
  !<
  subroutine charge_overlap_set_Hubbard_U(this, p, U, Z, error)
    implicit none

    type(charge_overlap_t), intent(inout) :: this  !> ChargeOverlap object
    type(particles_t),      intent(in)    :: p
    real(DP),               intent(in)    :: U(p%nel)
    real(DP),     optional, intent(in)    :: Z(p%nel)
    integer,      optional, intent(out)   :: error

    ! ---

    INIT_ERROR(error)

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

  endsubroutine charge_overlap_set_Hubbard_U


  !>
  !! Assign a Particles and a Neighbors object to this ChargeOverlap object
  !!
  !! Assign a Particles and a Neighbors object to this ChargeOverlap object. All subsequent operations
  !! will use this Particles object. Only a pointer to the object
  !! is copied, not the object itself.
  !!
  !! Note that this needs to be called *after* set_Hubbard_U!
  !<
  subroutine charge_overlap_bind_to(this, p, nl, ierror)
    implicit none

    type(charge_overlap_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    integer,      optional, intent(inout) :: ierror

    ! ---

    integer :: i, j, Z

    ! ---

    write (ilog, '(A)')     "- charge_overlap_bind_to -"
    write (ilog, '(5X,A,A)')       "elements       = ", trim(this%elements)

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
    call resize(this%Z, p%nel)

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
                call prlog("     " // ElementName(Z) // " - " // j)
                call prlog("     - U = " // this%db%U(i) // ", Z = " // this%db%Z(i))

                this%U(j)  = this%db%U(i) / (Hartree*Bohr)
                this%Z(j)  = this%db%Z(i)
             endif
          enddo
       enddo
    else
       if (this%db%nU > 0)  this%U(1:this%db%nU) = this%db%U(1:this%db%nU)
       if (this%db%nZ > 0)  this%Z(1:this%db%nU) = this%db%Z(1:this%db%nU)
    endif

    if (this%shape == SLATER) then

       ! Note: If shape == SLATER, then U_i has been converted such that the
       ! charge is f_i(r) ~ exp(-U_i r), i.e. U_i = tau_i from Elstner's paper
       ! or U_i = 2 zeta_i from Streitz-Mintmire's paper.

       this%U = 16*this%U/5
    endif

    this%cutoff_sq = this%cutoff**2

    if (this%shape /= SELF_ENERGY_ONLY) then

       call request_interaction_range(nl, this%cutoff)

       !
       ! Fixme!!! Does only work for square geometries
       !

       write (ilog, '(5X,A,F20.10)')  "cutoff    = ", this%cutoff
    else
       write (ilog, '(5X,A)')  "No overlap integrals, will only add self-energy."
    endif

    write (ilog, *)

  endsubroutine charge_overlap_bind_to


  !>
  !! Calculate the electrostatic potential of every atom (for variable charge
  !! models)
  !!
  !! Difference between the Ewald potential and the real one due to
  !! the Gaussian charge distribution. You always need an additional
  !! DirectCoulomb, Ewald, etc.
  !<
  subroutine charge_overlap_potential(this, p, nl, q, phi, ierror)
    implicit none

    type(charge_overlap_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    real(DP),               intent(in)    :: q(p%maxnatloc)
    real(DP),               intent(inout) :: phi(p%maxnatloc)
    integer,      optional, intent(inout) :: ierror

    !---

    real(DP) :: U_i_sq, abs_rij, hlp, expi, expj, src, fac, fac2, efac
    real(DP) :: avg, fi1, fj1, fi2, fj2, U_i, U_j, q_i, q_j, Z_i, Z_j

    integer :: i, ni, j

    integer, parameter :: sqrt_pi_2 = sqrt(PI/2)

    !---

    call timer_start('charge_overlap_potential')

    if (this%shape /= SELF_ENERGY_ONLY) then
       call update(nl, p, ierror)
       PASS_ERROR(ierror)
    endif

    if (this%shape == GAUSSIAN) then

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
                   j   = nl%neighbors(ni)

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
                   j   = nl%neighbors(ni)

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

    else if (this%shape == SLATER) then

       ! Note: If shape == SLATER, then U_i has been converted such that the
       ! charge is f_i(r) ~ exp(-U_i r), i.e. U_i = tau_i from Elstner's paper
       ! or U_i = 2 zeta_i from Streitz-Mintmire's paper.

       !$omp  parallel default(none) &
       !$omp& shared(nl, this, phi, p, q) &
       !$omp& private(i, j, q_i, q_j, U_i, U_j, Z_i, Z_j) &
       !$omp& private(ni, abs_rij, hlp, src, fac, avg) &
       !$omp& private(fac2, efac, fi1, fi2, fj1, fj2, expi, expj)

       call tls_init(size(phi), sca=1)  ! is called tls_sca1 (=phi)

       !$omp do 
       do i = 1, p%natloc

          if (IS_EL(this%els, p, i)) then

             q_i = q(i)
             U_i = this%U(p%el(i))
             Z_i = this%Z(p%el(i))

             !
             ! Atom i has a Gaussian charge cloud
             !

             Slater_ni_loop: do ni = nl%seed(i), nl%last(i)
                j = nl%neighbors(ni)

                if (i <= j .and. IS_EL(this%els, p, j)) then
                   abs_rij = GET_ABS_DRJ(p, nl, i, j, ni)
                   if (abs_rij < this%cutoff) then

                      q_j = q(j)
                      U_j = this%U(p%el(j))
                      Z_j = this%Z(p%el(j))

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
                         tls_sca1(i) = tls_sca1(i) + (q_j-Z_j)*hlp
                         tls_sca1(j) = tls_sca1(j) + (q_i-Z_i)*hlp

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
                         tls_sca1(i) = tls_sca1(i) + (q_j-Z_j)*hlp
                         tls_sca1(j) = tls_sca1(j) + (q_i-Z_i)*hlp

                      endif
                      
                   endif
                endif

             enddo Slater_ni_loop

             tls_sca1(i) = tls_sca1(i) + 5*q_i*U_i/16
             
          endif

       enddo

       call tls_reduce(p%natloc, sca1=phi)
       !$omp end parallel

    else if (this%shape == SELF_ENERGY_ONLY) then

       do i = 1, p%natloc
          if (IS_EL(this%els, p, i)) then
             phi(i) = phi(i) + q(i)*this%U(p%el(i))! - this%jellium_potential
          endif
       enddo

    else
      
       RAISE_ERROR("Fatal internal error: Unknown value (" // this%shape // ") for the shape parameter.", ierror)

    endif

    call timer_stop('charge_overlap_potential')

  endsubroutine charge_overlap_potential


  !>
  !! Calculate the electrostatic potential and the electric field
  !!
  !! Difference between the Ewald (point charge) force and the force
  !! due to the Gaussian charge distribution
  !<
  subroutine charge_overlap_energy_and_forces(this, p, nl, q, epot, f, wpot, error)
    implicit none

    type(charge_overlap_t), intent(inout) :: this
    type(particles_t),      intent(in)    :: p
    type(neighbors_t),      intent(inout) :: nl
    real(DP),               intent(in)    :: q(p%maxnatloc)
    real(DP),               intent(inout) :: epot
    real(DP),               intent(inout) :: f(3, p%maxnatloc)
    real(DP),               intent(inout) :: wpot(3, 3)
    integer,      optional, intent(inout) :: error

    !---

    real(DP) :: U_i_sq, q_i, q_j, rij(3), abs_rij, abs_rij_sq
    real(DP) :: c, df(3), hlp, sqrt_pi_2, src, fac, fac2, efac, expi, expj
    real(DP) :: avg, e, ffac, fi1, fj1, fi2, fj2, U_i, U_j, Z_i, Z_j

    integer :: i, ni, j

    !---

    call timer_start('charge_overlap_energy_and_forces')

    if (this%shape /= SELF_ENERGY_ONLY) then
       call update(nl, p, error)
       PASS_ERROR(error)
    end if

    if (this%shape == GAUSSIAN) then

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
                   j = nl%neighbors(ni)

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
                   j = nl%neighbors(ni)

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

    else if (this%shape == SLATER) then

       ! Note: If shape == SLATER, then U_i has been converted such that the
       ! charge is f_i(r) ~ exp(-U_i r), i.e. U_i = tau_i from Elstner's paper
       ! or U_i = 2 zeta_i from Streitz-Mintmire's paper.

       do i = 1, p%natloc

          if (IS_EL(this%els, p, i)) then

             q_i = q(i)
             U_i = this%U(p%el(i))
             Z_i = this%Z(p%el(i))

             !
             ! Atom i has a Gaussian charge cloud
             !

             Slater_ni_loop: do ni = nl%seed(i), nl%last(i)
                j   = nl%neighbors(ni)

                if (i <= j .and. IS_EL(this%els, p, j)) then
                   DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)

                   if (abs_rij_sq < this%cutoff_sq) then
                      abs_rij  = sqrt(abs_rij_sq)

                      q_j = q(j)
                      U_j = this%U(p%el(j))
                      Z_j = this%Z(p%el(j))

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

                      endif

                      fac = q_i*q_j-q_i*Z_j-Z_i*q_j
                      ffac = ffac + fac*fac2
                      df = ffac * rij/abs_rij
                      e = e + fac*hlp

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

          endif

       enddo

    else if (this%shape == SELF_ENERGY_ONLY) then

       do i = 1, p%natloc
          if (IS_EL(this%els, p, i)) then
             epot    = epot + 0.5_DP*q(i)*q(i)*this%U(p%el(i))
          endif
       enddo
          
    else
      
       RAISE_ERROR("Fatal internal error: Unknown value '" // this%shape // "' for the shape parameter.", error)

    endif

    call timer_stop('charge_overlap_energy_and_forces')

  endsubroutine charge_overlap_energy_and_forces


  !!>
  !! Registry
  !!
  !! Expose parameters to the user
  !!<
  subroutine charge_overlap_register(this, cfg, m)
    implicit none

    type(charge_overlap_t), target, intent(inout) :: this
    type(c_ptr),                    intent(in)    :: cfg
    type(c_ptr),                    intent(out)   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("ChargeOverlap"), &
         CSTR("This module assign a shape to each charge."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements."))

    call ptrdict_register_enum_property(m, c_loc(this%shape), &
         n_shapes, len_shape_str, shape_strs, &
         CSTR("shape"), &
         CSTR("Shape of the charges. Possibilites are 'self-energy-only', 'Gaussian' and 'Slater'."))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), &
         CSTR("cutoff"), &
         CSTR("Cutoff of the correction to the Coulomb potential."))

    call ptrdict_register_string_list_property(m, &
         c_loc11(this%db%el), 2, CHARGE_TRANSFER_MAX_EL, c_loc(this%db%nel), &
         CSTR("el"), CSTR("List of element symbols."))

    call ptrdict_register_list_property(m, &
         c_loc1(this%db%U), CHARGE_TRANSFER_MAX_EL, c_loc(this%db%nU), &
         CSTR("U"), CSTR("Hubbard U."))
    call ptrdict_register_list_property(m, &
         c_loc1(this%db%Z), CHARGE_TRANSFER_MAX_EL, c_loc(this%db%nZ), &
         CSTR("Z"), CSTR("Nuclear charge."))

  endsubroutine charge_overlap_register

endmodule charge_overlap
