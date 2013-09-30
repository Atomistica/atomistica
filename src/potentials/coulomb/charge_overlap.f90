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
!! Assigns a shape to the charges on the atoms. Currently Gaussians and
!! Slater-type (exponential) shapes are supported.
!!
!! This module is required for both the TightBinding and VariableCharge modules.
!<

#include "macros.inc"
#include "filter.inc"

module charge_overlap
  use libAtoms_module

  use c_f

  use logging
  use timer
  use tls

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

  !
  ! The module for the computation of energies/potentials
  !

  public :: charge_overlap_t
  type charge_overlap_t

     integer                :: shape = NONE

     !
     ! Element on which to apply the force
     !

     character(MAX_EL_STR)  :: elements = "*"
     integer                :: els

     !
     ! real space sampling
     !

     real(DP)               :: cutoff = 5.0_DP
     real(DP)               :: cutoff_sq

     !
     ! other stuff
     !

!     real(DP)            :: q_tot
!     real(DP)            :: q_tot_gauss
!     real(DP)            :: jellium_potential

     !
     ! U in current units
     !

     real(DP), allocatable  :: U_in(:)
     real(DP), allocatable  :: U(:)

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

  public :: potential_and_field
  interface potential_and_field
     module procedure charge_overlap_potential_and_field
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

    type(charge_overlap_t), intent(inout)  :: this

    ! ---

    if (allocated(this%U)) then
       deallocate(this%U)
    endif

  endsubroutine charge_overlap_del


  !>
  !! Set the Hubbard U values object to this ChargeOverlap object
  !!
  !! Assign a Neighbors object to this ChargeOverlap object. All subsequent operations
  !! will use this Neighbors object. Only a pointer to the object
  !! is copied, not the object itself.
  !!
  !! Note that this needs to be called before *bind_to*!
  !<
  subroutine charge_overlap_set_Hubbard_U(this, p, U, error)
    implicit none

    type(charge_overlap_t), intent(inout) :: this  !> ChargeOverlap object
    type(particles_t),      intent(in)    :: p
    real(DP),               intent(in)    :: U(p%nel)
    integer,      optional, intent(out)   :: error

    ! ---

    INIT_ERROR(error)

    if (allocated(this%U_in))  deallocate(this%U_in)
    allocate(this%U_in(p%nel))

    this%U_in  = U

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

    type(charge_overlap_t), intent(inout)  :: this
    type(particles_t), intent(in)          :: p
    type(neighbors_t), intent(inout)       :: nl
    integer, intent(inout), optional       :: ierror

    ! ---

    integer   :: i

    ! ---

    write (ilog, '(A)')     "- charge_overlap_bind_to -"
    write (ilog, '(5X,A,A)')       "elements       = ", trim(this%elements)

    this%els  = filter_from_string(this%elements, p, ierror=ierror)
    PASS_ERROR(ierror)

    !
    ! Convert units of Hubbard U's
    !

    if (.not. allocated(this%U_in)) then
       RAISE_ERROR("charge_overlap_bind_to: No Hubbard-Us have been specified. Did you specify a tight-binding SCC or a classical variable charge model?", ierror)
    endif

    allocate(this%U(p%nel))

    this%U = this%U_in / (Hartree*Bohr)

    do i = 1, p%nel
       write (ilog, '(5X,A2,A,F20.10)')  &
            ElementName(p%el2Z(i)), ", U = ", this%U(i)
    enddo

    if (this%shape == SLATER) then
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

    type(charge_overlap_t), intent(inout)  :: this
    type(particles_t), intent(in)          :: p
    type(neighbors_t), intent(inout)       :: nl
    real(DP), intent(in)                   :: q(p%maxnatloc)
    real(DP), intent(inout)                :: phi(p%maxnatloc)
    integer, intent(inout), optional       :: ierror

    !---

    real(DP)  :: U_i_sq, abs_rij, hlp, expi, expj, src, fac, fac2, efac
    real(DP)  :: avg, fi1, fj1, fi2, fj2, U_i, U_j, q_i

    integer   :: i, ni, j

    integer, parameter :: sqrt_pi_2 = sqrt(PI/2)

    !---

    call timer_start('charge_overlap_potential')

    if (this%shape /= SELF_ENERGY_ONLY) then
       call update(nl, p, ierror)
       PASS_ERROR(ierror)
    end if

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

       !$omp  parallel default(none) &
       !$omp& shared(nl, this, phi, p, q) &
       !$omp& private(i, q_i, U_i, U_j, j, ni, abs_rij, hlp, src, fac, avg) &
       !$omp& private(fac2, efac, fi1, fi2, fj1, fj2, expi, expj)

       call tls_init(size(phi), sca=1)  ! is called tls_sca1 (=phi)

       !$omp do 
       do i = 1, p%natloc

          if (IS_EL(this%els, p, i)) then

             q_i = q(i)
             U_i = this%U(p%el(i))

             !
             ! Atom i has a Gaussian charge cloud
             !

             Slater_ni_loop: do ni = nl%seed(i), nl%last(i)
                j = nl%neighbors(ni)

                if (i <= j .and. IS_EL(this%els, p, j)) then
                   abs_rij = GET_ABS_DRJ(p, nl, i, j, ni)
                   if (abs_rij < this%cutoff) then

                      U_j = this%U(p%el(j))

                      if (abs(U_i - U_j) < 1d-6) then

                         src = 1.0_DP/(U_i+U_j)
                         fac = U_i*U_j*src
                         avg = 1.6_DP*(fac+fac*fac*src)
                         fac = avg*abs_rij
                         fac2 = fac*fac
                         efac = exp(-fac)/(48*abs_rij)

                         hlp = -(48 + 33*fac + fac2*(9+fac))*efac
                         tls_sca1(i) = tls_sca1(i) + q(j)*hlp
                         tls_sca1(j) = tls_sca1(j) + q_i*hlp

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
                         tls_sca1(i) = tls_sca1(i) + q(j)*hlp
                         tls_sca1(j) = tls_sca1(j) + q_i*hlp

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
  subroutine charge_overlap_potential_and_field(this, p, nl, q, phi, epot, E, &
       wpot, ierror)
    implicit none

    type(charge_overlap_t), intent(inout)  :: this
    type(particles_t), intent(in)          :: p
    type(neighbors_t), intent(inout)       :: nl
    real(DP), intent(in)                   :: q(p%maxnatloc)
!    real(DP), intent(inout)                :: epot
    real(DP), intent(inout)                :: phi(p%maxnatloc)
    real(DP), intent(inout)                :: epot
    real(DP), intent(inout)                :: E(3, p%maxnatloc)
    real(DP), intent(inout)                :: wpot(3, 3)
    integer, intent(inout), optional       :: ierror

    !---

    real(DP)  :: U_i_sq, q_i, rij(3), abs_rij, abs_rij_sq
    real(DP)  :: c, df(3), hlp, sqrt_pi_2, src, fac, fac2, efac, expi, expj
    real(DP)  :: avg, fi1, fj1, fi2, fj2, U_i, U_j

    integer   :: i, ni, j

    !---

    call timer_start('charge_overlap_potential_and_field')

    if (this%shape /= SELF_ENERGY_ONLY) then
       call update(nl, p, ierror)
       PASS_ERROR(ierror)
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
                            phi(i)   = phi(i) + q_i*hlp

                            epot     = epot + 0.5_DP*q_i*q_i*hlp
                            wpot     = wpot - outer_product(rij, 0.5_DP*q_i*q(j)*df)
                         else
                            phi(i)   = phi(i) + q(j)*hlp
                            phi(j)   = phi(j) + q_i*hlp

                            VEC3(E, i)  = VEC3(E, i) + q(j)*df
                            VEC3(E, j)  = VEC3(E, j) - q_i*df

                            if (j > p%natloc) then
                               epot  = epot + 0.5_DP*q_i*q(j)*hlp
                               wpot  = wpot - outer_product(rij, 0.5_DP*q_i*q(j)*df)
                            else
                               epot  = epot + q_i*q(j)*hlp
                               wpot  = wpot - outer_product(rij, q_i*q(j)*df)
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

                         if (abs_rij_sq < this%cutoff_sq) then
                            abs_rij  = sqrt(abs_rij_sq)

                            c   = PI/2*this%U(p%el(j))**2
                            hlp = -erfc(sqrt_pi_2*this%U(p%el(j))*abs_rij) &
                                 /(abs_rij)

                            df = -rij/(abs_rij**3) *  &
                                 ( erfc(sqrt(c) * abs_rij) &
                                 + 2*sqrt(c/PI)*exp(-c*abs_rij_sq)*abs_rij &
                                 )

                            phi(i)      = phi(i) + q(j)*hlp
                            phi(j)      = phi(j) + q(i)*hlp

                            VEC3(E, i)  = VEC3(E, i) + q(j)*df
                            VEC3(E, j)  = VEC3(E, j) - q_i*df

                            if (j > p%natloc) then
                               epot  = epot + 0.5_DP*q_i*q(j)*hlp
                               wpot  = wpot - outer_product(rij, 0.5_DP*q_i*q(j)*df)
                            else
                               epot  = epot + q_i*q(j)*hlp
                               wpot  = wpot - outer_product(rij, q_i*q(j)*df)
                            endif
                         endif

                      endif

                   endif
                enddo ni_loop2

             endif

             phi(i) = phi(i) + q(i)*this%U(p%el(i))! - this%jellium_potential
             epot   = epot + 0.5_DP*q(i)*q(i)*this%U(p%el(i)) 
                      
          endif

       enddo

    else if (this%shape == SLATER) then

       do i = 1, p%natloc

          if (IS_EL(this%els, p, i)) then

             q_i = q(i)
             U_i = this%U(p%el(i))

             !
             ! Atom i has a Gaussian charge cloud
             !

             Slater_ni_loop: do ni = nl%seed(i), nl%last(i)
                j   = nl%neighbors(ni)

                if (i <= j .and. IS_EL(this%els, p, j)) then
                   DIST_SQ(p, nl, i, ni, rij, abs_rij_sq)

                   if (abs_rij_sq < this%cutoff_sq) then
                      abs_rij  = sqrt(abs_rij_sq)

                      U_j = this%U(p%el(j))

                      if (abs(U_i - U_j) < 1d-6) then

                         src = 1.0_DP/(U_i+U_j)
                         fac = U_i*U_j*src
                         avg = 1.6_DP*(fac+fac*fac*src)
                         fac = avg*abs_rij
                         fac2 = fac*fac
                         efac = exp(-fac)/(48*abs_rij)
                         hlp = -(48 + 33*fac + fac2*(9+fac))*efac

                         df = ( hlp/abs_rij + avg*hlp &
                              + (33*avg + 18*fac*avg + 3*fac2*avg)*efac &
                              )*rij/abs_rij

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

                         df = ( expi*( &
                                U_i*(fi1+fi2/abs_rij) + fi2/(abs_rij_sq) &
                                ) &
                              + expj*( &
                                U_j*(fj1+fj2/abs_rij) + fj2/(abs_rij_sq) &
                                ) &
                              )*rij/abs_rij

                      endif

                      phi(i) = phi(i) + q(j)*hlp
                      phi(j) = phi(j) + q_i*hlp

                      VEC3(E, i) = VEC3(E, i) + q(j)*df
                      VEC3(E, j) = VEC3(E, j) - q_i*df

                      if (j > p%natloc) then
                         epot  = epot + 0.5_DP*q_i*q(j)*hlp
                         wpot  = wpot - outer_product(rij, 0.5_DP*q_i*q(j)*df)
                      else
                         epot  = epot + q_i*q(j)*hlp
                         wpot  = wpot - outer_product(rij, q_i*q(j)*df)
                      endif
                      
                   endif
                endif

             enddo Slater_ni_loop

             phi(i) = phi(i) + 5*q_i*U_i/16
             epot   = epot + 5*q_i*q_i*U_i/32

          endif

       enddo

    else if (this%shape == SELF_ENERGY_ONLY) then

       do i = 1, p%natloc
          if (IS_EL(this%els, p, i)) then
             phi(i)  = phi(i) + q(i)*this%U(p%el(i))! - this%jellium_potential
             epot    = epot + 0.5_DP*q(i)*q(i)*this%U(p%el(i))
          endif
       enddo
          
    else
      
       RAISE_ERROR("Fatal internal error: Unknown value (" // this%shape // ") for the shape parameter.", ierror)

    endif

    call timer_stop('charge_overlap_potential_and_field')

  endsubroutine charge_overlap_potential_and_field


  !!>
  !! Registry
  !!
  !! Expose parameters to the user
  !!<
  subroutine charge_overlap_register(this, cfg, m)
    implicit none

    type(charge_overlap_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

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

  endsubroutine charge_overlap_register

endmodule charge_overlap
