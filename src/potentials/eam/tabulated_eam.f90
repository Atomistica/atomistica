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
!   classtype:tabulated_eam_t classname:TabulatedEAM interface:potentials
!   features:per_at
! @endmeta

!#define AVOID_SQRT

!>
!! The Tabulated_eam EAM potential
!!
!! The Tabulated_eam EAM potential.
!!
!! Usage example:
!!
!!   type(tabulated_eam_t)  :: pot
!!   ...
!!   call init(pot, db = Cleri_PRB_48_22_Al_Ag_Au)
!!   ...
!!   call energy_and_forces(pot, ...)
!!   ...
!!   call del(pot)
!<

#include "macros.inc"
#include "filter.inc"
#include "spline.inc"

module tabulated_eam
  use supplib

  use particles
  use neighbors

  use filter

  implicit none

  private

  !>
  !! TabulatedEAM potential class
  !<
  public :: tabulated_eam_t
  type tabulated_eam_t

     !
     ! Element on which to apply the force
     !

     character(MAX_EL_STR)  :: elements = "*"
     integer                :: els

     character(100)         :: fn = "default.in"

     !
     ! Additional information
     !

     character(100)         :: comment
     integer                :: Z         !< Element number
     real(DP)               :: mass      !< Element mass
     real(DP)               :: a0        !< Element lattice constant
     character(100)         :: lattice   !< Element ground-state

     !
     ! Splines
     !

     type(simple_spline_t)  :: fF        !< Embedding function
     type(simple_spline_t)  :: fphi      !< phi(r) - repulsive part
     type(simple_spline_t)  :: frho      !< rho(r) - embedding density
     

     real(DP)               :: cutoff    !< Cut-off radius

  endtype tabulated_eam_t


  public :: init
  interface init
     module procedure tabulated_eam_init
  endinterface

  public :: del
  interface del
     module procedure tabulated_eam_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure tabulated_eam_bind_to
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure tabulated_eam_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure tabulated_eam_register
  endinterface register

  !
  ! Internal use
  !

  public :: energy_and_forces_kernel
  interface energy_and_forces_kernel
     module procedure tabulated_eam_energy_and_forces_kernel
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize the Tabulated_eam potential
  !<
  subroutine tabulated_eam_init(this, elements, fn, ierror)
    implicit none

    type(tabulated_eam_t),           intent(inout)  :: this
    character(*),          optional, intent(in)     :: elements
    character(*),          optional, intent(in)     :: fn
    integer,               optional, intent(out)    :: ierror

    ! ---

    integer        :: f
    character(100) :: line

    integer        :: nF, nr
    real(DP)       :: dF, dr, cutoff

    ! ---

    INIT_ERROR(ierror)

    call prlog("- tabulated_eam_init -")

    call del(this)

    ASSIGN_PROPERTY(elements)
    ASSIGN_PROPERTY(fn)

#ifdef AVOID_SQRT
    call prlog("     * Avoiding sqrt evaluation.")
#endif

    f = fopen(this%fn, mode=F_READ)

    this%comment = read_line(f)

    call prlog("     Comment: " // trim(this%comment))

    line = read_line(f)
    read (line, *)  this%Z, this%mass, this%a0, this%lattice
    line = read_line(f)
    read (line, *)  nF, dF, nr, dr, cutoff

    call prlog("     nF            = " // nF)
    call prlog("     dF            = " // dF)
    call prlog("     nr            = " // nr)
    call prlog("     dr            = " // dr)
    call prlog("     cutoff        = " // cutoff)

    this%cutoff     = cutoff

    call read(this%fF, f, nF, dF, dF)
    call read(this%fphi, f, nr, dr, dr)
    call read(this%frho, f, nr, dr, dr)

    call prlog("     cutoff(fphi)  = " // this%fphi%cut)
    call prlog("     cutoff(frho)  = " // this%frho%cut)

    call write(this%fF, "fF.out", 0.0001_DP)
    call write(this%fphi, "fphi.out", 0.01_DP)
    call write(this%frho, "frho.out", 0.01_DP)

    call scale_y_axis(this%fphi, 0.5_DP)

#ifdef AVOID_SQRT
    call square_x_axis(this%fphi, 10*this%fphi%n)
    call square_x_axis(this%frho, 10*this%frho%n)

    call write(this%fphi, "fphi_sq.out", 0.01_DP)
    call write(this%frho, "frho_sq.out", 0.01_DP)
#endif
         
    call fclose(f)

  endsubroutine tabulated_eam_init


  !>
  !! Destructor
  !!
  !! Free all resources occupied by this module
  !<
  subroutine tabulated_eam_del(this)
    implicit none

    type(tabulated_eam_t), intent(inout)  :: this

    ! ---

    call del(this%fF)
    call del(this%fphi)
    call del(this%frho)

  endsubroutine tabulated_eam_del


  !>
  !! Constructor
  !!
  !! Initialize the Tabulated_eam potential
  !<
  subroutine tabulated_eam_bind_to(this, p, nl, ierror)
    implicit none

    type(tabulated_eam_t),           intent(inout)  :: this
    type(particles_t),               intent(inout)  :: p
    type(neighbors_t),               intent(inout)  :: nl
    integer,               optional, intent(out)    :: ierror

    ! ---

    integer :: i, j

    ! ---

    INIT_ERROR(ierror)

    this%els = filter_from_string(this%elements, p)

    do i = 1, p%nel
       do j = i, p%nel
          if (IS_EL2(this%els, i) .and. IS_EL2(this%els, j)) then
             call request_interaction_range(nl, this%cutoff, i, j)
          endif
       enddo
    enddo

  endsubroutine tabulated_eam_bind_to


  !>
  !! Compute energy and force
  !!
  !! Wrapper for energy and force computation. Computes size of internal
  !! neighbors list from maximum coordination numbers and passes it to the
  !! kernel. Local (per-atom) neighbor list is kept on stack.
  !<
  subroutine tabulated_eam_energy_and_forces(this, p, nl, epot, &
       f, wpot, epot_per_at, wpot_per_at, ierror)
    implicit none

    type(tabulated_eam_t), intent(inout) :: this
    type(particles_t),     intent(in)    :: p
    type(neighbors_t),     intent(inout) :: nl
    real(DP),              intent(inout) :: epot
    real(DP),              intent(inout) :: f(3, p%nat)
    real(DP),              intent(inout) :: wpot(3, 3)
    real(DP),    optional, intent(inout) :: epot_per_at(p%nat)
#ifdef LAMMPS
    real(DP),    optional, intent(inout) :: wpot_per_at(6, p%nat)
#else
    real(DP),    optional, intent(inout) :: wpot_per_at(3, 3, p%nat)
#endif
    integer,     optional, intent(inout) :: ierror

    ! ---

    integer :: i, maxneb

    ! ---

    call timer_start("tabulated_eam_energy_and_forces")

    INIT_ERROR(ierror)

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    maxneb = 0
    !$omp  parallel do default(none) &
    !$omp& shared(nl, p, this) &
    !$omp& reduction(+:maxneb)
    do i = 1, p%natloc
       if (IS_EL2(this%els, p%el(i))) then
          maxneb = max(maxneb, nl%last(i)-nl%seed(i)+1)
       endif
    enddo

    call energy_and_forces_kernel(this, p, nl, epot, f, wpot, maxneb, &
       epot_per_at=epot_per_at, wpot_per_at=wpot_per_at, ierror=ierror)
    PASS_ERROR(ierror)

    call timer_stop("tabulated_eam_energy_and_forces")

  endsubroutine tabulated_eam_energy_and_forces


  !>
  !! Compute energy and force
  !!
  !! Compute energy and force
  !<
  subroutine tabulated_eam_energy_and_forces_kernel(this, p, nl, epot, &
       f, wpot, maxneb, epot_per_at, wpot_per_at, ierror)
    implicit none

    type(tabulated_eam_t), intent(inout) :: this
    type(particles_t),     intent(in)    :: p
    type(neighbors_t),     intent(inout) :: nl
    real(DP),              intent(inout) :: epot
    real(DP),              intent(inout) :: f(3, p%maxnatloc)
    real(DP),              intent(inout) :: wpot(3, 3)
    integer,               intent(in)    :: maxneb
    real(DP),    optional, intent(inout) :: epot_per_at(p%maxnatloc)
#ifdef LAMMPS
    real(DP),    optional, intent(inout) :: wpot_per_at(6, p%maxnatloc)
#else
    real(DP),    optional, intent(inout) :: wpot_per_at(3, 3, p%maxnatloc)
#endif
    integer,     optional, intent(inout) :: ierror

    ! ---

    integer   :: i, ni, j, eli, elj, els, seedi, lasti
    real(DP)  :: dr(3), abs_dr, ri(3), fori(3)
    real(DP)  :: rho, drho, Fi, dFi
    real(DP)  :: e, w(3, 3), cutoff_sq

    real(DP)  :: phi(maxneb), dphi(maxneb), fac(maxneb)
    real(DP)  :: df(1:3, maxneb)

    ! Internal local neighbor list cache
    integer   :: neb_n, neb(maxneb)
    real(DP)  :: neb_dr(3, maxneb), neb_abs_dr(maxneb)

    ! ---

    SPLINE_INLINE
    SPLINE_INLINE_ARRAY(maxneb)
    SPLINE_INLINE_DEFINE(F, this%fF)
    SPLINE_INLINE_DEFINE(phi, this%fphi)
    SPLINE_INLINE_DEFINE(rho, this%frho)

    INIT_ERROR(ierror)

    SPLINE_INLINE_PREPARE(F, this%fF)
    SPLINE_INLINE_PREPARE(phi, this%fphi)
    SPLINE_INLINE_PREPARE(rho, this%frho)

    els        = this%els
    cutoff_sq  = this%cutoff**2

    !
    ! Compute densities
    !

    e  = 0.0_DP
    w  = 0.0_DP

    !$omp  parallel default(none) &
    !$omp& firstprivate(cutoff_sq, els) &
    !$omp& private(abs_dr, eli, elj, df, Fi, dFi) &
    !$omp& private(dphi, dr, fori, rho, drho, fac) &
    !$omp& private(i, j, ni, phi, ri, seedi, lasti) &
    !$omp& private(neb_n, neb, neb_dr, neb_abs_dr) &
    !$omp& shared(nl, f, p) &
    !$omp& shared(epot_per_at, this) &
    !$omp& SPLINE_INLINE_OMP &
    !$omp& SPLINE_INLINE_ARRAY_OMP &
    !$omp& SPLINE_INLINE_OMP_DEFINE(F) &
    !$omp& SPLINE_INLINE_OMP_DEFINE(phi) &
    !$omp& SPLINE_INLINE_OMP_DEFINE(rho) &
    !$omp& reduction(+:e) reduction(+:w)

    call tls_init(p%nat, sca=1)

    !$omp do
    do i = 1, p%natloc
       eli  = p%el(i)

       if (IS_EL2(els, eli)) then
          ri     = PNC3(p, i)

          !
          ! Compute embedding density
          !

          rho    = 0.0_DP

          seedi  = nl%seed(i)
          lasti  = nl%last(i)

          neb_n  = 0
          do ni = seedi, lasti
             j    = GET_NEIGHBOR(nl, ni)
             elj  = p%el(j)

             if (IS_EL2(els, elj)) then
                dr      = GET_DRJ(p, nl, i, j, ni)
                abs_dr  = dot_product(dr, dr)

                if (abs_dr < cutoff_sq) then
#ifndef AVOID_SQRT
                   abs_dr  = sqrt(abs_dr)
#endif
                   SPLINE_FUNC(rho, abs_dr, drho)

                   rho     = rho + drho

                   ! Add to neighbor list cache
                   neb_n               = neb_n + 1
                   neb(neb_n)          = j
                   neb_dr(1:3, neb_n)  = dr
                   neb_abs_dr(neb_n)   = abs_dr
                endif

             endif
          enddo

          !
          ! Embedding energy
          !

          SPLINE_F_AND_DF(F, rho, Fi, dFi)
          tls_sca1(i)  = tls_sca1(i) + Fi
        
          !
          ! Repulsive energy and forces
          !

          SPLINE_F_AND_DF_ARRAY(phi, 1:neb_n, neb_abs_dr, phi, dphi)
          tls_sca1(i)  = tls_sca1(i) + sum(phi(1:neb_n))

          !
          ! Forces due to embedding
          !

          SPLINE_DFUNC_ARRAY(rho, 1:neb_n, neb_abs_dr, fac)
          df(1:3, 1:neb_n)  = &
               - spread( ( dFi * fac(1:neb_n) + dphi(1:neb_n) )/neb_abs_dr(1:neb_n), dim=1, ncopies=3 ) &
               * neb_dr(1:3, 1:neb_n)

          fori  = 0.0_DP
          do ni = 1, neb_n
             j           = neb(ni)
             fori        = fori       + df(1:3, ni)
             VEC3(f, j)  = VEC3(f, j) - df(1:3, ni)
             w           = w + (- outer_product(neb_dr(1:3, ni), df(1:3, ni)))
          enddo
          VEC3(f, i)  = VEC3(f, i) + fori

       endif
    enddo

    e  = e + sum(tls_sca1(1:p%natloc))

    if (present(epot_per_at)) then
       call tls_reduce(p%nat, sca1=epot_per_at)
    endif

    !$omp end parallel

    w = w/2

    epot  = epot + e
    wpot  = wpot + w

  endsubroutine tabulated_eam_energy_and_forces_kernel


  subroutine tabulated_eam_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(tabulated_eam_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)               :: cfg
    type(c_ptr), intent(out)              :: m

    ! ---

    this%elements  = "*"

    m = ptrdict_register_section(cfg, CSTR("TabulatedEAM"), &
         CSTR("General tabulated EAM potential, see S.M. Foiles, M.I. Baskes, M.S. Daw, Phys. Rev. B 33, 7983 (1986)."))

    call ptrdict_register_string_property(m, c_locs(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Element for which to use this potential."))

    call ptrdict_register_string_property(m, c_locs(this%fn), 100, CSTR("fn"), &
         CSTR("Configuration file."))

  endsubroutine tabulated_eam_register

endmodule tabulated_eam
