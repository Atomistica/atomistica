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
!   classtype:tabulated_alloy_eam_t classname:TabulatedAlloyEAM interface:potentials
!   features:mask,per_at
! @endmeta

!#define AVOID_SQRT

!>
!! The tabulated alloy EAM potential
!!
!! The tabulated alloy EAM potential.
!!
!! Usage example:
!!
!!   type(tabulated_alloy_eam_t)  :: pot
!!   ...
!!   call init(pot, fn="Au.eam")
!!   ...
!!   call energy_and_forces(pot, ...)
!!   ...
!!   call del(pot)
!<

#include "macros.inc"
#include "filter.inc"

module tabulated_alloy_eam
  use supplib

  use particles
  use neighbors

  use filter

  implicit none

  private

  integer, parameter  :: MAX_EAM_ELS = 10

  !>
  !! Tabulated_alloy_eam potential class
  !<
  public :: tabulated_alloy_eam_t
  type tabulated_alloy_eam_t

     !
     ! Element on which to apply the force
     !

     character(MAX_EL_STR)  :: elements = "*"
     integer                :: els

     character(100)         :: fn = "default.in"

     logical(BOOL)          :: dump = .false.

     !
     ! Additional information
     !

     integer                :: db_nel
     character(MAX_EL_STR)  :: db_elements(MAX_EAM_ELS)

     !
     ! Convert internal element numbers to db numbers
     !

     integer, allocatable   :: el2db(:)

     !
     ! Splines
     !

     type(simple_spline_t), allocatable  :: fF(:)      !< Embedding function
     type(simple_spline_t), allocatable  :: fphi(:, :) !< phi(r) - repulsive part
     type(simple_spline_t), allocatable  :: frho(:)    !< rho(r) - embedding density


     real(DP)               :: cutoff    !< Cut-off radius

  endtype tabulated_alloy_eam_t


  public :: init
  interface init
     module procedure tabulated_alloy_eam_init
  endinterface

  public :: del
  interface del
     module procedure tabulated_alloy_eam_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure tabulated_alloy_eam_bind_to
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure tabulated_alloy_eam_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure tabulated_alloy_eam_register
  endinterface register

  !
  ! Internal use
  !

  public :: energy_and_forces_kernel
  interface energy_and_forces_kernel
     module procedure tabulated_alloy_eam_energy_and_forces_kernel
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize the Tabulated_alloy_eam potential
  !<
  subroutine tabulated_alloy_eam_init(this, elements, fn, ierror)
    implicit none

    type(tabulated_alloy_eam_t), intent(inout)  :: this
    character(*),      optional, intent(in)     :: elements
    character(*),      optional, intent(in)     :: fn
    integer,           optional, intent(out)    :: ierror

    ! ---

    integer        :: f
    character(256) :: line, lattice

    integer        :: i, j, nF, nr, Z
    real(DP)       :: dF, dr, cutoff, mass, a0

    ! ---

    INIT_ERROR(ierror)

    call prlog("- tabulated_alloy_eam_init -")

    call del(this)

    ASSIGN_PROPERTY(elements)
    ASSIGN_PROPERTY(fn)

#ifdef AVOID_SQRT
    call prlog("     * Avoiding sqrt evaluation.")
#endif

    f = fopen(this%fn, mode=F_READ)

    call prlog("     Parameters from setfl file.")
    call prlog("     Comment:")
    call prlog("     " // trim(read_line(f)))
    call prlog("     " // trim(read_line(f)))
    call prlog("     " // trim(read_line(f)))

    line = read_line(f)
    read (line, *)  this%db_nel
    read (line, *)  this%db_nel, this%db_elements(1:this%db_nel)

    line = read_line(f)
    read (line, *)  nF, dF, nr, dr, cutoff

    call prlog("     nF     = " // nF)
    call prlog("     dF     = " // dF)
    call prlog("     nr     = " // nr)
    call prlog("     dr     = " // dr)
    call prlog("     cutoff = " // cutoff)

    if (allocated(this%fF)) then
       deallocate(this%fF)
    endif
    if (allocated(this%fphi)) then
       deallocate(this%fphi)
    endif
    if (allocated(this%frho)) then
       deallocate(this%frho)
    endif

    allocate(this%fF(this%db_nel))
    allocate(this%fphi(this%db_nel, this%db_nel))
    allocate(this%frho(this%db_nel))

    this%cutoff = cutoff

    do i = 1, this%db_nel
       line = read_line(f)
       read (line, *)  Z, mass, a0, lattice
       if (trim(ElementName(Z)) /= trim(this%db_elements(i))) then
          RAISE_ERROR("Atomic number "//Z//" is element "//trim(ElementName(Z))//" and does not match element name "//trim(this%db_elements(i))//" found in data file.", ierror)
       endif
       call read(this%fF(i), f, nF, 0.0_DP, dF)
       call read(this%frho(i), f, nr, 0.0_DP, dr, pad=(/0.0_DP,0.0_DP/))

       if (this%dump) then
          call write(this%fF(i), "fF_"//trim(this%db_elements(i))//".out", &
               0.0001_DP)
          call write(this%frho(i), "frho_"//trim(this%db_elements(i))//".out", &
               0.01_DP)
       endif

#ifdef AVOID_SQRT
       call square_x_axis(this%frho(i), 10*this%frho(i)%n)
#endif
    enddo

    do i = 1, this%db_nel
       do j = 1, i
          call read(this%fphi(i, j), f, nr, 0.0_DP, dr, pad=(/0.0_DP,0.0_DP/))

          if (this%dump) then
             call write(this%fphi(i, j), "fphi_"//trim(this%db_elements(i))// &
                  "-"//trim(this%db_elements(j))//".out", 0.01_DP)
          endif

          call scale_y_axis(this%fphi(i, j), 0.5_DP)
#ifdef AVOID_SQRT
          call square_x_axis(this%fphi(i, j), 10*this%fphi(i, j)%n)
#endif
          if (i /= j) then
             call associate(this%fphi(j, i), this%fphi(i, j))
          endif
       enddo
    enddo

    call fclose(f)

    call prlog

  endsubroutine tabulated_alloy_eam_init


  !>
  !! Destructor
  !!
  !! Free all resources occupied by this module
  !<
  subroutine tabulated_alloy_eam_del(this)
    implicit none

    type(tabulated_alloy_eam_t), intent(inout)  :: this

    ! ---

    if (allocated(this%fF)) then
       call del(this%fF)
       deallocate(this%fF)
    endif

    if (allocated(this%fphi)) then
       call del(this%fphi)
       deallocate(this%fphi)
    endif

    if (allocated(this%frho)) then
       call del(this%frho)
       deallocate(this%frho)
    endif

  endsubroutine tabulated_alloy_eam_del


  !>
  !! Constructor
  !!
  !! Initialize the Tabulated_alloy_eam potential
  !<
  subroutine tabulated_alloy_eam_bind_to(this, p, nl, ierror)
    implicit none

    type(tabulated_alloy_eam_t), intent(inout)  :: this
    type(particles_t),           intent(inout)  :: p
    type(neighbors_t),           intent(inout)  :: nl
    integer,           optional, intent(out)    :: ierror

    ! ---

    integer      :: i, j
    character(3) :: sym

    ! ---

    INIT_ERROR(ierror)

    call prlog("- tabulated_alloy_eam_bind_to -")

    this%els = filter_from_string(this%elements, p)

    if (allocated(this%el2db)) then
       deallocate(this%el2db)
    endif
    allocate(this%el2db(p%nel))
    this%el2db = -1

    do i = 1, p%nel
       if (p%el2Z(i) <= MAX_Z) then
          sym = ElementName(p%el2Z(i))
          do j = 1, this%db_nel
             if (trim(sym) == trim(this%db_elements(j))) then
                this%el2db(i) = j
             endif
          enddo
          if (this%el2db(i) /= -1) then
             call prlog("     Found parameter set for element '" // sym // "'.")
          endif
       else
          RAISE_ERROR("Unknown element with atomic number " // p%el2Z(i) // " encountered.", ierror)
       endif
       do j = i, p%nel
          if (IS_EL2(this%els, j)) then
             call request_interaction_range(nl, this%cutoff, i, j)
#ifdef LAMMPS
             call set_interaction_range(p, 2*this%cutoff, i, j)
#endif
          endif
       enddo
    enddo

    call prlog

  endsubroutine tabulated_alloy_eam_bind_to


  !>
  !! Compute energy and force
  !!
  !! Wrapper for energy and force computation. Computes size of internal
  !! neighbors list from maximum coordination numbers and passes it to the
  !! kernel. Local (per-atom) neighbor list is kept on stack.
  !<
  subroutine tabulated_alloy_eam_energy_and_forces(this, p, nl, epot, &
       f, wpot, mask, epot_per_at, wpot_per_at, ierror)
    implicit none

    type(tabulated_alloy_eam_t), intent(inout) :: this
    type(particles_t),           intent(in)    :: p
    type(neighbors_t),           intent(inout) :: nl
    real(DP),                    intent(inout) :: epot
    real(DP),                    intent(inout) :: f(3, p%nat)
    real(DP),                    intent(inout) :: wpot(3, 3)
    integer,           optional, intent(in)    :: mask(p%nat)
    real(DP),          optional, intent(inout) :: epot_per_at(p%nat)
#ifdef LAMMPS
    real(DP),          optional, intent(inout) :: wpot_per_at(6, p%nat)
#else
    real(DP),          optional, intent(inout) :: wpot_per_at(3, 3, p%nat)
#endif
    integer,           optional, intent(inout) :: ierror

    ! ---

    integer :: i, eli, maxneb

    ! ---

    call timer_start("tabulated_alloy_eam_energy_and_forces")

    INIT_ERROR(ierror)

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    maxneb = 0
#ifndef __GFORTRAN__
    !$omp  parallel do default(none) &
    !$omp  private(eli) &
    !$omp& shared(mask, nl, p, this) &
    !$omp& reduction(max:maxneb)
#endif
    do i = 1, p%natloc
       if (.not. present(mask) .or. mask(i) /= 0) then
          eli = p%el(i)
          if (IS_EL2(this%els, eli) .and. this%el2db(eli) > 0) then
             maxneb = max(maxneb, int(nl%last(i)-nl%seed(i)+1))
          endif
       endif
    enddo

    call energy_and_forces_kernel(this, p, nl, epot, f, wpot, maxneb, &
       mask=mask, epot_per_at=epot_per_at, wpot_per_at=wpot_per_at, &
       ierror=ierror)
    PASS_ERROR(ierror)

    call timer_stop("tabulated_alloy_eam_energy_and_forces")

  endsubroutine tabulated_alloy_eam_energy_and_forces


  !>
  !! Compute energy and forces
  !!
  !! Compute energy and forces
  !<
  subroutine tabulated_alloy_eam_energy_and_forces_kernel(this, p, nl, epot, &
       f, wpot, maxneb, mask, epot_per_at, wpot_per_at, ierror)
    implicit none

    type(tabulated_alloy_eam_t), intent(inout) :: this
    type(particles_t),           intent(in)    :: p
    type(neighbors_t),           intent(inout) :: nl
    real(DP),                    intent(inout) :: epot
    real(DP),                    intent(inout) :: f(3, p%nat)
    real(DP),                    intent(inout) :: wpot(3, 3)
    integer,                     intent(in)    :: maxneb
    integer,           optional, intent(in)    :: mask(p%nat)
    real(DP),          optional, intent(inout) :: epot_per_at(p%nat)
#ifdef LAMMPS
    real(DP),          optional, intent(inout) :: wpot_per_at(6, p%nat)
#else
    real(DP),          optional, intent(inout) :: wpot_per_at(3, 3, p%nat)
#endif
    integer,           optional, intent(inout) :: ierror

    ! ---

    integer             :: i, j, eli, elj, dbi, dbj, els
    integer(NEIGHPTR_T) :: ni, seedi, lasti

    real(DP)  :: dr(3), abs_dr, r_abs_dr, ri(3), fori(3)
    real(DP)  :: rho, drho, Fi, dFi
    real(DP)  :: e, w(3, 3), wij(3, 3), cutoff_sq

    real(DP)  :: phi, dphi, fac
    real(DP)  :: df(3)

    ! Internal local neighbor list cache
    integer   :: neb_n, neb(maxneb)
    real(DP)  :: neb_dr(3, maxneb), neb_abs_dr(maxneb)

    ! ---

    INIT_ERROR(ierror)

    els        = this%els
    cutoff_sq  = this%cutoff**2

    !
    ! Compute densities
    !

    e  = 0.0_DP
    w  = 0.0_DP

    !$omp  parallel default(none) &
    !$omp& firstprivate(cutoff_sq, els) &
    !$omp& private(abs_dr, r_abs_dr, eli, elj, df, Fi, dFi) &
    !$omp& private(dphi, dr, fori, rho, drho, fac) &
    !$omp& private(i, j, ni, phi, ri, seedi, lasti) &
    !$omp& private(neb_n, neb, neb_dr, neb_abs_dr) &
    !$omp& private(dbi, dbj, wij) &
    !$omp& shared(mask, nl, f, p) &
    !$omp& shared(epot_per_at, wpot_per_at, this) &
    !$omp& reduction(+:e) reduction(+:w)

    call tls_init(p%nat, sca=1, vec=1)

    !$omp do
    do i = 1, p%natloc
       if (.not. present(mask) .or. mask(i) /= 0) then

          eli  = p%el(i)
          dbi  = this%el2db(eli)

          if (IS_EL2(els, eli) .and. dbi > 0) then
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
                dbj  = this%el2db(elj)

                if (IS_EL2(els, elj) .and. dbj > 0) then
                   dr      = GET_DRJ(p, nl, i, j, ni)
                   abs_dr  = dot_product(dr, dr)

                   if (abs_dr < cutoff_sq) then
#ifndef AVOID_SQRT
                      abs_dr  = sqrt(abs_dr)
#endif
#ifdef _OPENMP
                      drho    = func(this%frho(dbj), abs_dr)
#else
                      drho    = func(this%frho(dbj), abs_dr, ierror=ierror)
                      PASS_ERROR_WITH_INFO("Error while computing density for atom " // i // ".", ierror)
#endif
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

#ifdef _OPENMP
             call f_and_df(this%fF(dbi), rho, Fi, dFi)
#else
             call f_and_df(this%fF(dbi), rho, Fi, dFi, ierror=ierror)
             PASS_ERROR_WITH_INFO("Error evaluating the embedding energy for atom " // i // ". Is the density too large?", ierror)
#endif
             tls_sca1(i)  = tls_sca1(i) + Fi

             !
             ! Loop over neighbors and compute pair terms
             !

             fori  = 0.0_DP
             do ni = 1, neb_n
                ! Pull from neighbor list cache
                j       = neb(ni)
                elj     = p%el(j)
                dbj     = this%el2db(elj)
                dr      = neb_dr(1:3, ni)
                abs_dr  = neb_abs_dr(ni)

                !
                ! Repulsive energy and forces
                !

#ifdef _OPENMP
                call f_and_df(this%fphi(dbi, dbj), abs_dr, phi, dphi)
#else
                call f_and_df(this%fphi(dbi, dbj), abs_dr, phi, dphi, &
                  ierror=ierror)
                PASS_ERROR(ierror)
#endif
                r_abs_dr     = 1.0_DP/abs_dr
                tls_sca1(i)  = tls_sca1(i) + phi*r_abs_dr

                !
                ! Forces due to embedding
                !

#ifdef _OPENMP
                fac  = dfunc(this%frho(dbj), abs_dr)
#else
                fac  = dfunc(this%frho(dbj), abs_dr, ierror=ierror)
                PASS_ERROR(ierror)
#endif
                df   = - ( dFi * fac + (dphi-phi*r_abs_dr)*r_abs_dr )* &
                     r_abs_dr * dr

                fori              = fori              + df
                VEC3(tls_vec1, j) = VEC3(tls_vec1, j) - df
                wij               = - outer_product(dr, df)
                w                 = w + wij

                if (present(wpot_per_at)) then
                   wij = wij/2
                   SUM_VIRIAL(wpot_per_at, i, wij)
                   SUM_VIRIAL(wpot_per_at, j, wij)
                endif
             enddo
             VEC3(tls_vec1, i)  = VEC3(tls_vec1, i) + fori

          endif
       endif
    enddo

    e  = e + sum(tls_sca1(1:p%natloc))

    if (present(epot_per_at)) then
       call tls_reduce(p%nat, sca1=epot_per_at, vec1=f)
    else
       call tls_reduce(p%nat, vec1=f)
    endif

    !$omp end parallel

    epot  = epot + e
    wpot  = wpot + w

  endsubroutine tabulated_alloy_eam_energy_and_forces_kernel


  !>
  !! Registry
  !!
  !! Queries parameters of the potential before initialization.
  !<
  subroutine tabulated_alloy_eam_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(tabulated_alloy_eam_t), target, intent(inout)  :: this
    type(c_ptr),                         intent(in)     :: cfg
    type(c_ptr),                         intent(out)    :: m

    ! ---

    this%elements  = "*"

    m = ptrdict_register_section(cfg, CSTR("TabulatedAlloyEAM"), &
         CSTR("General tabulated EAM potential for many component systems (alloys), see S.M. Foiles, M.I. Baskes, M.S. Daw, Phys. Rev. B 33, 7983 (1986)."))

    call ptrdict_register_string_property(m, c_locs(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Element for which to use this potential."))

    call ptrdict_register_string_property(m, c_locs(this%fn), 100, CSTR("fn"), &
         CSTR("Configuration file."))

    call ptrdict_register_boolean_property(m, c_loc(this%dump), CSTR("dump"), &
         CSTR("Dump interatomic potential to disk."))

  endsubroutine tabulated_alloy_eam_register

endmodule tabulated_alloy_eam
