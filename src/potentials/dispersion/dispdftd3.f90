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
!   classtype:dispdftd3_t classname:DFTD3 interface:potentials
! @endmeta
!>
!! Grimme's DFT-D3 dispersion potential
!!
!! Reference:
!! Grimme et al., J. Chem. Phys. 132, 154104 (2010)
!! Grimme et al., J. Comp. Chem. 32, 1456 (2011)
!!
!! Requirement: 
!! dftd3-lib: git@github.com:dftbplus/dftd3-lib.git
!! 
!<

#include "macros.inc"
#include "filter.inc"

module dispdftd3
  use libAtoms_module 
  use ptrdict

  use logging
  use timer

  use particles
  use neighbors
  use filter

#ifdef HAVE_DFTD3
  use dftd3_api
#endif

  implicit none

  private

  public :: dispdftd3_t 
  type dispdftd3_t

     real(DP)  :: a1       = 0.5719_DP
     real(DP)  :: a2       = 3.6017_DP
     real(DP)  :: s6       = 1.0000_DP
     real(DP)  :: s8       = 0.5883_DP
     real(DP)  :: sr6      = 0.7461_DP
     real(DP)  :: sr8      = 1.0000_DP
     real(DP)  :: alpha6   = 14.000_DP

     real(DP)  :: cutoff   = 80.0_DP
     real(DP)  :: cutoffCN = 40.0_DP

     logical :: BeckeJohnson = .true.
     logical :: threebody    = .false.
                           
#ifdef HAVE_DFTD3
     type(dftd3_calc), allocatable  ::  calculator
#endif

  endtype dispdftd3_t

  public :: init
  interface init
     module procedure dispdftd3_init
  endinterface

  public :: del
  interface del
     module procedure dispdftd3_del
  endinterface

  public :: bind_to
  interface bind_to
     module procedure dispdftd3_bind_to
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure dispdftd3_energy_and_forces
  endinterface

  public :: register
  interface register
     module procedure dispdftd3_register
  endinterface

contains


  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine dispdftd3_init(this)
    implicit none

    type(dispdftd3_t), intent(inout)     :: this

#ifdef HAVE_DFTD3
    type(dftd3_input) :: input
#endif


    ! ---
   
#ifdef HAVE_DFTD3
    allocate(this%calculator)

    input%threebody = this%threebody 
    input%numgrad = .false.

    input%cutoff = this%cutoff * length_to_Bohr
    if (this%threebody) then
      input%cutoff_cn = this%cutoffCN * length_to_Bohr
    endif 

    call dftd3_init(this%calculator, input)

    if (this%BeckeJohnson) then
       call dftd3_set_params(this%calculator, [this%s6, this%a1, this%s8, this%a2, 0.0_DP], 4)
    else
       call dftd3_set_params(this%calculator, [this%s6, this%sr8, this%s8, this%sr8, this%alpha6], 3)
    endif
#endif

  endsubroutine dispdftd3_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine dispdftd3_del(this)
    implicit none

    type(dispdftd3_t), intent(inout)  :: this

    ! ---

#ifdef HAVE_DFTD3
    deallocate(this%calculator)
#endif

  endsubroutine dispdftd3_del


  !>
  !! Initialization
  !!
  !! Constructs the parameter sets
  !<
  subroutine dispdftd3_bind_to(this, p, nl, ierror)
    implicit none

    type(dispdftd3_t), intent(inout) :: this
    type(particles_t), intent(in)    :: p
    type(neighbors_t), intent(inout) :: nl
    integer, optional, intent(inout) :: ierror

    ! ---

#ifndef HAVE_DFTD3
    RAISE_ERROR("This version of Atomistica was not compiled with DFT-D3 support.", ierror)
#endif

    write (ilog, '(A)')  "- dftd3_bind_to -"

    write (ilog, '(5X,A,L)')      "BeckeJohnson = ", this%BeckeJohnson
    write (ilog, '(5X,A,L)')      "threebody    = ", this%threebody

    write (ilog, '(5X,A,F7.3)')  "a1       = ", this%a1
    write (ilog, '(5X,A,F7.3)')  "a2       = ", this%a2
                                           
    write (ilog, '(5X,A,F7.3)')  "s6       = ", this%s6
    write (ilog, '(5X,A,F7.3)')  "s8       = ", this%s8
                                           
    write (ilog, '(5X,A,F7.3)')  "sr6      = ", this%sr6
    write (ilog, '(5X,A,F7.3)')  "sr8      = ", this%sr8
    write (ilog, '(5X,A,F7.3)')  "alpha6   = ", this%alpha6
                                           
    write (ilog, '(5X,A,F7.3)')  "cutoff   = ", this%cutoff
    write (ilog, '(5X,A,F7.3)')  "cutoffCN = ", this%cutoffCN

    write (ilog, *)

  endsubroutine dispdftd3_bind_to


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine dispdftd3_energy_and_forces(this, p, nl, epot, f, wpot, mask, &
       epot_per_at, wpot_per_at, ierror)
    implicit none

    type(dispdftd3_t),  intent(inout) :: this
    type(particles_t),  intent(in)    :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(inout) :: epot
    real(DP),           intent(inout) :: f(3, p%maxnatloc)  !< forces
    real(DP),           intent(inout) :: wpot(3, 3)
    integer,  optional, intent(in)    :: mask(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_at(p%maxnatloc)
#ifdef LAMMPS
    real(DP), optional, intent(inout) :: wpot_per_at(6, p%maxnatloc)
#else
    real(DP), optional, intent(inout) :: wpot_per_at(3, 3, p%maxnatloc)
#endif
    integer,  optional, intent(inout) :: ierror

    ! ---

    real(DP)  :: edisp, vbox
    real(DP)  :: unit_energy, unit_forces, unit_virial

    integer   :: elem(p%nat)
    real(DP)  :: stress(3, 3), latvecs(3, 3)

    real(DP), allocatable :: coords(:,:), grads(:,:)

    ! ---

    call timer_start("dftd3_energy_and_forces")

    allocate(coords(3, p%nat))
    allocate(grads(3, p%nat))

#ifdef HAVE_DFTD3
    elem(1:p%nat) = p%Z(1:p%nat)
    coords(1:3, 1:p%nat) = POS3(p, 1:p%nat) * length_to_Bohr
    latvecs(1:3, 1:3) = p%Abox(1:3, 1:3) * length_to_Bohr

    vbox = volume(p)

#if !defined(PYTHON) && !defined(LAMMPS)
    ! system_of_units only exists for the standalone code, LAMMPS and Python are eV/A by default
    if (system_of_units == eV_A .or. system_of_units == eV_A_fs) then
#endif

       ! Energy: a.u. -> eV
       unit_energy = Hartree

       ! Force: a.u. -> eV/A
       unit_forces = Hartree * Bohr

       ! Pressure: a.u. -> eV/A^3*
       unit_virial = Hartree / Bohr**3 * vbox

#if !defined(PYTHON) && !defined(LAMMPS)
    else 

       unit_energy = 1.0_DP
       unit_forces = 1.0_DP
       unit_virial = vbox

    endif  
#endif

    !> 
    !> For periodic system
    !> 

    call dftd3_pbc_dispersion(this%calculator, coords, elem, latvecs, edisp, grads, stress)

    epot = epot + edisp * unit_energy

    f = f - grads * unit_forces

    wpot = wpot - stress * unit_virial

#endif

    deallocate(coords)
    deallocate(grads)

    call timer_stop("dftd3_energy_and_forces")


  endsubroutine dispdftd3_energy_and_forces

  subroutine dispdftd3_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(dispdftd3_t), target, intent(inout) :: this
    type(c_ptr),               intent(in)    :: cfg
    type(c_ptr),               intent(out)   :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("DFTD3"), &
         CSTR("The DFT-D3 dispersion potential"))

    call ptrdict_register_boolean_property(m, c_loc(this%BeckeJohnson), &
         CSTR("BeckeJohnson"), &
         CSTR("Enable Becke-Johnson damping in DFT-D3."))
    call ptrdict_register_boolean_property(m, c_loc(this%threebody), &
         CSTR("threebody"), &
         CSTR("Enable three-bondy therm in DFT-D3."))

    call ptrdict_register_real_property(m, c_loc(this%a1), "a1" // char(0), &
         "Becke-Johnson-damping parameter." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%a2), "a2" // char(0), &
         "Becke-Johnson-damping parameter." // char(0))

    call ptrdict_register_real_property(m, c_loc(this%sr6), "sr6" // char(0), &
         "Zero-damping parameter." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%sr8), "sr8" // char(0), &
         "Zero-damping parameter." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%alpha6), "alpha6" // char(0), &
         "Zero-damping parameter." // char(0))

    call ptrdict_register_real_property(m, c_loc(this%s6), "s6" // char(0), &
         "Functional-dependent coefficient." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%s8), "s8" // char(0), &
         "Functional-dependent coefficient." // char(0))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), "cutoff" // char(0), &
         "Potential cutoff: If smaller than zero, the cutoff is set to a default value (90 Ang)." // char(0))
    call ptrdict_register_real_property(m, c_loc(this%cutoffCN), "cutoffCN" // char(0), &
         "Potential cutoff for three-body term: If smaller than zero, the cutoff is set to a default value (40 Ang)." // char(0))

  endsubroutine dispdftd3_register

endmodule dispdftd3
