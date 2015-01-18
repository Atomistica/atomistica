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
!   classtype:confinement_t classname:Confinement interface:potentials
! @endmeta

!>
!! A Lennard-Jones like confinement potential
!! The particles will be confined within a slab in x-y direction
!! between z1 and z2.
!<

#include "macros.inc"
#include "filter.inc"

module confinement
  use libAtoms_module

  use io
  use logging
  use timer

  use particles
  use neighbors
  use filter

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  integer, parameter  :: n_dims       = 3
  integer, parameter  :: len_dim_str  = 15
  integer, parameter  :: ALL_DIMS     = 0

  ! This is need for xlf
  character(len_dim_str), parameter  :: STR_x            = CSTR("x")
  character(len_dim_str), parameter  :: STR_y            = CSTR("y")
  character(len_dim_str), parameter  :: STR_z            = CSTR("z")
  character(len_dim_str), parameter  :: dim_strs(n_dims) = &
       (/ STR_x, STR_y, STR_z /)


  real(DP), parameter  :: fac  = -1.0_DP/((2.0_DP/5)**(5.0_DP/3)-(2.0_DP/5)**(2.0_DP/3))

  public :: confinement_t
  type confinement_t

     !
     ! Elements on which to act
     !

     character(MAX_EL_STR)  :: elements  = "*"
     integer                :: els  = 0

     !
     ! Potential parameters
     !

     integer   :: dir  = 3
     integer   :: dir2
     integer   :: dir3

     real(DP)  :: epsilon  = 0.001_DP
     real(DP)  :: sigma    = 1.0_DP

     real(DP)  :: z1  = -1.0_DP
     real(DP)  :: z2  = -1.0_DP

     logical(BOOL)   :: output_force  = .false.

     !
     ! The cutoff
     !

     real(DP)  :: cutoff  = -1.0_DP

     !
     ! Shift
     !

     real(DP)  :: shift

     !
     ! Output file
     !

     integer   :: un
     
  endtype confinement_t


  public :: init
  interface init
     module procedure confinement_init
  endinterface

  public :: del
  interface del
     module procedure confinement_del
  endinterface

  public :: energy_and_forces
  interface energy_and_forces
     module procedure confinement_energy_and_forces
  endinterface

  public :: register
  interface register
    module procedure confinement_register
  endinterface

contains

  !**********************************************************************
  ! Initialize the confinement module
  !**********************************************************************
  subroutine confinement_init(this)
    implicit none

    type(confinement_t), intent(inout)  :: this

    ! ---

    write (ilog, '(A)')  "- confinement_init -"

    write (ilog, '(5X,A)')  trim(this%elements)

    this%dir   = mod(this%dir,   3)+1
    this%dir2  = mod(this%dir+1, 3)+1
    this%dir3  = mod(this%dir+2, 3)+1

    write (ilog, '(5X,A,I1)')      "dir      = ", this%dir
    write (ilog, '(5X,A,F20.10)')  "z1       = ", this%z1
    write (ilog, '(5X,A,F20.10)')  "z2       = ", this%z2
    write (ilog, '(5X,A,F7.3)')    "epsilon  = ", this%epsilon
    write (ilog, '(5X,A,F7.3)')    "sigma    = ", this%sigma

    if (this%cutoff < 0.0_DP) then
!       this%cutoff = (2.0d0**(1.0d0/6))*this%sigma
       this%cutoff  = (2.5_DP**(1.0_DP/6))*this%sigma
    endif

    write (ilog, '(5X,A,F7.3)')  "cutoff   = ", this%cutoff

    this%shift = (this%sigma/this%cutoff)**10-(this%sigma/this%cutoff)**4

#ifdef _MP
    if (mod_communicator%mpi%my_proc == ROOT) then
#endif

    if (this%output_force) then
       this%un  = fopen("confinement_pressure.out", F_WRITE)
    endif

#ifdef _MP
    endif
#endif

    this%els  = 0

    write (ilog, *)
    
  endsubroutine confinement_init


  !**********************************************************************
  ! Delete the confinement module
  !**********************************************************************
  subroutine confinement_del(this)
    implicit none

    type(confinement_t), intent(inout)  :: this

    ! ---

#ifdef _MP
    if (mod_communicator%mpi%my_proc == ROOT) then
#endif

    if (this%output_force) then
       call fclose(this%un)
    endif

#ifdef _MP
    endif
#endif

  endsubroutine confinement_del


  !**********************************************************************
  ! Compute the force
  !**********************************************************************
  subroutine confinement_energy_and_forces(this, p, nl, epot, for, wpot, epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)
    implicit none

    type(confinement_t), intent(inout)  :: this
    type(particles_t), intent(inout)    :: p
    type(neighbors_t), intent(in)       :: nl
    real(DP), intent(inout)             :: epot
    real(DP), intent(inout)             :: for(3, p%maxnatloc)
    real(DP), intent(inout)             :: wpot(3, 3)
    real(DP), intent(inout), optional   :: epot_per_at(p%maxnatloc)
    real(DP), intent(inout), optional   :: epot_per_bond(nl%neighbors_size)
    real(DP), intent(inout), optional   :: f_per_bond(3, nl%neighbors_size)
    real(DP), intent(inout), optional   :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), intent(inout), optional   :: wpot_per_bond(3, 3, nl%neighbors_size)
    integer, intent(inout), optional    :: ierror

    ! ---

    integer   :: i, d
    real(DP)  :: r(3), dr, s_r, f, e, f_tot

    ! ---

    call timer_start("confinement_force")

    if (this%els == 0) then
       this%els  = filter_from_string(this%elements, p)
    endif

    if (this%z1 < 0.0_DP) then
       this%z1  = 0.0_DP
    endif

    if (this%z2 < 0.0_DP) then
       this%z2  = p%Abox(this%dir, this%dir)
    endif


    d      = this%dir
    e      = 0.0_DP
    f_tot  = 0.0_DP

    !$omp  parallel do default(none) &
    !$omp& private(dr, f, r, s_r) &
    !$omp& firstprivate(d) &
    !$omp& shared(for, p, this) &
    !$omp& reduction(+:e) reduction(+:f_tot)
    do i = 1, p%natloc

       if (p%g(i) > 0 .and. IS_EL(this%els, p, i)) then

          r  = POS3(p, i)
       
          if (r(d) >= this%z1 .and. r(d) < this%z1+this%cutoff) then
       
             dr   = r(d)-this%z1
             s_r  = this%sigma/dr
          
             e    = e + fac*this%epsilon*(s_r**10 - s_r**4 - this%shift)
             f    = fac*this%epsilon*(10*s_r**10 - 4*s_r**4)/dr
          
             VEC(for, i, d)  = VEC(for, i, d) + f

             f_tot  = f_tot + f
       
          else if (r(d) > this%z2-this%cutoff .and. r(d) <= this%z2) then

             dr   = this%z2-r(d)
             s_r  = this%sigma/dr
          
             e    = e + fac*this%epsilon*(s_r**10 - s_r**4 - this%shift)
             f    = fac*this%epsilon*(10*s_r**10 - 4*s_r**4)/dr
          
             VEC(for, i, d)  = VEC(for, i, d) - f

             f_tot  = f_tot - f

          endif
       
       endif

    enddo

    epot = epot + e

#ifdef _MP
    if (this%output_force) then
       call sum_in_place(mod_communicator%mpi, f_tot)
    endif
  
    if (mod_communicator%mpi%my_proc == ROOT) then
#endif

    if (this%output_force) then
       write (this%un, '(ES20.10)')  f_tot/(p%Abox(this%dir2, this%dir2)*p%Abox(this%dir3, this%dir3))
    endif

#ifdef _MP
    endif
#endif

    call timer_stop("confinement_force")

  endsubroutine confinement_energy_and_forces


  subroutine confinement_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(confinement_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                :: cfg
    type(c_ptr), intent(out)               :: m

    ! ---
 
    m = ptrdict_register_section(cfg, CSTR("Confinement"), &
         CSTR("A Lennard-Jones-like confinement potential. Particles are confined in the x-y plane."))

    call ptrdict_register_string_property(m, c_loc(this%elements), MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("List of elements on which this potential should act."))

    call ptrdict_register_enum_property(m, c_loc(this%dir), &
         n_dims-1, len_dim_str, dim_strs(:), &
         CSTR("d"), &
         CSTR("Direction in which to confinge the simulation: 'x', 'y', 'z'"))

    call ptrdict_register_real_property(m, c_loc(this%epsilon), &
         CSTR("epsilon"), &
         CSTR("Interaction energy."))
    call ptrdict_register_real_property(m, c_loc(this%sigma), CSTR("sigma"), &
         CSTR("Interaction diameter."))

    call ptrdict_register_real_property(m, c_loc(this%z1), CSTR("z1"), &
         CSTR("Lower bound (z-direction)."))
    call ptrdict_register_real_property(m, c_loc(this%z2), CSTR("z2"), &
         CSTR("Upper bound (z-direction)."))

    call ptrdict_register_real_property(m, c_loc(this%cutoff), CSTR("cutoff"), &
         CSTR("Potential cutoff: If smaller than zero, the cutoff is set such that the potential is only repulsive."))

    call ptrdict_register_boolean_property(m, c_loc(this%output_force), &
         CSTR("output_force"), &
         CSTR("Output pressure on the sidewalls."))

  endsubroutine confinement_register

endmodule confinement
