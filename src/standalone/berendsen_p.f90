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
!   classtype:berendsen_p_t classname:BerendsenP interface:callables
! @endmeta

!>
!! Berendsen pressure control
!!
!! Berendsen pressure control according to
!! H. J. C. Berendsen et al., J. Chem. Phys. 81, 3684 (1984).
!!
!! The pressure is controlled by multiplying box vectors and positions with
!! the factor
!! \f[
!!   s = \left( 1 + \frac{\Delta t \beta}{\tau} \left( \frac{1}{3}\textrm{tr} P - P_0 \right) \right)^\frac{1}{3}
!! \f]
!! in the hydrostatic case where $P$ is the pressure and $P_0$ is the target
!! pressure. It is also possible to control the pressure in each of the three
!! cartesian directions individually.
!<

#include "macros.inc"

module berendsen_p
  use supplib

  use particles
  use neighbors
  use dynamics
  
  implicit none

  private

  integer, parameter  :: n_dims       = 6
  integer, parameter  :: len_dim_str  = 15
  integer, parameter  :: ALL_DIMS     = 0
  integer, parameter  :: DIMS_XY      = 5
  integer, parameter  :: HYDROSTATIC  = 4

  ! This is need for xlf
  character(len_dim_str), parameter  :: STR_all          = CSTR("all")
  character(len_dim_str), parameter  :: STR_x            = CSTR("x")
  character(len_dim_str), parameter  :: STR_y            = CSTR("y")
  character(len_dim_str), parameter  :: STR_z            = CSTR("z")
  character(len_dim_str), parameter  :: STR_hydrostatic  = CSTR("hydrostatic")
  character(len_dim_str), parameter  :: STR_xy           = CSTR("xy")
  character(len_dim_str), parameter  :: dim_strs(n_dims) = &
       (/ STR_all, STR_x, STR_y, STR_z, STR_hydrostatic, STR_xy /)

  public :: berendsen_p_t
  type berendsen_p_t

     !
     ! Hydrostatic pressure components
     !

     real(DP)  :: P(3)   = (/ 0.0_DP, 0.0_DP, 0.0_DP /)

     real(DP)  :: tau(3) = (/ 100.0_DP, 100.0_DP, 100.0_DP /)

     real(DP)  :: beta   = 1.0_DP

     integer   :: d      = ALL_DIMS

     !
     ! Off diagonal pressure components (using Lees-Edward BCs)
     !

     logical(BOOL)   :: shear_stress  = .false.

     real(DP)  :: shear_sigma   = 0.0_DP
     
     integer   :: shear_d       = 1
     
     !
     ! Else
     !

     logical(BOOL)  :: log  = .false.          !< Log box size and volume to a file
     integer  :: un

  endtype berendsen_p_t


  public :: init
  interface init
     module procedure berendsen_p_init
  endinterface

  public :: set
  interface set
     module procedure berendsen_p_init
  endinterface

  public :: adjust_pressure
  interface adjust_pressure
     module procedure berendsen_p_adjust_pressure
  endinterface

  public :: del
  interface del
     module procedure berendsen_p_del
  endinterface

  public :: invoke
  interface invoke
     module procedure berendsen_p_invoke
  endinterface

  public :: register
  interface register
    module procedure berendsen_p_register
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Initialize a berendsen_p object
  !<
  subroutine berendsen_p_init(this, P, tau, beta, d, shear_sigma, shear_d, log)
    implicit none

    type(berendsen_p_t), intent(inout)   :: this
    real(DP), intent(in), optional       :: P(3)
    real(DP), intent(in), optional       :: tau(3)
    real(DP), intent(in), optional       :: beta
    integer, intent(in), optional        :: d
    integer, intent(in), optional        :: shear_sigma
    integer, intent(in), optional        :: shear_d
    logical, intent(in), optional     :: log

    ! ---

    if (present(P)) then
       this%P  = P
    endif
    if (present(tau)) then
       this%tau  = tau
    endif
    if (present(beta)) then
       this%beta  = beta
    endif
    if (present(d)) then
       this%d  = d
    endif

    if (present(shear_sigma)) then
       this%shear_stress  = .true.
       this%shear_sigma   = shear_sigma
    endif
    if (present(shear_d)) then
       this%shear_stress  = .true.
       this%shear_d       = shear_d
    endif

    if (present(log)) then
       this%log  = log
    endif

    call prlog("- berendsen_p_init -")
    call prlog("     Using Berendsen pressure control with parameters")
    call prlog("     log  = " // logical(this%log))
    call prlog("     P    = " // this%P(1) // " " // this%P(2) // " " // this%P(3))
    call prlog("     tau  = " // this%tau(1) // " " // this%tau(2) // " " // this%tau(3))
    call prlog("     beta = " // this%beta)
    call prlog("     d    = " // this%d)
    call prlog
    call prlog("     shear_stress = " // logical(this%shear_stress))
    call prlog("     shear_sigma  = " // this%shear_sigma)
    call prlog("     shear_d      = " // this%shear_d)
    call prlog

    if (this%log) then
       this%un  = fopen("berendsen_p.out", F_WRITE)
       write (this%un, '(A6,14X,4A20)')  "#01:ti", "02:cellx", "03:celly", "04:cellz", "05:V"
    endif

  endsubroutine berendsen_p_init


  !>
  !! Pressure control
  !!
  !! Adjust the current cell size to match the target pressure.
  !! To be called after the second Verlet step, and after
  !! BerendsenT
  !<
  subroutine berendsen_p_adjust_pressure(this, p, pressure, dt, ti, ierror)
    implicit none

    type(berendsen_p_t), intent(inout)  :: this
    type(particles_t), intent(inout)    :: p
    real(DP), intent(in)                :: pressure(3, 3)
    real(DP), intent(in)                :: dt
    real(DP), intent(in), optional      :: ti
    integer, intent(inout), optional    :: ierror

    ! ---

    integer   :: i
    real(DP)  :: r, s(3), t(3, 3), Abox(3, 3)

    ! ---

    call timer_start("berendsen_p_adjust_pressure")

    s  = 1.0_DP
    if (this%d == HYDROSTATIC) then
       s          = ( 1 + dt*this%beta*( tr(3, pressure)/3 - this%P(1) )/this%tau(1) )**(1.d0/3)
    else if (this%d == ALL_DIMS) then
       s(1)       = ( 1 + dt*this%beta*( pressure(1, 1) - this%P(1) )/this%tau(1) )  !**(1.d0/3)
       s(2)       = ( 1 + dt*this%beta*( pressure(2, 2) - this%P(2) )/this%tau(2) )  !**(1.d0/3)
       s(3)       = ( 1 + dt*this%beta*( pressure(3, 3) - this%P(3) )/this%tau(3) )  !**(1.d0/3)
    else if (this%d == DIMS_XY) then
       s(1)       = ( 1 + dt*this%beta*( pressure(1, 1) - this%P(1) )/this%tau(1) )  !**(1.d0/3)
       s(2)       = ( 1 + dt*this%beta*( pressure(2, 2) - this%P(2) )/this%tau(2) )  !**(1.d0/3)
    else if (this%d == 1 .or. this%d == 2 .or. this%d == 3) then
       s(this%d)  = ( 1 + dt*this%beta*( pressure(this%d, this%d) - this%P(this%d) )/this%tau(this%d) )  !**(1.d0/3)
    else
       RAISE_ERROR("BerendsenP does not support '" // dim_strs(this%d+1) // "' mode.", ierror)
    endif

    if (.not. this%shear_stress) then

       Abox  = 0.0_DP
       do i = 1, 3
#ifndef IMPLICIT_R
          POS(p, :, i)  = POS(p, :, i) * s(i)
#endif
          PNC(p, :, i)  = PNC(p, :, i) * s(i)
          Abox(i, i)    = p%Abox(i, i) * s(i)
       enddo

       call set_cell(p, Abox, error=ierror)
       PASS_ERROR(ierror)
       
    else

       r  = dt*this%beta*( pressure(this%shear_d, 3) - this%shear_sigma )/this%tau(1)

       t                   = 0.0_DP
       t(1, 1)             = s(1)
       t(2, 2)             = s(2)
       t(3, 3)             = s(3)
       t(this%shear_d, 3)  = r

       do i = 1, p%nat
#ifndef IMPLICIT_R
          POS3(p, i)  = matmul(t, POS3(p, i))
#endif
          PNC3(p, i)  = matmul(t, PNC3(p, i))
       enddo

       Abox        = 0.0_DP
       Abox(1, 1)  = s(1)*p%Abox(1, 1)
       Abox(2, 2)  = s(2)*p%Abox(2, 2)
       Abox(3, 3)  = s(3)*p%Abox(3, 3)

       p%shear_dx(this%shear_d)  = r*p%Abox(3, 3) + s(1)*p%shear_dx(this%shear_d)
       call set_cell(p, Abox, error=ierror)
       PASS_ERROR(ierror)

    endif

    if (this%log) then
       if (present(ti)) then
          write (this%un, '(5ES20.10)')  ti, p%Abox(1, 1), p%Abox(2, 2), p%Abox(3, 3), P%Abox(1, 1) * P%Abox(2, 2) * P%Abox(3, 3)
       else
          write (this%un, '(4ES20.10)')  p%Abox(1, 1), p%Abox(2, 2), p%Abox(3, 3), P%Abox(1, 1) * P%Abox(2, 2) * P%Abox(3, 3)
       endif
    endif

    call timer_stop("berendsen_p_adjust_pressure")

  endsubroutine berendsen_p_adjust_pressure

  !>
  !! Destructor
  !!
  !! Destroy a berendsen_p object
  !<
  subroutine berendsen_p_del(this)
    implicit none

    type(berendsen_p_t), intent(inout)   :: this

    ! ---

    if (this%log) then
       call fclose(this%un)
    endif

    ! ---

  end subroutine berendsen_p_del


  !>
  !! Adjuste the pressure
  !!
  !! Adjuste the pressure
  !<
  subroutine berendsen_p_invoke(this, dyn, nl, ierror)
    implicit none

    type(berendsen_p_t), intent(inout)  :: this
    type(dynamics_t), intent(inout)     :: dyn
    type(neighbors_t), intent(in)       :: nl
    integer, intent(inout), optional    :: ierror

    ! ---

    call adjust_pressure(this, dyn%p, dyn%pressure, dyn%dt, dyn%ti, ierror)
    PASS_ERROR(ierror)

  endsubroutine berendsen_p_invoke


  subroutine berendsen_p_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(berendsen_p_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)             :: cfg
    type(c_ptr), intent(out)            :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("BerendsenP"), &
         CSTR("Berendsen barostat (H.J.C. Berendsen, J.P.M. Postma, W.F. van Gunsteren, A. DiNola, J.R. Haak, J. Chem. Phys. 81, 3684 (1984)."))

    call ptrdict_register_point_property(m, c_loc(this%P(1)), CSTR("P"), &
         CSTR("Target pressure (in x-, y- and z-direction)."))

    call ptrdict_register_point_property(m, c_loc(this%tau(1)), CSTR("tau"), &
         CSTR("Pressure coupling time constants (in x-, y- and z-direction)."))

    call ptrdict_register_real_property(m, c_loc(this%beta), CSTR("beta"), &
         CSTR("Isothermal compressibility."))

    call ptrdict_register_enum_property(m, c_loc(this%d), &
         n_dims, len_dim_str, dim_strs, &
         CSTR("d"), &
         CSTR("Dimension for pressure equilization: 'x', 'y', 'z', 'xy', 'all' or 'hydrostatic'"))

    call ptrdict_register_boolean_property(m, c_loc(this%shear_stress), CSTR("shear_stress"), &
         CSTR("In addition to hydrostatic pressure apply a shear stress."))

    call ptrdict_register_real_property(m, c_loc(this%shear_sigma), CSTR("shear_sigma"), &
         CSTR("Target shear stress."))

    call ptrdict_register_enum_property(m, c_loc(this%shear_d), &
         3, len_dim_str, dim_strs, &
         CSTR("shear_d"), &
         CSTR("Shearing direction: 'x' or 'y'"))

    call ptrdict_register_boolean_property(m, c_loc(this%log), CSTR("log"), &
         CSTR("Log cell size and volume of cell to 'berendsen_p.out'."))

  endsubroutine berendsen_p_register

endmodule berendsen_p
