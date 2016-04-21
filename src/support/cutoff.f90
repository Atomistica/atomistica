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

!>
!! Different cutoff functions
!<

module cutoff
  use system_module
  use units_module

  implicit none

  private

  public :: trig_on, trig_off, exp_cutoff

  public :: trig_on_t
  type trig_on_t
     real(DP) :: r1, r2, fac
  endtype trig_on_t

  public :: trig_off_t
  type trig_off_t
     real(DP) :: r1, r2, fac
  endtype trig_off_t

  public :: exp_cutoff_t
  type exp_cutoff_t
     real(DP) :: r1, r2, fac1, fac2, c, d, off
  endtype exp_cutoff_t

  public :: init
  interface init
     module procedure trig_on_init, trig_off_init, exp_cutoff_init
  endinterface

  public :: fc
  interface fc
     module procedure trig_on_f, trig_off_f, exp_cutoff_f
  endinterface fc

contains

  !>
  !! Initialize trigonometric cutoff
  !<
  elemental subroutine trig_on_init(this, r1, r2)
    implicit none

    type(trig_on_t), intent(out) :: this
    real(DP),        intent(in)  :: r1, r2

    ! ---

    this%r1  = r1
    this%r2  = r2
    this%fac = PI/( this%r2 - this%r1 )

  endsubroutine trig_on_init


  !>
  !! This if f(x)=0.5(1+cos(pi*x)).
  !! Function is differentiable once.
  !<
  subroutine trig_on_f(this, r, val, dval)
    implicit none

    type(trig_on_t), intent(in)  :: this
    real(DP),        intent(in)  :: r
    real(DP),        intent(out) :: val, dval

    ! ---

    real(DP) :: x

    ! ---

    if (r <= this%r1) then
       val  = 0.0_DP
       dval = 0.0_DP
    else if (r >= this%r2) then
       val  = 1.0_DP
       dval = 0.0_DP
    else
       x    = this%fac*( r - this%r1 )
       val  = 0.5_DP *  ( 1.0_DP - cos( x ) )
       dval = 0.5_DP * this%fac * sin( x )
    endif

  endsubroutine trig_on_f


  !>
  !! This if f(x)=0.5(1+cos(pi*x)).
  !! Function is differentiable once.
  !<
  subroutine trig_on(r1, r2, r, val, dval)
    implicit none

    real(DP), intent(in)  :: r1, r2, r
    real(DP), intent(out) :: val, dval

    ! ---

    type(trig_on_t) :: this

    ! ---

    if (r <= r1) then
       val  = 0.0_DP
       dval = 0.0_DP
    else if (r >= r2) then
       val  = 1.0_DP
       dval = 0.0_DP
    else
       call init(this, r1, r2)
       call fc(this, r, val, dval)
    endif

  endsubroutine trig_on


  !>
  !! Initialize trigonometric cutoff
  !<
  elemental subroutine trig_off_init(this, r1, r2)
    implicit none

    type(trig_off_t), intent(out) :: this
    real(DP),         intent(in)  :: r1, r2

    ! ---

    this%r1  = r1
    this%r2  = r2
    this%fac = PI/( this%r2 - this%r1 )

  endsubroutine trig_off_init


  !>
  !! This if f(x)=0.5(1+cos(pi*x)).
  !! Function is differentiable once.
  !<
  subroutine trig_off_f(this, r, val, dval)
    implicit none

    type(trig_off_t), intent(in)  :: this
    real(DP),         intent(in)  :: r
    real(DP),         intent(out) :: val, dval

    ! ---

    real(DP) :: x

    ! ---

    if (r <= this%r1) then
       val  = 1.0_DP
       dval = 0.0_DP
    else if (r >= this%r2) then
       val  = 0.0_DP
       dval = 0.0_DP
    else
       x    = this%fac*( r - this%r1 )
       val  = 0.5_DP *  ( 1.0_DP + cos( x ) )
       dval = -0.5_DP * this%fac * sin( x )
    endif

  endsubroutine trig_off_f


  !>
  !! This if f(x)=0.5(1+cos(pi*x)).
  !! Function is differentiable once.
  !<
  subroutine trig_off(r1, r2, r, val, dval)
    implicit none

    real(DP), intent(in)  :: r1, r2, r
    real(DP), intent(out) :: val, dval

    ! ---

    type(trig_off_t) :: this

    ! ---

    if (r <= r1) then
       val  = 1.0_DP
       dval = 0.0_DP
    else if (r >= r2) then
       val  = 0.0_DP
       dval = 0.0_DP
    else
       call init(this, r1, r2)
       call fc(this, r, val, dval)
    endif

  endsubroutine trig_off


  !>
  !! Initialize exponential cutoff
  !<
  elemental subroutine exp_cutoff_init(this, r1, r2)
    implicit none

    type(exp_cutoff_t), intent(out) :: this
    real(DP),           intent(in)  :: r1, r2

    ! ---

    real(DP) :: val1, dval1, ddval1

    ! ---

    this%r1   = r1
    this%r2   = r2
    this%fac1 = 1.0_DP/( this%r2 - this%r1 )
    val1      = exp(-8.0_DP)
    dval1     = -24*val1
    ddval1    = -48*val1-24*dval1
    this%c    = (-3*dval1+ddval1)/3
    this%d    = (2*dval1-ddval1)/4
    this%fac2 = 1.0_DP/(1-val1-this%c-this%d)
    this%off  = val1+this%c+this%d

  endsubroutine exp_cutoff_init


  !>
  !! This is f(x)=exp(-8*x**3), but corrected such that function, first
  !! and second derivative go to zero at x=1.
  !! Function is differentiable twice.
  !<
  subroutine exp_cutoff_f(this, r, val, dval)
    implicit none

    type(exp_cutoff_t), intent(in)  :: this
    real(DP),           intent(in)  :: r
    real(DP),           intent(out) :: val, dval

    ! ---

    real(DP) :: x, x2

    ! ---

    if (r <= this%r1) then
       val  = 1.0_DP
       dval = 0.0_DP
    else if (r >= this%r2) then
       val  = 0.0_DP
       dval = 0.0_DP
    else
       x      = this%fac1*( r-this%r1 )
       x2     = x*x
       val    = exp(-8*x*x2)
       dval   = -24*x2*val
       ! The following two lines are the correction that forces the cutoff to
       ! zero.
       dval   = this%fac1*this%fac2*(dval+3*this%c*x2+4*this%d*x*x2)
       val    = this%fac2*(val+this%c*x*x2+this%d*x2*x2-this%off)
    endif

  endsubroutine exp_cutoff_f


  !>
  !! This is f(x)=exp(-8*x**3), but corrected such that function, first
  !! and second derivative go to zero at x=1.
  !! Function is differentiable twice.
  !<
  subroutine exp_cutoff(r1, r2, r, val, dval)
    implicit none

    real(DP), intent(in)  :: r1, r2, r
    real(DP), intent(out) :: val, dval

    ! ---

    type(exp_cutoff_t) :: this

    ! ---

    if (r <= r1) then
       val  = 1.0_DP
       dval = 0.0_DP
    else if (r >= r2) then
       val  = 0.0_DP
       dval = 0.0_DP
    else
       call init(this, r1, r2)
       call fc(this, r, val, dval)
    endif

  endsubroutine exp_cutoff

endmodule cutoff
