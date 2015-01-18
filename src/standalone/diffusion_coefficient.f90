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
!   classtype:diffusion_coefficient_t classname:DiffusionCoefficient
!   interface:callables
! @endmeta

!>
!! Computation of diffusion coefficients
!!
!! Computation of diffusion coefficients
!<

#include "macros.inc"
#include "filter.inc"

module diffusion_coefficient
  use supplib

  use particles
  use neighbors
  use dynamics
  use filter

  implicit none

  private

  character(MAX_NAME_STR), parameter  :: NORM_STR       = "diffusion_coefficient_norm"
  character(MAX_NAME_STR), parameter  :: ABS_DR_SQ_STR  = "diffusion_coefficient_abs_dr_sq"
  character(MAX_NAME_STR), parameter  :: DR_SQ_STR      = "diffusion_coefficient_dr_sq"

  character(MAX_NAME_STR), parameter  :: R0_STR         = "diffusion_coefficient_r0"

  character(MAX_NAME_STR), parameter  :: CUR_DR_SQ_STR  = "diffusion_coefficient_cur_abs_dr"
  character(MAX_NAME_STR), parameter  :: DC_STR         = "diffusion_coefficient"
  character(MAX_NAME_STR), parameter  :: DC3_STR        = "diffusion_coefficient_cart"

  public :: diffusion_coefficient_t
  type diffusion_coefficient_t

     !
     ! From input file
     !

     character(MAX_EL_STR)  :: element
     integer                :: el

     real(DP)               :: freq

     real(DP)               :: tau

     !
     ! Helper stuff
     !

     real(DP)               :: ti, tot_ti

     integer                :: un

     real(DP), pointer      :: r0(:, :)

     !
     ! Running averages
     !

     real(DP), pointer      :: norm => NULL()
     real(DP), pointer      :: abs_dr_sq => NULL()
     real(DP), pointer      :: dr_sq(:) => NULL()

     !
     ! Per particle diffusion coefficient
     !

     real(DP), pointer      :: cur_dr_sq(:, :) => NULL()

     real(DP), pointer      :: abs_Dc(:) => NULL()
     real(DP), pointer      :: Dc(:, :) => NULL()

  endtype diffusion_coefficient_t

  public :: register_data
  interface register_data
     module procedure diffusion_coefficient_register_data
  endinterface

  public :: init
  interface init
     module procedure diffusion_coefficient_init
  endinterface

  public :: del
  interface del
     module procedure diffusion_coefficient_del
  endinterface

  public :: invoke
  interface invoke
     module procedure diffusion_coefficient_invoke
  endinterface

  public :: register
  interface register
    module procedure diffusion_coefficient_register
  endinterface

contains


  !**********************************************************************
  ! Initialize a diffusion_coefficient object
  !**********************************************************************
  subroutine diffusion_coefficient_register_data(this, p)
    implicit none

    type(diffusion_coefficient_t), intent(inout)  :: this
    type(particles_t), intent(inout)              :: p

    ! ---

    call add_real_attr(p%data, NORM_STR)
    call add_real_attr(p%data, ABS_DR_SQ_STR)
    call add_real3_attr(p%data, DR_SQ_STR)

    call add_real3(p%data, R0_STR, F_RESTART + F_VERBOSE_ONLY)
    call add_real3(p%data, CUR_DR_SQ_STR, F_RESTART + F_VERBOSE_ONLY)
    call add_real (p%data, DC_STR)
    call add_real3(p%data, DC3_STR)

  endsubroutine diffusion_coefficient_register_data
  

  !**********************************************************************
  ! Initialize a diffusion_coefficient object
  !**********************************************************************
  subroutine diffusion_coefficient_init(this)
    implicit none

    type(diffusion_coefficient_t), intent(inout)  :: this

    ! ---

    write (ilog, '(A)')  "- diffusion_coefficient_init -"

    write (ilog, '(5X,A,A)')       "element   = ", this%element

    this%ti         = 0.0_DP
    this%tot_ti     = 0.0_DP

    this%un  = fopen("diffusion_coefficient.out", F_WRITE)
    write (this%un, '(A)')  "# 1:ti 2:<D> 3:<Dx> 4:<Dy> 5:<Dz> 6:<dr^2> 7:<dx^2> 8:<dy^2> 9:<dz^2> 10:D 11:Dx 12:Dy 13:Dz 14:dr^2 15:dx^2 16:dy^2 17:dz^2"

    write (ilog, *)

  endsubroutine diffusion_coefficient_init


  !**********************************************************************
  ! Delete a diffusion_coefficient object
  !**********************************************************************
  subroutine diffusion_coefficient_del(this)
    implicit none

    type(diffusion_coefficient_t), intent(inout)  :: this

    ! ---

    call fclose(this%un)

  endsubroutine diffusion_coefficient_del


  !**********************************************************************
  ! Perform the measurement
  !**********************************************************************
  subroutine diffusion_coefficient_invoke(this, dyn, nl, error)
    implicit none

    type(diffusion_coefficient_t), intent(inout)  :: this
    type(dynamics_t), intent(inout)               :: dyn
    type(neighbors_t), intent(in)                 :: nl
    integer, intent(inout), optional              :: error

    ! ---

    integer                :: i, n
    real(DP)               :: dr(3), dr_sq(3), Dc(3)
!    real(DP)               :: dx_dy, dy_dz, dz_dx, Dxy, Dyz, Dzx
    real(DP)               :: fac, dfac, dfac2

    real(DP)               :: old_norm

    ! ---

    call timer_start("diffusion_coefficient_invoke")

    if (.not. associated(this%norm)) then
       this%el = filter_from_string(this%element, dyn%p)

       call attr_by_name(dyn%p%data, NORM_STR, this%norm)
       call attr_by_name(dyn%p%data, ABS_DR_SQ_STR, this%abs_dr_sq)
       call attr_by_name(dyn%p%data, DR_SQ_STR, this%dr_sq)

       call ptr_by_name(dyn%p%data, R0_STR, this%r0)
       call ptr_by_name(dyn%p%data, CUR_DR_SQ_STR, this%cur_dr_sq)
       call ptr_by_name(dyn%p%data, DC_STR, this%abs_Dc)
       call ptr_by_name(dyn%p%data, DC3_STR, this%Dc)

       do i = 1, dyn%p%nat
          if (IS_EL(this%el, dyn%p, i)) then
             VEC3(this%r0, i) = PCN3(dyn%p, i)
          endif
       enddo
    endif

    this%ti               = this%ti + dyn%dt
    this%tot_ti           = this%tot_ti + dyn%dt

    old_norm              = this%norm

    fac                   = exp(-dyn%dt/this%tau)
    if (old_norm > 0.0_DP) then
       dfac               = ( fac/this%norm - 1.0_DP/old_norm )/(2*dyn%dt)
       dfac2              = 1.0_DP/(2*this%norm)
    else
       dfac               = 0.0_DP
       dfac2              = 0.0_DP
    endif

    this%norm             = fac*this%norm + dyn%dt

    dr_sq(:)              = 0.0_DP
!    dx_dy                 = 0.0_DP
!    dy_dz                 = 0.0_DP
!    dz_dx                 = 0.0_DP
    
    n = 0
    !$omp  parallel do default(none) &
    !$omp& shared(dfac, dfac2, dyn, fac, this) &
    !$omp& private(Dc, dr) &
!!  !$omp& private(Dxy, Dyz, Dzx, dx_dy, dy_dz, dz_dx) &
    !$omp& reduction(+:dr_sq) reduction(+:n)
    do i = 1, dyn%p%nat
       if (IS_EL(this%el, dyn%p, i)) then
          n                         = n+1
          dr(:)                     = PCN3(dyn%p, i) - VEC3(this%r0, i)
!          dx_dy                     = dx_dy + dr(1)*dr(2)
!          dy_dz                     = dy_dz + dr(2)*dr(3)
!          dz_dx                     = dz_dx + dr(3)*dr(1)
          dr(:)                     = dr(:)*dr(:)
          dr_sq(:)                  = dr_sq(:) + dr(:)

          Dc(:)                     = dfac*VEC3(this%cur_dr_sq, i) + dfac2*dr(:)
!          Dxy                       = dfac*this%cur_dx_dy(gi) + dfac2*dx_dy
!          Dyz                       = dfac*this%cur_dy_dz(gi) + dfac2*dy_dz
!          Dzx                       = dfac*this%cur_dz_dx(gi) + dfac2*dz_dx

          VEC3(this%cur_dr_sq, i)   = fac*VEC3(this%cur_dr_sq, i) + dyn%dt*dr(:)
!          this%cur_dx_dy(gi)        = fac*this%cur_dx_dy(gi) + dt*dx_dy
!          this%cur_dy_dz(gi)        = fac*this%cur_dy_dz(gi) + dt*dy_dz
!          this%cur_dz_dx(gi)        = fac*this%cur_dz_dx(gi) + dt*dz_dx

          VEC3(this%Dc, i)          = Dc(:)
!          this%Dxy(i)               = Dxy
!          this%Dyz(i)               = Dyz
!          this%Dzx(i)               = Dzx
          this%abs_Dc(i)            = ( Dc(1) + Dc(2) + Dc(3) )/3
       endif
    enddo

    dr_sq(:)        = dr_sq(:)/n

    if (this%ti > this%freq) then

       write (this%un, '(18ES20.10)')  &
            dyn%ti, &
            ( dfac*this%abs_dr_sq + dfac2*(dr_sq(1)+dr_sq(2)+dr_sq(3)) )/3, &
            dfac*this%dr_sq(:) + dfac2*dr_sq(:), &
            this%abs_dr_sq/this%norm, &
            this%dr_sq(:)/this%norm, &
            (dr_sq(1)+dr_sq(2)+dr_sq(3))/(6*this%tot_ti), &
            dr_sq(:)/(2*this%tot_ti), &
            (dr_sq(1)+dr_sq(2)+dr_sq(3)), &
            dr_sq(:)

       this%ti = 0.0_DP

    endif

    this%abs_dr_sq  = fac*this%abs_dr_sq + dyn%dt*(dr_sq(1)+dr_sq(2)+dr_sq(3))
    this%dr_sq(:)   = fac*this%dr_sq(:)  + dyn%dt*dr_sq(:)

    call timer_stop("diffusion_coefficient_invoke")

  endsubroutine diffusion_coefficient_invoke


  subroutine diffusion_coefficient_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(diffusion_coefficient_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)                       :: cfg
    type(c_ptr), target, intent(inout)                    :: m

    ! ---

    this%element   = "*"

    this%freq      = -1.0
    this%tau       = 100.0

    m = ptrdict_register_section(cfg, CSTR("DiffusionCoefficient"), &
         CSTR("Compute diffusion coefficient by tracing particle motion, i.e. from D = (<x^2> - <x0^2>)/2t."))

    call ptrdict_register_string_property(m, c_loc(this%element), MAX_EL_STR, &
         CSTR("element"), &
         CSTR("Element for which to compute diffusion."))

    call ptrdict_register_real_property(m, c_loc(this%freq), CSTR("freq"), &
         CSTR("Output frequency."))

    call ptrdict_register_real_property(m, c_loc(this%tau), CSTR("tau"), &
         CSTR("Time constant for running averages."))

  endsubroutine diffusion_coefficient_register

endmodule diffusion_coefficient
