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
!   dependencies:verlet_support.f90
!   classtype:fire_t classname:FIRE interface:integrators
! @endmeta

#include "macros.inc"

!>
!! The FIRE algorithm
!!
!! The Fast Interial Relaxation Engine, see
!! Bitzek, Koskinen, Gumbsch, Moseler paper
!<
module fire
  use supplib

  use particles

  use verlet_support

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  public :: fire_t
  type fire_t

     !
     ! Configuration
     !

     integer   :: minsteps  = 10
     real(DP)  :: incfac    = 1.2_DP
     real(DP)  :: decfac    = 0.5_DP
     real(DP)  :: mix_in    = 0.1_DP
     real(DP)  :: mixdec    = 0.99_DP
     real(DP)  :: max_dt    = 0.01_DP

     real(DP)  :: fmax      = 0.01_DP

     !
     ! Local variables
     !

     real(DP)  :: mix
     integer   :: cut
     integer   :: cuts

     !
     ! Verbose mode
     !

     logical(C_BOOL) :: log

  endtype fire_t


  public :: init
  interface init
     module procedure fire_init
  endinterface

  public :: del
  interface del
     module procedure fire_del
  endinterface

  public :: step1
  interface step1
     module procedure fire_step1
  endinterface

  public :: step2
  interface step2
     module procedure fire_step2
  endinterface

  public :: register
  interface register
    module procedure fire_register
  endinterface

contains

  !**********************************************************************
  ! Constructor
  !**********************************************************************
  subroutine fire_init(this, p, minsteps, incfac, decfac, mix_in, mixdec, max_dt)
    implicit none

    type(fire_t), intent(inout)     :: this
    type(particles_t), intent(in)   :: p
    integer, intent(in), optional   :: minsteps
    real(DP), intent(in), optional  :: incfac
    real(DP), intent(in), optional  :: decfac
    real(DP), intent(in), optional  :: mix_in
    real(DP), intent(in), optional  :: mixdec
    real(DP), intent(in), optional  :: max_dt

    ! ---

    ASSIGN_PROPERTY(minsteps)
    ASSIGN_PROPERTY(incfac)
    ASSIGN_PROPERTY(decfac)
    ASSIGN_PROPERTY(mix_in)
    ASSIGN_PROPERTY(mixdec)
    ASSIGN_PROPERTY(max_dt)

  endsubroutine fire_init


  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine fire_del(this)
    implicit none

    type(fire_t), intent(inout)  :: this

    ! ---

  endsubroutine fire_del


  !**********************************************************************
  ! Position update and velocity estimation
  !**********************************************************************
  subroutine fire_step1(this, p, v, f, dt, max_dt, max_dr, max_dr_sq)
    implicit none

    type(fire_t), intent(inout)       :: this
    type(particles_t), intent(inout)  :: p
    real(DP), intent(inout)           :: v(3, p%maxnatloc)
    real(DP), intent(in)              :: f(3, p%maxnatloc)
    real(DP), intent(inout)           :: dt
    real(DP), intent(in), optional    :: max_dt
    real(DP), intent(in), optional    :: max_dr
    real(DP), intent(inout), optional :: max_dr_sq

    ! ---

    real(DP)           :: d2t

    integer            :: i

    real(DP)           :: dr(3), l_max_dr_sq

    ! ---

    call timer_start("fire_step1")

    !
    ! Adaptive time stepping
    !

    call timestep(p, v, f, dt, max_dt, max_dr)
        
    d2t = dt**2

    !
    ! Integrate
    !

    l_max_dr_sq  = 0.0_DP

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, v) &
    !$omp& firstprivate(dt, d2t) &
    !$omp& private(dr) &
    !$omp& reduction(max:l_max_dr_sq)
    do i = 1, p%natloc

       if (p%g(i) > 0) then

          dr               = VEC3(v, i) * dt + 0.5_DP * VEC3(f, i) * d2t
#ifndef IMPLICIT_R
          POS3(p, i)       = POS3(p, i) + dr
#endif
          PNC3(p, i)       = PNC3(p, i) + dr
          PCN3(p, i)       = PCN3(p, i) + dr
          VEC3(v, i)       = VEC3(v, i) + 0.5_DP * VEC3(f, i) * dt

          l_max_dr_sq      = max(l_max_dr_sq, dot_product(dr, dr))

       endif

    enddo

    !
    ! Maximum particle displacement
    !

    p%accum_max_dr  = p%accum_max_dr + sqrt(l_max_dr_sq)

    if (present(max_dr_sq)) then
       max_dr_sq  = max(max_dr_sq, l_max_dr_sq)
    endif

    call I_changed_positions(p)

    call timer_stop("fire_step1")

  endsubroutine fire_step1


  !**********************************************************************
  ! Velocity correction
  !**********************************************************************
  subroutine fire_step2(this, p, v, f, dt)
    implicit none

    type(fire_t), intent(inout)       :: this
    type(particles_t), intent(inout)  :: p
    real(DP), intent(inout)           :: v(3, p%maxnatloc)
    real(DP), intent(in)              :: f(3, p%maxnatloc)
    real(DP), intent(inout)           :: dt

    ! ---

    integer   :: i

    real(DP)  :: vf, vg_dot_vg, Fg_dot_Fg, help

    ! ---

    call timer_start("fire_step2")


    !
    ! Communicate forces back to this processor if required
    !

#ifdef _MP
    if (mod_communicator%communicate_forces) then
       DEBUG_WRITE("- communicate_forces -")
       call communicate_forces(mod_communicator, p)
    endif
#endif


    !
    ! Integrate
    !

    !$omp parallel do default(none) &
    !$omp shared(f, p, v) &
    !$omp firstprivate(dt)
    do i = 1, p%natloc

       if (p%g(i) > 0) &
            VEC3(v, i) = VEC3(v, i) + 0.5_DP * VEC3(f, i) * dt

    enddo


    ! 
    ! Turn the global velocity a little bit 
    ! more along the global force...
    !

    vf         = 0.0_DP
    vg_dot_vg  = 0.0_DP
    Fg_dot_Fg  = 0.0_DP

    !$omp  parallel do default(none) &
    !$omp& shared(f, p, this, v) &
    !$omp& reduction(+:vf) reduction(+:vg_dot_vg) reduction(+:Fg_dot_Fg)
    do i = 1, p%natloc
       if (p%g(i) > 0) then
          vf = vf + dot_product(VEC3(v, i), VEC3(f, i))
          vg_dot_vg = vg_dot_vg + dot_product(VEC3(v, i), VEC3(v, i))
          Fg_dot_Fg = Fg_dot_Fg + dot_product(VEC3(f, i), VEC3(f, i))
       endif
    enddo

#ifdef _MP
    call sum_in_place(mod_communicator%mpi, vf)
    call sum_in_place(mod_communicator%mpi, vg_dot_vg)
    call sum_in_place(mod_communicator%mpi, Fg_dot_Fg)
#endif

    help = this%mix*sqrt(vg_dot_vg/Fg_dot_Fg)

    if (this%log) then
       call prlog("{FIRE} v.f = " // vf // ", v.v = " // vg_dot_vg // ", f.f = " // Fg_dot_Fg)
    endif

    !$omp  parallel do default(none) &
    !$omp& shared(f, help, p, this, v)
    do i = 1, p%natloc
       if (p%g(i) > 0) then
          VEC3(v, i) = (1-this%mix)*VEC3(v, i) + help*VEC3(f, i)
       else
          VEC3(v, i) = 0.0_DP
       endif
    enddo

    !-------------------------------------
    !
    ! Cut the velocities if the total
    ! power done by forces is negative
    !
    !-------------------------------------
    if (vf < 0.0_DP) then
       v          = 0.0_DP
       this%cut   = this%minsteps
       dt         = dt*this%decfac
       this%mix   = this%mix_in 
       this%cuts  = this%cuts + 1
    else 
       ! mixing is important only right after cut...
       if (this%cut < 0) then
          dt        = min(dt*this%incfac, this%max_dt)
          this%mix  = this%mix*this%mixdec 
       else
          this%cut  = this%cut - 1
       endif
    endif

    call timer_stop("fire_step2")
    
  endsubroutine fire_step2


  subroutine fire_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(fire_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)      :: cfg
    type(c_ptr), intent(out)     :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("FIRE"), &
         CSTR("The FIRE optimizer."))

    call ptrdict_register_integer_property(m, c_loc(this%minsteps), &
         CSTR("minsteps"), &
         CSTR("Minimum number of steps for mixing directly after a cut (i.e., after the velocities have been reset.)"))
    call ptrdict_register_real_property(m, c_loc(this%incfac), &
         CSTR("dt_inc"), &
         CSTR("Factor for increase of time step if no cuts occure."))
    call ptrdict_register_real_property(m, c_loc(this%decfac), CSTR("dt_dec"), &
         CSTR("Factor for decrease of time step directly after cut."))
    call ptrdict_register_real_property(m, c_loc(this%mix_in), CSTR("mix_in"), &
         CSTR("Start value for mixing parameter. Also used after a cut."))
    call ptrdict_register_real_property(m, c_loc(this%mixdec), CSTR("mix_dec"), &
         CSTR("Factor for decrease of the mixing parameter if no cut occurs."))
    call ptrdict_register_real_property(m, c_loc(this%max_dt), CSTR("max_dt"), &
         CSTR("Upper limit to the time step."))
    call ptrdict_register_real_property(m, c_loc(this%fmax), CSTR("fmax"), &
         CSTR("Convergence criterion."))

    call ptrdict_register_boolean_property(m, c_loc(this%log), CSTR("log"), &
         CSTR("Verbose logging."))

  endsubroutine fire_register

endmodule fire
