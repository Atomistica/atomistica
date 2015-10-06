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
!   private broken
!   classtype:heatflux_t classname:Heatflux
! @endmeta

!>
!! Muller-Plathe method for imposing a heat flux on an atomistic system
!!
!! Florian Muller-Plathe's method for imposing a heat flux on an atomistic
!! (for the computation of thermal conductivities).
!! See: F. Muller-Plathe, J. Chem. Phys. 106, 6082 (1997)
!>

#include "macros.inc"

module heatflux
  use supplib

  use particles
  use neighbors

  implicit none

  type heatflux_t

     integer            :: n_bins
     real(DP)           :: freq
     real(DP)           :: out_freq

     real(DP)           :: dx

     real(DP)           :: t
     real(DP)           :: out_t

     real(DP)           :: etransfer

     real(DP)           :: etransfer_tot

     integer            :: un

     !
     ! Velocities
     !

     real(DP), pointer  :: v(:, :)

  endtype heatflux_t

contains

  !**********************************************************************
  ! Initialize a heatflux object
  !**********************************************************************
  subroutine heatflux_particles_initialized(this, p)
    implicit none

    type(heatflux_t), intent(inout)  :: this
    type(particles_t), intent(in)    :: p

    ! ---

    integer  :: i

    ! ---

    write (ilog, '(A)')  "- heatflux_particles_initialized -"

    if (mod(this%n_bins, 2) /= 0) then
       EXIT_ON_ERROR("Number of bins must be even.", i)
    endif

    this%dx             = p%Abox(3, 3)/this%n_bins
    this%etransfer      = 0.0_DP
    this%etransfer_tot  = 0.0_DP

    this%t              = 0.0_DP
    this%out_t          = 0.0_DP

    this%un  = fopen("heatflux.out")

    write (this%un, '(A1,A9,3A20)')  "#", "it", "ti", "etransfer", "etransfer/t"

    write (ilog, '(5X,A,I10)')      "n_bins    = ", this%n_bins
    write (ilog, '(5X,A,F20.10)')   "freq      = ", this%freq
    write (ilog, '(5X,A,F20.10)')   "out_freq  = ", this%out_freq
    write (ilog, '(5X,A,F20.10)')   "dx        = ", this%dx
    write (ilog, *)

    call ptr_by_name(p%data, V_STR, this%v)

  endsubroutine heatflux_particles_initialized


  !**********************************************************************
  ! Delete a heatflux object
  !**********************************************************************
  subroutine heatflux_del(this, it, ti, dt, p, nl)
    implicit none

    type(heatflux_t), intent(inout)  :: this
    integer, intent(in)              :: it
    real(DP), intent(in)             :: ti
    real(DP), intent(in)             :: dt
    type(particles_t), intent(in)    :: p
    type(neighbors_t), intent(in)    :: nl

    ! ---

    integer  :: un

    ! ---

    write (ilog, '(A)')  "- heatflux_del -"
    write (ilog, '(5X,A,ES20.10)')  "Amount of transfered energy: ", this%etransfer_tot

    call fclose(this%un)

    write (ilog, *)

  endsubroutine heatflux_del


  !**********************************************************************
  ! Perform the measurement
  !**********************************************************************
  subroutine heatflux_invoke(this, dyn, nl)
    implicit none

    type(heatflux_t), intent(inout)  :: this
    type(dynamics_t), intent(inout)  :: dyn
    type(neighbors_t), intent(in)    :: nl

    ! ---

    integer   :: i, b, coldest, hottest
    real(DP)  :: ekin_coldest, ekin_hottest, ekin, v(3)

    ! ---

    call timer_start("heatflux_invoke")

    this%t      = this%t + dt
    this%out_t  = this%out_t + dt

    if (this%t >= this%freq) then

       coldest = -1
       hottest = -1

       ekin_coldest = -1
       ekin_hottest = -1

       do i = 1, p%nat
          b = int(POS(p, i, 3)/this%dx)

          if (b == 0) then
             ekin = p%m(i)*dot_product(VEC3(this%v, i), VEC3(this%v, i))

             if (ekin > ekin_hottest) then
                hottest = i
                ekin_hottest = ekin
             endif
          else if (b == this%n_bins/2) then
             ekin = p%m(i)*dot_product(VEC3(this%v, i), VEC3(this%v, i))

             if (ekin < ekin_coldest .or. ekin_coldest < 0) then
                coldest = i
                ekin_coldest = ekin
             endif

          endif

       enddo

       if (coldest == -1) then
          EXIT_ON_ERROR("No coldest particle found. Weird.", i)
       endif
       if (hottest == -1) then
          EXIT_ON_ERROR("No hottest particle found. Weird.", i)
       endif

       this%etransfer = this%etransfer + ekin_hottest - ekin_coldest

       v = VEC3(this%v, coldest)
       VEC3(this%v, coldest) = VEC3(this%v, hottest)
       VEC3(this%v, hottest) = v

       this%t  = 0.0_DP
    endif

    if (this%out_t >= this%out_freq) then
       write (this%un, '(I10,3ES20.10)')  it, ti, this%etransfer, this%etransfer/this%out_t

       this%etransfer_tot  = this%etransfer_tot + this%etransfer
       this%etransfer      = 0.0_DP
       this%out_t          = 0.0_DP
    endif

    call timer_stop("heatflux_invoke")

  endsubroutine heatflux_invoke


#ifdef MDCORE_MONOLITHIC

  !****************************************************************
  ! Initialize the property list
  !****************************************************************
  subroutine heatflux_register(this, cfg, m)
    implicit none

    type(heatflux_t), intent(inout)  :: this
    CPOINTER, intent(in)             :: cfg
    CPOINTER, intent(inout)          :: m

    ! ---

    this%n_bins      = 100
    this%freq        = 100.0_DP
    this%out_freq    = 100.0_DP

    m = ptrdict_register_section(cfg, CSTR("HeatFlux"), &
         CSTR("Impose a heat flux upon the system. This is Florian Muller-Plathe's method for the computation of thermal conductivity. See: F. Muller-Plathe, J. Chem. Phys. 106, 6083 (1997)"))

    call ptrdict_register_integer_property(m, this%n_bins, CSTR("n_bins"), &
         CSTR("Number of bins."))

    call ptrdict_register_real_property(m, this%freq, CSTR("freq"), &
         CSTR("Interval in which to transfer energy."))

    call ptrdict_register_real_property(m, this%out_freq, CSTR("out_freq"), &
         CSTR("Interval in which to output the amount of energy transfered to 'heatflux.out'."))

  endsubroutine heatflux_register

#endif

endmodule heatflux
