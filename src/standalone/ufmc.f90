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
!   classtype:ufmc_t classname:UFMC interface:integrators
! @endmeta
!
! Uniform-acceptance force-bias Monte Carlo (UFMC)
! Mees et al. PRB 85, 134301 (2012)
!

#include "macros.inc"

module ufmc
  use supplib
  use rng

  use particles

  implicit none

  private

  public :: ufmc_t
  type ufmc_t

    logical               :: firstcall
    real(DP)              :: T            ! [K]       temperature
    real(DP)              :: max_disp     ! [Ang]     maximum displacement of lightest atom
    real(DP), allocatable :: delta(:)     ! [Ang]     maximum displacement of each atom
    real(DP), allocatable :: gamma_F(:)   ! [Ang/eV]
    real(DP)              :: beta         ! [1/eV]    1/(kB*T)
    real(DP), allocatable :: r_prev(:,:)  ! [Ang]     positions in previous iteration 

  endtype ufmc_t


  public :: init
  interface init
     module procedure ufmc_init
  endinterface

  public :: del
  interface del
     module procedure ufmc_del
  endinterface

  public :: step1
  interface step1
     module procedure ufmc_step1
  endinterface

  public :: step2
  interface step2
     module procedure ufmc_step2
  endinterface

  public :: register
  interface register
    module procedure ufmc_register
  endinterface

contains

  !>
  !! Constructor
  !<
  subroutine ufmc_init(this, p)
    implicit none

    type(ufmc_t),      intent(inout) :: this
    type(particles_t), intent(inout) :: p

    ! ---


    this%beta = 1.0_DP / (this%T * BOLTZMANN_K)

    this%firstcall = .true.


  endsubroutine ufmc_init


  !>
  !! Destructor
  !<
  subroutine ufmc_del(this)
    implicit none

    type(ufmc_t), intent(inout) :: this

    ! ---

    if (allocated(this%delta  ))  deallocate(this%delta  )
    if (allocated(this%gamma_F))  deallocate(this%gamma_F)
    if (allocated(this%r_prev ))  deallocate(this%r_prev )

  endsubroutine ufmc_del


  subroutine ufmc_step1(this, p, v, f, dt, max_dt, max_dr, max_dr_sq)
    implicit none

    type(ufmc_t),      intent(inout)           :: this
    type(particles_t), intent(inout)           :: p
    real(DP),          intent(inout)           :: v(3, p%maxnatloc)
    real(DP),          intent(in)              :: f(3, p%maxnatloc)
    real(DP),          intent(inout)           :: dt
    real(DP),          intent(in),    optional :: max_dt
    real(DP),          intent(in),    optional :: max_dr
    real(DP),          intent(inout), optional :: max_dr_sq

    ! ---

    integer  :: i
    real(DP) :: m_min

    ! ---

    ! This has to be done here because the integrators are
    ! initialized before the atomic configuration is read.
    if (this%firstcall .eqv. .true.) then

      allocate(this%delta(p%nat))
      allocate(this%gamma_F(p%nat))
      allocate(this%r_prev(3,p%nat))

      this%r_prev = p%r_non_cyc

      m_min = minval(p%m(1:p%nat))

      do i=1,p%nat
        this%delta(i)   = this%max_disp * sqrt(m_min / p%m(i))
        this%gamma_F(i) = 0.5_DP * this%delta(i) * this%beta
      end do

      ! Estimated(!) timestep
      dt = this%max_disp / 3.0_DP * sqrt(0.5_DP*pi*m_min*this%beta)   ! [10.2fs]

      this%firstcall = .false.
    end if

  endsubroutine ufmc_step1


  subroutine ufmc_step2(this, p, v, f, dt)
    implicit none

    type(ufmc_t),      intent(inout) :: this
    type(particles_t), intent(inout) :: p
    real(DP),          intent(inout) :: v(3, p%maxnatloc)
    real(DP),          intent(in)    :: f(3, p%maxnatloc)
    real(DP),          intent(in)    :: dt

    ! ---

    integer  :: i, j
    real(DP) :: xi
    real(DP) :: P0
    real(DP) :: Pij
    real(DP) :: gamma
    real(DP) :: exp_g, exp_mg

    ! ---

    call timer_start("ufmc")

    do i=1,p%nat
      if (p%g(i) > 0) then
        j = 1
        do
          xi = rng_uniform(-1.0_DP, 1.0_DP)
          P0 = rng_uniform( 0.0_DP, 1.0_DP)

          gamma = f(j,i) * this%gamma_F(i)
        
          ! Approximation for Pij in case gamma
          ! is small (avoids division by 0).
          if (abs(gamma) .lt. 1e-10_DP) then
            if (xi .lt. 0.0_DP) then
              Pij = xi + 1.0_DP
            else if (xi .gt. 0.0_DP) then
              Pij = 1.0_DP - xi
            else
              Pij = 1.0_DP
            end if
          else
            exp_g  = exp( gamma)
            exp_mg = exp(-gamma)

            if (xi .lt. 0.0_DP) then
              Pij = (exp( gamma*(2.0_DP*xi + 1.0_DP) ) - exp_mg)/(exp_g - exp_mg)
            else if (xi .gt. 0.0_DP) then
              Pij = (exp_g - exp( gamma*(2.0_DP*xi - 1.0_DP) ))/(exp_g - exp_mg)
            else
              Pij = 1.0_DP
            end if
          end if

          ! Check if displacement is accepted.
          if (P0 .lt. Pij) then
            p%r_non_cyc(j,i) = p%r_non_cyc(j,i) + this%delta(i)*xi

            ! Estimate velocities from displacements.
            v(j,i) = (p%r_non_cyc(j,i) - this%r_prev(j,i)) / dt
            this%r_prev(j,i) = p%r_non_cyc(j,i)

            j = j+1
          else
            cycle
          end if

          if (j .gt. 3) exit
        end do
      end if
    end do

    call timer_stop("ufmc")
    
  endsubroutine ufmc_step2


  subroutine ufmc_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ufmc_t), target, intent(inout)  :: this
    type(c_ptr),          intent(in)     :: cfg
    type(c_ptr),          intent(out)    :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("UFMC"), &
         CSTR("Uniform-acceptance force-bias Monte Carlo."))

    call ptrdict_register_real_property(m, c_loc(this%T), CSTR("T"), &
         CSTR("Temperature"))
    call ptrdict_register_real_property(m, c_loc(this%max_disp), CSTR("max_disp"), &
         CSTR("Maximum displacement of the lightest atom."))

  endsubroutine ufmc_register

endmodule ufmc
