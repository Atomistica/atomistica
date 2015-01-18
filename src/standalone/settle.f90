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
!   classtype:settle_t classname:SETTLE interface:integrators
! @endmeta

!>
!! SETTLE method for constraint water molecules (i.e., 3 rigid point charges)
!!
!! SETTLE method for constraint water molecules (i.e., 3 rigid point charges)
!<


#include "macros.inc"
#include "filter.inc"

module settle
  use supplib

  use particles
  use filter
  use molecules
  use dynamics

  use verlet_support

#ifdef _MP
  use communicator
#endif

  implicit none

  private

  integer, parameter :: H_ = 1
  integer, parameter :: O_ = 8

  public :: settle_t
  type settle_t

     type(particles_t), pointer  :: p  => NULL()

     character(MAX_EL_STR)       :: elements  = "*"
     integer                     :: els

     real(DP)                    :: d_OH  = 1.0_DP
     real(DP)                    :: d_HH  = 1.6329931618554521_DP

     real(DP)                    :: mO, mH
     real(DP)                    :: mOrmT, mHrmT
     real(DP)                    :: ra, rb, rc, rra

     logical, allocatable        :: done(:)

     type(molecules_t)           :: molecules

  endtype settle_t


  public :: init
  interface init
     module procedure settle_init
  endinterface

  public :: del
  interface del
     module procedure settle_del
  endinterface

  public :: step1_with_dyn
  interface step1_with_dyn
     module procedure settle_step1
  endinterface

  public :: step2_with_dyn
  interface step2_with_dyn
     module procedure settle_step2
  endinterface

  public :: register
  interface register
    module procedure settle_register
  endinterface

!--- Internal

  interface set_particles
     module procedure settle_set_particles
  endinterface


!  private ssqrt

contains

  !>
  !! Normalize a vector
  !<
  pure subroutine normalize_vector(a)
    implicit none

    real(DP), intent(inout) :: a(3)

    ! ---

    a  = a / sqrt(dot_product(a, a))

  endsubroutine normalize_vector


  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine settle_init(this, p, error)
    implicit none

    type(settle_t), intent(inout)     :: this
    type(particles_t), intent(inout)  :: p
    integer, intent(out), optional    :: error

    ! ---

    real(DP)         :: rmT, t1

    ! ---

    INIT_ERROR(error)

    this%mO     = ElementMass(O_)
    this%mH     = ElementMass(H_)

    rmT         = 1.0_DP / (this%mO+2*this%mH)
    this%mOrmT  = this%mO * rmT
    this%mHrmT  = this%mH * rmT
    t1          = 0.5_DP*this%mO/this%mH
    this%rc     = 0.5_DP*this%d_HH
    this%ra     = sqrt(this%d_OH*this%d_OH-this%rc*this%rc)/(1.0+t1)
    this%rb     = t1*this%ra
    this%rra    = 1.0 / this%ra

    call init(this%molecules, p)

  endsubroutine settle_init


  !>
  !! Set a particles object
  !!
  !! Set a particles object
  !<
  subroutine settle_set_particles(this, p)
    implicit none

    type(settle_t), intent(inout)  :: this
    type(particles_t), target      :: p

    ! ---

    integer       :: i

    ! ---

    write (ilog, '(A)') "- settle_set_particles -"

    this%p  => p

#ifdef _MP
    ! The forces need to be correct for a water right at the boundary
    call request_border(mod_communicator, p, max(this%d_OH, this%d_HH))
    mod_communicator%communicate_forces  = .true.
#endif

    !
    ! Adjust the degrees of freedom
    !

    do i = 1, p%natloc
       if (p%Z(i) == O_ .and. this%molecules%next(i) /= 0) then !p%Z(i+1) == H_ .and. p%Z(i+2) == H_) then
          p%dof = p%dof - 3
       endif
    enddo

    allocate(this%done(p%maxnatloc))

    this%els  = filter_from_string(this%elements, p)
    
    write (ilog, '(5X,A,I10)')  "dof = ", p%dof
    write (ilog, *)

  endsubroutine settle_set_particles


  !>
  !! Destructor
  !!
  !! Deallocates the done-pointer and the molecules object
  !<
  subroutine settle_del(this)
    implicit none

    type(settle_t), intent(inout)  :: this

    ! ---

    if (allocated(this%done)) then
       deallocate(this%done)
    endif

    call del(this%molecules)

  endsubroutine settle_del


  !>
  !! Save sqrt
  !!
  !! Return 0 if x < 0
  !<
!!$  real(DP) function ssqrt(x)
!!$    implicit none
!!$
!!$    real(DP), intent(in)  :: x
!!$
!!$    ! ---
!!$
!!$    if (x <= 1d-12) then
!!$       ssqrt  = 1d-6
!!$    else
!!$       ssqrt  = sqrt(x)
!!$    endif
!!$
!!$  endfunction ssqrt

#define ssqrt(x) sqrt(x)


  !>
  !! First integration step
  !!
  !! First integration step
  !<
  recursive subroutine settle_step1(this, dyn, max_dt, max_dr, max_dr_sq) 
    implicit none

    type(settle_t), intent(inout)      :: this
    type(dynamics_t), intent(inout)    :: dyn
    real(DP), intent(in), optional     :: max_dt
    real(DP), intent(in), optional     :: max_dr
    real(DP), intent(inout), optional  :: max_dr_sq

    ! ---

    type(particles_t), pointer :: p

    integer   :: i0, i1, i2

    real(DP)  :: p0(3), p1(3), p2(3), q0(3), q1(3), q2(3)
    real(DP)  :: b0(3), c0(3), d0(3), a1(3), b1(3), c1(3), a2(3), b2(3), c2(3), a3(3), b3(3), c3(3)
    real(DP)  :: dv0(3), dv1(3), dv2(3), dr0(3), dr1(3), dr2(3)
    real(DP)  :: n(3, 3), n0(3), n1(3), n2(3)
    real(DP)  :: A1Z, sinphi, tmp, tmp1, tmp2, cosphi, sinpsi, cospsi
    real(DP)  :: alpha, beta, gamma, a2b2, sintheta, costheta

    real(DP)  :: w(3, 3)

    real(DP)  :: d2t, l_max_dr_sq

    logical   :: is_water

    ! ---

    call timer_start("settle_step1")

    p => dyn%p

    if (.not. associated(this%p, p)) then
       call set_particles(this, p)
    endif

    d2t          = dyn%dt**2

    l_max_dr_sq  = 0.0_DP

    ! Taken from NAMD! Check license

    !
    ! Here we assume that each oxygen starts a water molecule,
    ! i.e. the following two atoms need to be hydrogens
    !

    this%done  = .false.
    w          = 0.0_DP

!!    !$omp do
    do i0 = 1, p%nat
       ! Find water molecules and apply constraint

       if (IS_EL(this%els, p, i0)) then

          is_water = .false.
          if (p%Z(i0) == O_ .and. this%molecules%next(i0) > 0) then
             i1 = p%global2local(this%molecules%next(i0))
             if (i1 > 0) then
                if (this%molecules%next(i1) > 0) then
                   i2 = p%global2local(this%molecules%next(i1))
                   if (i2 > 0) then
                      is_water = .true.
                   endif
                endif
             endif
          endif

          if (is_water) then

             ! vectors in the plane of the original positions
             d0     = PNC3(p, i1) - PNC3(p, i0)
             b0     = in_bounds(p, d0)
             d0     = PNC3(p, i2) - PNC3(p, i0)
             c0     = in_bounds(p, d0)

             q0     = VEC3(dyn%v, i0) + 0.5_DP * VEC3(dyn%f, i0) / p%m(i0) * dyn%dt
             q1     = VEC3(dyn%v, i1) + 0.5_DP * VEC3(dyn%f, i1) / p%m(i1) * dyn%dt
             q2     = VEC3(dyn%v, i2) + 0.5_DP * VEC3(dyn%f, i2) / p%m(i2) * dyn%dt

             p0     = PNC3(p, i0)       + q0 * dyn%dt
             p1     = PNC3(p, i0) + b0  + q1 * dyn%dt
             p2     = PNC3(p, i0) + c0  + q2 * dyn%dt

             ! new center of mass
             d0     = p0*this%mOrmT &
                  + ((p1+p2)*this%mHrmT)

             a1     = p0 - d0
             b1     = p1 - d0
             c1     = p2 - d0

             ! Vectors describing transformation from original coordinate system to
             ! the 'primed' coordinate system as in the diagram.
             n0     = cross_product(b0, c0)
             b2     = b0+c0
             n1     = cross_product(n0, b2)
             n2     = cross_product(n0, n1)

             call normalize_vector(n0)
             call normalize_vector(n1)
             call normalize_vector(n2)

             n(1, :)  = n1
             n(2, :)  = n2
             n(3, :)  = n0

             b0       = matmul(n, b0)
             c0       = matmul(n, c0)

             A1Z      = dot_product(n0, a1)
             b1       = matmul(n, b1)
             c1       = matmul(n, c1)

             ! now we can compute positions of canonical water
             sinphi    = A1Z * this%rra
             tmp       = 1.0_DP - sinphi*sinphi
             cosphi    = ssqrt(tmp)
             sinpsi    = (b1(3) - c1(3))/(2.0_DP*this%rc*cosphi)
             tmp       = 1.0_DP - sinpsi*sinpsi
             if (tmp < 0.0_DP) then
                call particles_dump_info(p, i0)
                call particles_dump_info(p, i1)
                call particles_dump_info(p, i2)
             endif
             cospsi    = ssqrt(tmp)

             tmp1      = this%rc*sinpsi*sinphi
             tmp2      = this%rc*sinpsi*cosphi

             a2        = (/          0.0_DP,   this%ra*cosphi,          this%ra*sinphi        /)
             b2        = (/ -this%rc*cospsi,  -this%rb*cosphi - tmp1,  -this%rb*sinphi + tmp2 /)
             c2        = (/  this%rc*cosphi,  -this%rb*cosphi + tmp1,  -this%rb*sinphi - tmp2 /)

             ! there are no a0 terms because we've already subtracted the term off 
             ! when we first defined b0 and c0.
             alpha     = b2(1)*(b0(1) - c0(1)) + b0(2)*b2(2) + c0(2)*c2(2)
             beta      = b2(1)*(c0(2) - b0(2)) + b0(1)*b2(2) + c0(1)*c2(2)
             gamma     = b0(1)*b1(2) - b1(1)*b0(2) + c0(1)*c1(2) - c1(1)*c0(2)

             a2b2      = alpha*alpha + beta*beta
             sintheta  = (alpha*gamma - beta*ssqrt(a2b2 - gamma*gamma))/a2b2
             costheta  = ssqrt(1.0_DP - sintheta*sintheta)

!          write (*, *)  acos(cosphi)*180/PI, acos(cospsi)*180/PI, acos(costheta)*180/PI

             a3        = (/                  -a2(2)*sintheta,                    a2(2)*costheta,  a2(3) /)
             b3        = (/  b2(1)*costheta - b2(2)*sintheta,   b2(1)*sintheta + b2(2)*costheta,  b2(3) /)
             c3        = (/ -b2(1)*costheta - c2(2)*sintheta,  -b2(1)*sintheta + c2(2)*costheta,  c2(3) /)

             n         = transpose(n)
             dr0       = matmul(n, a3) + d0 - p0
             dr1       = matmul(n, b3) + d0 - p1
             dr2       = matmul(n, c3) + d0 - p2

             dv0       = dr0 / dyn%dt
             dv1       = dr1 / dyn%dt
             dv2       = dr2 / dyn%dt

             dr0       = dr0 + q0 * dyn%dt
             dr1       = dr1 + q1 * dyn%dt
             dr2       = dr2 + q2 * dyn%dt

#ifndef IMPLICIT_R
             POS3(p, i0)  = POS3(p, i0) + dr0
             POS3(p, i1)  = POS3(p, i1) + dr1
             POS3(p, i2)  = POS3(p, i2) + dr2
#endif

             PNC3(p, i0)  = PNC3(p, i0) + dr0
             PNC3(p, i1)  = PNC3(p, i1) + dr1
             PNC3(p, i2)  = PNC3(p, i2) + dr2

             VEC3(dyn%v, i0)  = q0 + dv0
             VEC3(dyn%v, i1)  = q1 + dv1
             VEC3(dyn%v, i2)  = q2 + dv2

             dv0          = dv0*p%m(i0)/dyn%dt
             dv1          = dv1*p%m(i1)/dyn%dt
             dv2          = dv2*p%m(i2)/dyn%dt

             w            = w &
                  + outer_product(p0, dv0) &
                  + outer_product(p1, dv1) &
                  + outer_product(p2, dv2)

             l_max_dr_sq = max(l_max_dr_sq, &
                  maxval( &
                  (/ dot_product(dr0, dr0), dot_product(dr1, dr1), dot_product(dr2, dr2) /) &
                  ) &
                  )

             this%done(i0) = .true.
             this%done(i1) = .true.
             this%done(i2) = .true.

          endif

       endif

    enddo

    ! Fixme!!! Parallelize top as well

    !$omp  parallel default(none) &
    !$omp& private(A1Z, alpha, beta, gamma) &
    !$omp& private(a1, a2, a2b2, a3, b0, b1, b2, b3) &
    !$omp& private(c0, c1, c2, c3, d0) &
    !$omp& private(dr0, dr1, dr2, dv0, dv1, dv2) &
    !$omp& private(i0, i1, i2, is_water) &
    !$omp& private(n, n0, n1, n2) &
    !$omp& private(p0, p1, p2, q0, q1, q2) &
    !$omp& private(sinphi, cosphi, sinpsi, cospsi, sintheta, costheta) &
    !$omp& private(tmp, tmp1, tmp2) &
    !$omp& shared(dyn, d2t, p, this) &
    !$omp& reduction(max:l_max_dr_sq) reduction(+:w)

    !$omp do
    do i0 = 1, p%natloc
          
       if (.not. this%done(i0)) then

          if (p%g(i0) > 0 .and. IS_EL(this%els, p, i0)) then
                
             dr0          = VEC3(dyn%v, i0) * dyn%dt + 0.5_DP * VEC3(dyn%f, i0) / p%m(i0) * d2t
#ifndef IMPLICIT_R
             POS3(p, i0)  = POS3(p, i0) + dr0
#endif
             PNC3(p, i0)  = PNC3(p, i0) + dr0
             VEC3(dyn%v, i0)  = VEC3(dyn%v, i0) + 0.5_DP * VEC3(dyn%f, i0) / p%m(i0) * dyn%dt

             l_max_dr_sq  = max(l_max_dr_sq, dot_product(dr0, dr0))

          endif

       endif

    enddo

    !$omp end parallel

    !
    ! Maximum particle displacement
    !

    dyn%wpot  = dyn%wpot + w

    p%accum_max_dr  = p%accum_max_dr + sqrt(l_max_dr_sq)

    if (present(max_dr_sq)) then
       max_dr_sq  = max(max_dr_sq, l_max_dr_sq)
    endif

    call I_changed_positions(p)

    call timer_stop("settle_step1")

  endsubroutine settle_step1


  !>
  !! Second integration step
  !!
  !! Second integration step, velocity constraint
  !<
  recursive subroutine settle_step2(this, dyn)
    implicit none

    type(settle_t),   intent(inout) :: this
    type(dynamics_t), intent(inout) :: dyn

    ! ---

    type(particles_t), pointer :: p

    integer   :: i0, i1, i2

    real(DP)  :: rAB(3), rBC(3), rCA(3), AB(3), BC(3), CA(3)
    real(DP)  :: vab, vbc, vca, ga(3), gb(3), gc(3)
    real(DP)  :: cosA, cosB, cosC, ma, mb, mab, d, tab, tbc, tca, dt2
    real(DP)  :: w(3, 3)

    logical   :: is_water

    ! ---

    call timer_start("settle_step2")

    p => dyn%p

    !
    ! Communicate forces back to this processor if required
    !

#ifdef _MP
    if (mod_communicator%communicate_forces) then
       DEBUG_WRITE("- communicate_forces -")
       call communicate_forces(mod_communicator, p)
    endif
#endif


    dt2 = 2*dyn%dt

    ! Taken from NAMD! Check license

    !
    ! Here we assume that each oxygen starts a water molecule,
    ! i.e. the following two atoms need to be hydrogens
    !

    this%done  = .false.
    w          = 0.0_DP

    !$omp  parallel default(none) &
    !$omp& private(i0, i1, i2, is_water) &
    !$omp& private(AB, BC, CA, rAB, rBC, rCA) &
    !$omp& private(cosa, cosb, cosc, d) &
    !$omp& private(ga, gb, gc) &
    !$omp& private(ma, mb, mab) &
    !$omp& private(tab, tbc, tca) &
    !$omp& private(vab, vbc, vca) &
    !$omp& shared(dyn, dt2, p, this) &
    !$omp& reduction(+:w)

    !$omp do
    do i0 = 1, p%nat
       ! Find water molecules and apply constraint

       if (IS_EL(this%els, p, i0)) then

          is_water = .false.
          if (p%Z(i0) == O_ .and. this%molecules%next(i0) > 0) then
             i1 = p%global2local(this%molecules%next(i0))
             if (i1 > 0) then
                if (this%molecules%next(i1) > 0) then
                   i2 = p%global2local(this%molecules%next(i1))
                   if (i2 > 0) then
                      is_water = .true.
                   endif
                endif
             endif
          endif

          if (is_water) then
             VEC3(dyn%v, i0)  = VEC3(dyn%v, i0) + 0.5_DP * VEC3(dyn%f, i0) / p%m(i0) * dyn%dt
             VEC3(dyn%v, i1)  = VEC3(dyn%v, i1) + 0.5_DP * VEC3(dyn%f, i1) / p%m(i1) * dyn%dt
             VEC3(dyn%v, i2)  = VEC3(dyn%v, i2) + 0.5_DP * VEC3(dyn%f, i2) / p%m(i2) * dyn%dt

             AB    = POS3(p, i1) - POS3(p, i0)
             BC    = POS3(p, i2) - POS3(p, i1)
             CA    = POS3(p, i0) - POS3(p, i2)

             rAB   = in_bounds(p, AB)
             rBC   = in_bounds(p, BC)
             rCA   = in_bounds(p, CA)

             AB    = rAB
             BC    = rBC
             CA    = rCA

             call normalize_vector(AB)
             call normalize_vector(BC)
             call normalize_vector(CA)

             cosA    = -dot_product(AB, CA)
             cosB    = -dot_product(BC, AB)
             cosC    = -dot_product(CA, BC)
             
             vab     = dot_product((VEC3(dyn%v, i1) - VEC3(dyn%v, i0)), AB)
             vbc     = dot_product((VEC3(dyn%v, i2) - VEC3(dyn%v, i1)), BC)
             vca     = dot_product((VEC3(dyn%v, i0) - VEC3(dyn%v, i2)), CA)

             ma      = this%mO
             mb      = this%mH
             mab     = ma+mb

             d       = (2*mab*mab + 2*ma*mb*cosA*cosB*cosC - 2*mb*mb*cosA*cosA &
                  - ma*mab*(cosB*cosB + cosC*cosC))*0.5_DP/mb

             tab     = (vab*(2*mab - ma*cosC*cosC) + &
                  vbc*(mb*cosC*cosA - mab*cosB) + &
                  vca*(ma*cosB*cosC - 2*mb*cosA))*ma/d

             tbc     = (vbc*(mab*mab - mb*mb*cosA*cosA) + &
                  vca*ma*(mb*cosA*cosB - mab*cosC) + &
                  vab*ma*(mb*cosC*cosA - mab*cosB))/d

             tca     = (vca*(2*mab - ma*cosB*cosB) + &
                  vab*(ma*cosB*cosC - 2*mb*cosA) + &
                  vbc*(mb*cosA*cosB - mab*cosC))*ma/d

             AB   = tab*AB
             BC   = tbc*BC
             CA   = tca*CA

             ga   = AB - CA
             gb   = BC - AB
             gc   = CA - BC

             VEC3(dyn%v, i0)  = VEC3(dyn%v, i0) + (0.5_DP/ma)*ga
             VEC3(dyn%v, i1)  = VEC3(dyn%v, i1) + (0.5_DP/mb)*gb
             VEC3(dyn%v, i2)  = VEC3(dyn%v, i2) + (0.5_DP/mb)*gc

             AB   = AB/dt2
             BC   = BC/dt2
             CA   = CA/dt2

             w    = w &
                  - outer_product(rAB, AB) &
                  - outer_product(rBC, BC) &
                  - outer_product(rCA, CA)

             this%done(i0) = .true.
             this%done(i1) = .true.
             this%done(i2) = .true.

          endif

       endif

    enddo

    !$omp do
    do i0 = 1, p%natloc

       if (.not. this%done(i0)) then

          if (p%g(i0) > 0 .and. IS_EL(this%els, p, i0)) &
               VEC3(dyn%v, i0) = VEC3(dyn%v, i0) + 0.5_DP * VEC3(dyn%f, i0) / p%m(i0) * dyn%dt

       endif

    enddo

    !$omp end parallel

    dyn%wpot  = dyn%wpot + w


    !
    ! Update virial and kinetic energy
    !

!    call compute_kinetic_energy_and_virial(p)

    call timer_stop("settle_step2")

  endsubroutine settle_step2


  subroutine settle_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(settle_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)        :: cfg
    type(c_ptr), intent(out)       :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("SETTLE"), &
         CSTR("SETTLE method for constraint simulation of water molecules (see: S. Miyamoto and P. A. Kollman, J. Comput. Chem. 13, 952 (1992))."))

    call ptrdict_register_string_property(m, c_loc(this%elements(1:1)), &
         MAX_EL_STR, &
         CSTR("elements"), &
         CSTR("Elements for which to enable this integrator."))

    call ptrdict_register_real_property(m, c_loc(this%d_OH), CSTR("d_OH"), &
         CSTR("O-H distance in water molecule."))
    call ptrdict_register_real_property(m, c_loc(this%d_HH), CSTR("d_HH"), &
         CSTR("H-H distance in water molecule."))

  endsubroutine settle_register

endmodule settle
