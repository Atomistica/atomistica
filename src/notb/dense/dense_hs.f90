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

#include "macros.inc"
#include "filter.inc"

module dense_hs
  use, intrinsic :: iso_c_binding

  use supplib

  use nonuniform_spline
  use timer
  use tls

  use particles
  use neighbors

  use materials

  use dense_hamiltonian_type
  use dense_hamiltonian

  implicit none

  private

  public :: hs_setup

contains

  !************************************************************
  !
  !                 function transform_orb
  !
  ! transform_orbs the (Hamiltonian or overlap) matrix elements
  ! (located in e(:)) for orbitals a and b, where c(:)
  ! are the direction cosines
  !
  ! Notation for the orbitals:
  !  s px py pz xy yz zx x2-y2 3z2-r2
  !  1  2  3  4  5  6  7   8      9
  !  0  1  1  1  2  2  2   2      2     (angular momentum=lo(:))
  !
  ! Notation for the orbital-integrals:
  !  dds ddp ddd pds pdp pps ppp sds sps sss
  !   1   2   3   4   5   6   7   8   9   10 
  !
  !************************************************************ 
  function transform_orb( a0, b0, c, e ) result(r)
    implicit none
    integer :: a0,b0
    real(8)  :: c(3),e(10)
    real(8)  :: r

    real(8)  :: l,m,n,ll,mm,nn
    real(8)  :: dds, ddp, ddd, pds, pdp, pps, ppp, sds, sps, sss
    integer :: f,a,b
    integer, save :: lo(9)
    real(8),  save :: s3
    logical, save :: first=.true.
    !$omp threadprivate(lo, s3, first)

    if( first ) then
       lo  = (/0,1,1,1,2,2,2,2,2/)
       s3  = sqrt(3.0) 
       first = .false.
    end if


    l   = c(1)
    m   = c(2)
    n   = c(3)
    ll  = l**2
    mm  = m**2
    nn  = n**2
    dds = e(1); ddp=e(2); ddd=e(3); pds=e(4); pdp=e(5)
    pps = e(6); ppp=e(7); sds=e(8); sps=e(9); sss=e(10)

    !--------------------------------------
    !
    ! The pre-factor comes from the 
    ! interchange of the orbitals
    ! with respect to tabulated values
    ! [parity = (-1)^(angular momentum)];
    ! tabulated MELs have always the lower
    ! angular momentum as the FIRST (bra-)
    ! state
    !
    !---------------------------------------
    f = 1
    a = a0
    b = b0
    if( a0>b0 ) then
       f = (-1)**( lo(a)+lo(b) )
       a = b0
       b = a0
    end if

    !-------------------------------------------
    !
    ! choose the transform_orbation rule according
    ! to the first orbital (a) and the second
    ! orbital (b) (now a<=b). This is the table
    ! of Slater & Koster, but filled with
    ! the missing rules (by permuting
    ! coordinates and dir. cosines, indicated
    ! by star) 
    !
    !-------------------------------------------
    select case(a)
    case(1)
       select case(b)
       case(1)
          r = sss
       case(2)
          r = l*sps
       case(3)
          r = m*sps
       case(4)
          r = n*sps
       case(5)
          r = s3*l*m*sds
       case(6) !(*)
          r = s3*m*n*sds
       case(7) !(*)
          r = s3*n*l*sds
       case(8)
          r = 0.5_DP*s3*(ll-mm)*sds
       case(9)
          r = (nn-0.5_DP*(ll+mm))*sds
       case default
          stop 'transf. not defined'
       end select
    case(2)
       select case(b)
       case(2)
          r = ll*pps + (1-ll)*ppp
       case(3)
          r = l*m*pps - l*m*ppp
       case(4)
          r = l*n*pps - l*n*ppp
       case(5)
          r = s3*ll*m*pds + m*(1-2*ll)*pdp
       case(6)
          r = s3*l*m*n*pds - 2*l*m*n*pdp
       case(7)
          r = s3*ll*n*pds + n*(1-2*ll)*pdp
       case(8)
          r = 0.5_DP*s3*l*(ll-mm)*pds + l*(1-ll+mm)*pdp
       case(9)
          r = l*(nn-0.5_DP*(ll+mm))*pds - s3*l*nn*pdp
       case default
          stop 'transf. not defined'
       end select
    case(3)
       select case(b)
       case(3) !(*)
          r = mm*pps + (1-mm)*ppp
       case(4) !(*)
          r = m*n*pps - m*n*ppp
       case(5) !(*)
          r = s3*mm*l*pds + l*(1-2*mm)*pdp
       case(6) !(*)
          r = s3*mm*n*pds + n*(1-2*mm)*pdp
       case(7) !(*)
          r = s3*m*n*l*pds - 2*m*n*l*pdp
       case(8)
          r = 0.5_DP*s3*m*(ll-mm)*pds - m*(1+ll-mm)*pdp
       case(9)
          r = m*(nn-0.5_DP*(ll+mm))*pds - s3*m*nn*pdp
       case default
          stop 'transf. not defined'
       end select
    case(4)
       select case(b)
       case(4) !(*)
          r = nn*pps + (1-nn)*ppp
       case(5) !(*)
          r = s3*l*m*n*pds - 2*m*n*l*pdp
       case(6) !(*)
          r = s3*nn*m*pds + m*(1-2*nn)*pdp
       case(7) !(*)
          r = s3*nn*l*pds + l*(1-2*nn)*pdp
       case(8)
          r = 0.5_DP*s3*n*(ll-mm)*pds - n*(ll-mm)*pdp
       case(9)
          r = n*(nn-0.5_DP*(ll+mm))*pds + s3*n*(ll+mm)*pdp
       case default
          stop 'transf. not defined'
       end select
    case(5)
       select case(b)
       case(5)
          r = 3*ll*mm*dds + (ll+mm-4*ll*mm)*ddp + (nn+ll*mm)*ddd
       case(6)
          r = 3*l*mm*n*dds + l*n*(1-4*mm)*ddp + l*n*(mm-1)*ddd
       case(7)
          r = 3*ll*m*n*dds + m*n*(1-4*ll)*ddp + m*n*(ll-1)*ddd
       case(8)
          r = 1.5*l*m*(ll-mm)*dds + 2*l*m*(mm-ll)*ddp + 0.5_DP*l*m*(ll-mm)*ddd
       case(9)
          r = s3*l*m*(nn-0.5_DP*(ll+mm))*dds - 2*s3*l*m*nn*ddp + 0.5_DP*s3*l*m*(1+nn)*ddd
       case default
          stop 'transf. not defined'
       end select
    case(6)
       select case(b)
       case(6) !(*)
          r = 3*mm*nn*dds + (mm+nn-4*mm*nn)*ddp + (ll+mm*nn)*ddd
       case(7) !(*)
          r = 3*m*nn*l*dds + m*l*(1-4*nn)*ddp + m*l*(nn-1)*ddd
       case(8)
          r = 1.5*m*n*(ll-mm)*dds - m*n*(1+2*(ll-mm))*ddp + m*n*(1+0.5_DP*(ll-mm))*ddd
       case(9)
          r = s3*m*n*(nn-0.5_DP*(ll+mm))*dds + s3*m*n*(ll+mm-nn)*ddp - 0.5_DP*s3*m*n*(ll+mm)*ddd
       case default
          stop 'transf. not defined'
       end select
    case(7)
       select case(b)
       case(7) !(*)
          r = 3*nn*ll*dds + (nn+ll-4*nn*ll)*ddp + (mm+nn*ll)*ddd
       case(8)
          r = 1.5*n*l*(ll-mm)*dds + n*l*(1-2*(ll-mm))*ddp - n*l*(1-0.5_DP*(ll-mm))*ddd
       case(9)
          r = s3*l*n*(nn-0.5_DP*(ll+mm))*dds + s3*l*n*(ll+mm-nn)*ddp - 0.5_DP*s3*l*n*(ll+mm)*ddd
       case default
          stop 'transf. not defined'
       end select
    case(8)
       select case(b)
       case(8)
          r = 0.75_DP*(ll-mm)**2*dds + (ll+mm-(ll-mm)**2)*ddp + (nn+0.25_DP*(ll-mm)**2)*ddd
       case(9)
          r = 0.5_DP*s3*(ll-mm)*(nn-0.5_DP*(ll+mm))*dds + s3*nn*(mm-ll)*ddp + 0.25_DP*s3*(1+nn)*(ll-mm)*ddd
       case default
          stop 'transf. not defined'
       end select
    case(9)
       select case(b)
       case(9)
          r = (nn-0.5_DP*(ll+mm))**2*dds + 3*nn*(ll+mm)*ddp + 0.75_DP*(ll+mm)**2*ddd
       case default
          stop 'transf. not defined'
       end select
    case default
       stop 'transf. not defined'
    end select

    r = r*f
  end function transform_orb


  subroutine hs_setup(this, db, p, nl, error)
    implicit none

    type(dense_hamiltonian_t), intent(inout)  :: this
    type(materials_t), intent(in)         :: db
    type(particles_t), intent(in)         :: p
    type(neighbors_t), intent(in)         :: nl
    integer, intent(inout), optional      :: error

    ! ---

    integer  :: k

    WF_T(DP), pointer  :: this_H(:, :, :), this_S(:, :, :)

    ! ---

    INIT_ERROR(error)

    call timer_start("hs_setup")

    call c_f_pointer(this%H, this_H, [this%norb, this%norb, this%nk])
    call c_f_pointer(this%S, this_S, [this%norb, this%norb, this%nk])

    this_H = 0.0_DP
    this_S = 0.0_DP

    do k = 1, this%nk

       call hs_setup_single_k( &
            this, db, p, nl, &
            this_H(:, :, k), this_S(:, :, k), &
            error)
       PASS_ERROR_WITH_INFO_AND_STOP_TIMER("H and S setup for k-point/spin number " // k // ".", "hs_setup", error)

    enddo

    call timer_stop("hs_setup")

  endsubroutine hs_setup


  !*******************************************************************
  !
  !                    subroutine setup_HS
  !
  ! calculates the Hamiltonian and overlap matrix elements
  ! when it is given the current positions of the atoms
  !
  !
  ! Modified for d-orbitals: Pekka Koskinen, May 2004
  !
  !*******************************************************************
  subroutine hs_setup_single_k(this, db, p, nl, H, S, error)
    implicit none

    type(dense_hamiltonian_t), intent(inout)  :: this
    type(materials_t),         intent(in)     :: db
    type(particles_t),         intent(in)     :: p
    type(neighbors_t),         intent(in)     :: nl
    WF_T(DP),                  intent(out)    :: H(this%norb, this%norb)
    WF_T(DP),                  intent(out)    :: S(this%norb, this%norb)
    integer,         optional, intent(inout)  :: error

    ! ---

    integer   :: i,ia0,ia,a,a0,j,noi,ni
    integer   :: eli,elj,x,y,z
    real(DP)  :: dr,c(3)
    real(DP)  :: vec(3)

    integer   :: error_loc

    type(notb_element_t), pointer  :: this_at(:)

    ! ---

    call c_f_pointer(this%at, this_at, [this%nat])

    H  = 0.0_DP   ! tls_mat1
    S  = 0.0_DP   ! tls_mat2

    error_loc = ERROR_NONE

    !$omp  parallel default(none) &
    !$omp& private(a, c, dr, eli, elj, i, ia, ia0, ni, noi, j, vec, x, y, z) &
    !$omp& shared(db, nl, H, p, S, this, this_at) &
    !$omp& reduction(+:error_loc)

#ifdef _OPENMP
    call tls_init(this%norb, mat=2)
#else
#define tls_mat1 H
#define tls_mat2 S
#endif
    
    !$omp do
    i_loop: do i = 1, p%nat

       if (IS_EL(this%f, p, i) .and. error_loc == ERROR_NONE) then

          noi = this_at(i)%no    ! number of atomic orbitals
          ia0 = this_at(i)%o1    ! first orbital (in global list)

          do a0 = 1, noi
             ia = ia0 + a0-1
             a  = get_orbital(noi, a0)

             tls_mat1(ia, ia) = this_at(i)%e(a)    ! orbital energy
             tls_mat2(ia, ia) = 1.0_DP             ! orbitals are normalized
          enddo

          eli = this_at(i)%enr   ! internal element number

          ni_loop: do ni = nl%seed(i), nl%last(i)

             j = GET_NEIGHBOR(nl, ni)

             if (IS_EL(this%f, p, j) .and. error_loc == ERROR_NONE) then

                DIST(p, nl, i, ni, vec, dr)

#ifdef SPARSE
                if (i <= j .and. dr < this%rho_cutoff) then
#else
#ifdef LAMMPS
                if (p%tag(i) < p%tag(j)) then
#else
                if (i <= j) then
#endif
#endif

                   elj = this_at(j)%enr   ! internal element number

                   if (dr < db%cut(eli, elj)) then
                      c = -vec / dr

#ifdef LAMMPS
                      x = 0
                      y = 0
                      z = 0
#else
                      x = VEC(nl%dc, ni, 1)
                      y = VEC(nl%dc, ni, 2)
                      z = VEC(nl%dc, ni, 3)
#endif

                      call calc_matrix_elements( &
                           this, db, &
                           this_at(i), this_at(j), dr, c, &
                           tls_mat1, tls_mat2, &
                           error_loc)
                      TRACE_DELAYED_ERROR_WITH_INFO("Computation of matrix element between atom " // i // " and atom " // j // ".", error_loc)
                   endif

                endif

             endif

          enddo ni_loop

       endif

    enddo i_loop

#ifdef _OPENMP
    call tls_reduce(this%norb, mat1=H, mat2=S)
#else
#undef H
#undef S
#endif

    !$omp end parallel

    INVOKE_DELAYED_ERROR(error_loc, error)

  endsubroutine hs_setup_single_k


  subroutine calc_matrix_elements(this, db, el_i, el_j, my_r_ij, &
       n_ij, H, S, error)
    implicit none

    ! ---

    type(dense_hamiltonian_t), intent(in)  :: this
    type(materials_t), intent(in)          :: db

    type(notb_element_t), intent(in)       :: el_i         ! element information for first element
    type(notb_element_t), intent(in)       :: el_j         ! element information for second element
    real(DP), intent(in)                   :: my_r_ij         ! distance of atom i and atom j (|r_i - r_j|)
    real(DP), intent(in)                   :: n_ij(3)      ! normal vector (r_i - r_j/|r_i - r_j|), direction cosines
    WF_T(DP), intent(inout)                :: H(this%norb, this%norb)
    WF_T(DP), intent(inout)                :: S(this%norb, this%norb)
    integer, intent(inout), optional       :: error

    ! ---
    
    integer  :: list(0:10,9) = -1

    integer  :: noi, noj  ! number of orbitals
    integer  :: eli, elj  ! element numbers
    integer  :: nr        ! number of contributing orbital combinations
    integer  :: bo        ! bond (orbital combination) index
    integer  :: ia, ia0, jb, jb0  ! orbital numbers (in global matrix)
    integer  :: a, b, a0, b0, q, m
    real(8)  :: he_ij(10), se_ij(10), he_ji(10), se_ji(10)

    WF_T(DP)  :: H_el, S_el

    ! ---

    ! list tells the bond integrals needed when you know the maximum of number
    ! of orbitals of the atom pair the zeroth position tells the number of MELs
    !                   dds ddp ddd pds pdp pps ppp sds sps sss
    !                    1   2   3   4   5   6   7   8   9   10
    list(0:1,  1) = [1,                                      10] ! s
    list(0:2,  3) = [2,                      6,  7             ] ! p
    list(0:4,  4) = [4,                      6,  7,      9,  10] ! sp
    list(0:3,  5) = [3,  1,  2,  3                             ] ! d
    list(0:5,  6) = [5,  1,  2,  3,                  8,      10] ! sd
    list(0:7,  8) = [7,  1,  2,  3,  4,  5,  6,  7             ] ! pd
    list(0:10, 9) = [10, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10] ! spd

    eli = el_i%enr
    elj = el_j%enr
    noi = el_i%no
    noj = el_j%no
    ia0 = el_i%o1
    jb0 = el_j%o1

    m  = max(noi, noj)
    nr = list(0, m)                    ! number of non-vanishing matrix elements
    a  = interval(db%HS(eli, elj), my_r_ij, error)  ! determine the interval dr can be found in (within xs)
    PASS_ERROR(error)
    b  = interval(db%HS(elj, eli), my_r_ij, error)  ! determine the interval dr can be found in (within xs)
    PASS_ERROR(error)
    do q = 1, nr
       bo = list(q, m)   ! bond order
       if (bo <= 0) then
          RAISE_ERROR("bo <= 0!", error)
       endif
       he_ij(bo) = f(db%HS(eli, elj), bo,          my_r_ij, a)
       se_ij(bo) = f(db%HS(eli, elj), bo+MAX_NORB, my_r_ij, a)
       he_ji(bo) = f(db%HS(elj, eli), bo,          my_r_ij, b)
       se_ji(bo) = f(db%HS(elj, eli), bo+MAX_NORB, my_r_ij, b)
    enddo

!    call timer('setup_hs_1',2)

    !---------------------------------------
    ! 
    ! Now, calculate matrix element (ia,jb)
    ! using the Slater-Koster
    ! transform_orbation rules. This is fast
    ! procedure.
    !
    !---------------------------------------
!    call timer('setup_hs_2',1)
    a_loop: do a0 = 1, noi
       ia = ia0 + a0-1
       a  = get_orbital(noi, a0)
       b_loop: do b0 = 1, noj
          jb = jb0 + b0-1
          b  = get_orbital(noj, b0)
          !------------------------------------------------------- 
          ! if b>a (i.e. ang.momenta l_b>l_a), we must use
          ! the other table
          !------------------------------------------------------- 
          if (a <= b) then
             H_el = transform_orb( a, b, n_ij, he_ij )
             S_el = transform_orb( a, b, n_ij, se_ij )
          else
             H_el = transform_orb( a, b, n_ij, he_ji )
             S_el = transform_orb( a, b, n_ij, se_ji )
          endif

          ! Now the order becomes important!!!
          H(ia, jb) = H(ia, jb) + H_el
          S(ia, jb) = S(ia, jb) + S_el

          ! Fixme!!! It would suffice if this is done once at the end
          ! of setup_HS when H and S is setup
          if (ia0 /= jb0) then
             H(jb, ia) = H(ia, jb)
             S(jb, ia) = S(ia, jb)
          endif

       enddo b_loop
    enddo a_loop

!    call timer('setup_hs_2',2)
  endsubroutine calc_matrix_elements

endmodule dense_hs
