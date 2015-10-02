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

!**********************************************************************
! Tight-binding force calculation
!**********************************************************************

#include "macros.inc"
#include "filter.inc"

module dense_force_calc
  use, intrinsic :: iso_c_binding

  use supplib
  
  use particles
  use neighbors

  use materials

  use dense_hamiltonian_type
  use dense_hamiltonian

  implicit none

  private

  public  :: forces

contains

  !**********************************************************************
  ! Force calculation including k-space summation
  !**********************************************************************
  subroutine forces(p, nl, tb, db, f, wpot, wpot_per_bond, error)
    implicit none

    type(particles_t), intent(inout)   :: p
    type(neighbors_t), intent(in)      :: nl
    type(dense_hamiltonian_t), intent(in)  :: tb
    type(materials_t), intent(in)      :: db
    real(DP), intent(inout)            :: f(3, p%maxnatloc)
    real(DP), intent(inout)            :: wpot(3, 3)
#ifdef LAMMPS
    real(DP), intent(inout), optional  :: wpot_per_bond(6, nl%neighbors_size)
#else
    real(DP), intent(inout), optional  :: wpot_per_bond(3, 3, nl%neighbors_size)
#endif
    integer, intent(inout), optional   :: error

    ! ---

    integer   :: list(0:10,9) = -1

    integer   :: ni, I, J, mu, nu, mu0, nu0, kk
    integer   :: elI, elJ, noI, noJ, Imu, Jnu
    real(DP)  :: rIJ(3), abs_rIJ, l(3)
    real(DP)  :: Hdiff_ij(3,9,9),Sdiff_ij(3,9,9),Hdiff_ji(3,9,9),Sdiff_ji(3,9,9)
    real(DP)  :: dH_ij(-10:10),dS_ij(-10:10),dH_ji(-10:10),dS_ji(-10:10)
    real(DP)  :: w(3, 3), wij(3, 3)
    WF_T(DP)  :: Ft(3), F2(3), rho_ImuJnu, e_ImuJnu
    integer   :: lo(9), m, nr, q, a, b

    integer   :: k!, tmpun

    integer   :: error_loc

    WF_T(DP),             pointer  :: tb_rho(:, :, :), tb_e(:, :, :)
    type(notb_element_t), pointer  :: tb_at(:)

    ! ---

    INIT_ERROR(error)

    call timer_start("forces")

    call c_f_pointer(tb%rho, tb_rho, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%e, tb_e, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%at, tb_at, [tb%nat])

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

    lo     = (/0,1,1,1,2,2,2,2,2/)

    w      = 0.0_DP

    error_loc  = ERROR_NONE

    !$omp  parallel default(none) &
    !$omp& private(a, abs_rIJ, b, dH_ij, dH_ji, dS_ij, dS_ji, e_ImuJnu) &
    !$omp& private(elI, elJ, F2, Ft, Hdiff_ij, Hdiff_ji) &
    !$omp& private(I, Imu, mu, J, Jnu, nu, kk, l, m, ni) &
    !$omp& private(noI, noJ, nr, q, rho_ImuJnu, rIJ) &
    !$omp& private(Sdiff_ij, Sdiff_ji, wij) &
    !$omp& shared(db, f, list, lo, nl, p, tb, wpot_per_bond) &
    !$omp& shared(tb_at, tb_rho, tb_e) &
    !$omp& reduction(+:w) reduction(+:error_loc)

    call tls_init(p%nat, vec=1)

    !$omp do
    I_loop: do I = 1, p%natloc

       I_is_el: if (IS_EL(tb%f, p, I)) then

          ni_loop: do ni = nl%seed(I), nl%last(I)
             J = GET_NEIGHBOR(nl, ni)

             J_is_el: if (IS_EL(tb%f, p, J)) then

                I_lt_J: if (I <= J) then

                   elI  = tb_at(I)%enr  ! Internal element number
                   elJ  = tb_at(J)%enr
                   noI  = tb_at(I)%no   ! Number of orbitals
                   noJ  = tb_at(J)%no

!                   write (*, *)  p%nat, I, J, elI, elJ

                   DIST(p, nl, I, ni, rIJ, abs_rIJ)
                   l    = rIJ / abs_rIJ

                   !
                   ! If we use multiprocessing, do only add the repulsive force once.
                   !
                   Ft = 0
                   if (abs_rIJ < db%R(elI, elJ)%cut) then
                      a = interval(db%R(elI, elJ), abs_rIJ, error_loc)
                      TRACE_DELAYED_ERROR_WITH_INFO("elI = " // elI // ", elJ = " // elJ // ", abs_rIJ = " // abs_rIJ, error_loc)
                      Ft = -l * df(db%R(elI, elJ), 1, abs_rIJ, a)
                   endif

                   !-------------------------------
                   !    the force from the
                   !  derivatives of S and H   
                   !-------------------------------
                   atoms_within_cutoff: if (abs_rIJ <= db%cut(elI, elJ)) then

                      !---------------------------------------
                      ! Create vector of spline inter-
                      ! polated coeffs at abs_rIJ. That is,
                      ! the matrix elements of H and
                      ! S, and their derivatives 
                      ! with respect to |abs_rIJ|. 
                      ! Do not calculate elements not needed.
                      !---------------------------------------
                      m  = max(noI, noJ)
                      nr = list(0, m)

                      a = interval(db%HS(elI, elJ), abs_rIJ, error_loc)
                      TRACE_DELAYED_ERROR_WITH_INFO("elI = " // elI // ", elJ = " // elJ // ", abs_rIJ = " // abs_rIJ, error_loc)
                      b = interval(db%HS(elJ, elI), abs_rIJ, error_loc)
                      TRACE_DELAYED_ERROR_WITH_INFO("elI = " // elI // ", elJ = " // elJ // ", abs_rIJ = " // abs_rIJ, error_loc)

                      dH_ij = 0.0_DP
                      dS_ij = 0.0_DP
                      do q = 1,nr
                         kk = list(q,m)
                         call f_and_df(db%HS(elI, elJ), kk,          abs_rIJ, a, dH_ij(kk), dH_ij(-kk))
                         call f_and_df(db%HS(elI, elJ), kk+MAX_NORB, abs_rIJ, a, dS_ij(kk), dS_ij(-kk))
                         call f_and_df(db%HS(elJ, elI), kk,          abs_rIJ, b, dH_ji(kk), dH_ji(-kk))
                         call f_and_df(db%HS(elJ, elI), kk+MAX_NORB, abs_rIJ, b, dS_ji(kk), dS_ji(-kk))
                      enddo

                      !---------------------------------
                      ! Use Slater-Koster transformation
                      ! rules to get the derivatives
                      ! with respect to all components
                      !---------------------------------
                      call mdiff(lo,noI,noJ,abs_rIJ,l,dH_ij,Hdiff_ij)
                      call mdiff(lo,noI,noJ,abs_rIJ,l,dS_ij,Sdiff_ij)
                      call mdiff(lo,noJ,noI,abs_rIJ,l,dH_ji,Hdiff_ji)
                      call mdiff(lo,noJ,noI,abs_rIJ,l,dS_ji,Sdiff_ji)

                      !------------------------------------------------------------
                      !  The force vector to orbital a at atom j due to 
                      !  orbital b at atom i, when an electrons occupy eigenstates
                      !  1:norb. The code is quite well optimized.
                      !------------------------------------------------------------

                      F2   = 0
                      k_loop: do k = 1, tb%nk
                         mu_loop: do mu0 = 1, noI
                            mu = get_orbital(noI, mu0)
                            nu_loop: do nu0 = 1, noJ
                               nu = get_orbital(noJ, nu0)

                               Imu = tb_at(I)%o1 + mu0 - 1 
                               Jnu = tb_at(J)%o1 + nu0 - 1
                      
                               rho_ImuJnu = tb_rho(Imu, Jnu, k)
                               e_ImuJnu   = tb_e(Imu, Jnu, k)

                               !------------------------------------------------------- 
                               ! if b>a (i.e. ang.momenta l_b>l_a), we must use
                               ! the other table
                               !------------------------------------------------------- 
                               if (nu > mu) then
                                  F2 = F2 + (rho_ImuJnu * Hdiff_ij(:, nu, mu) &
                                       - e_ImuJnu * Sdiff_ij(:, nu, mu))
                               else
                                  F2 = F2 + (rho_ImuJnu * Hdiff_ji(:, nu, mu) &
                                       - e_ImuJnu * Sdiff_ji(:, nu, mu))
                               endif

                            enddo nu_loop
                         enddo mu_loop
                      enddo k_loop

                      Ft = Ft - 2*F2
                   endif atoms_within_cutoff

                   !
                   ! Note: The I == J term does not contribute to the forces
                   ! but does give a significant contribution to the stress tensor!
                   ! 

                   if (J > p%natloc) then
                      VEC3(tls_vec1, I) = VEC3(tls_vec1, I) + 0.5_DP*Ft
                      VEC3(tls_vec1, J) = VEC3(tls_vec1, J) - 0.5_DP*Ft

                      wij  = -0.5_DP*outer_product(real(Ft, DP), rIJ)
                   else if (I /= J) then
                      VEC3(tls_vec1, I) = VEC3(tls_vec1, I) + Ft
                      VEC3(tls_vec1, J) = VEC3(tls_vec1, J) - Ft

                      wij  = -outer_product(real(Ft, DP), rIJ)
                   else
                      wij  = -0.5_DP*outer_product(real(Ft, DP), rIJ)
                   endif

                   w  = w + wij

                   if (present(wpot_per_bond)) then
                      SUM_VIRIAL(wpot_per_bond, ni, wij)
                   endif

                endif I_lt_J

             endif J_is_el

          enddo ni_loop

       endif I_is_el

    enddo I_loop

    call tls_reduce(p%nat, vec1=f)

    !$omp end parallel

    INVOKE_DELAYED_ERROR(error_loc, error)

    wpot = wpot + w

    call timer_stop("forces")

  endsubroutine forces



  !******************************************************************
  !
  !                    function mdiff
  !
  !  calculates the derivative of the H or S matrix element 
  !  <a|H|b> at dr (with direction cosine co) in direction i.
  !  R_vec = r*co(:); d = d (<a|H|b>(R_vec)) / d x_i
  !
  !******************************************************************
  subroutine mdiff(lo,noi,noj,r,co,ci,diff) 
    implicit none

    integer, intent(in)    :: lo(9),noi,noj
    real(DP), intent(in)   :: r,co(3),ci(-10:10)
    real(DP), intent(out)  :: diff(3,9,9)

    ! ---

!    real(DP), parameter  :: s3  = sqrt(3.0_DP)
    real(DP)  :: s3

    real(DP)  :: l,m,n,ll,mm,nn,li,mi,ni,g,d,lli,mmi,nni,c(-10:10)
    integer   :: a,b,i,mx

    real(DP)  :: dds, ddp, ddd, pds, pdp, pps, ppp, sds, sps, sss
    real(DP)  :: ddsi,ddpi,dddi,pdsi,pdpi,ppsi,pppi,sdsi,spsi,sssi

!    call timer('mdiff',1)
    s3  = sqrt(3.0_DP)

    diff = 0.0_DP

    !-------------------------------------------
    ! initial arrangements...
    ! using mx this way is an easy way to speed
    ! up if most of the atoms are e.g. carbon
    !-------------------------------------------
    l   = co(1)
    m   = co(2)
    n   = co(3)
    ll  = l**2
    mm  = m**2
    nn  = n**2
    dds = ci(1); ddp=ci(2); ddd=ci(3); pds=ci(4); pdp=ci(5)
    pps = ci(6); ppp=ci(7); sds=ci(8); sps=ci(9); sss=ci(10)
    mx  = max(get_orbital(noi, noi), get_orbital(noj, noj))

    !---------------------------------------------------------
    !e.g. spsi is the derivative of (sps) with respect to x_i
    !     li         - " -          l        - " -
    !     lli        - " -          l^2      - " -
    !     etc...
    !---------------------------------------------------------

    i_loop: do i=1,3
       c(-10:-1) = ci(-10:-1) * co(i)
       ddsi = c(-1); ddpi=c(-2); dddi=c(-3); pdsi=c(-4); pdpi=c(-5 )
       ppsi = c(-6); pppi=c(-7); sdsi=c(-8); spsi=c(-9); sssi=c(-10)

       li = ( k_delta(i,1) - l*co(i) )/r
       mi = ( k_delta(i,2) - m*co(i) )/r
       ni = ( k_delta(i,3) - n*co(i) )/r
       lli = 2*l*li
       mmi = 2*m*mi
       nni = 2*n*ni


       a_loop: do a = 1, mx
          b_loop: do b = a, mx

             !-------------------------------------------
             !
             ! choose the transformation rule according
             ! to the first orbital (a) and the second
             ! orbital (b) (now a<=b). This is the table
             ! of Slater & Koster, but filled with
             ! the missing rules (by permuting
             ! coordinates and dir. cosines, indicated
             ! by star) and derivated with respect
             ! to r_i
             !
             !-------------------------------------------
             select case(a)
             case(1)
                select case(b)
                case(1)
                   d = sssi
                   g = 0
                case(2)
                   d = l*spsi
                   g = li*sps
                case(3)
                   d = m*spsi
                   g = mi*sps
                case(4)
                   d = n*spsi
                   g = ni*sps
                case(5)
                   d = s3*l*m*sdsi
                   g = s3*(li*m+l*mi)*sds
                case(6) !(*)
                   d = s3*m*n*sdsi
                   g = s3*(mi*n+m*ni)*sds
                case(7) !(*)
                   d = s3*n*l*sdsi
                   g = s3*(ni*l+n*li)*sds
                case(8)
                   d = 0.5_DP*s3*(ll-mm)*sdsi
                   g = 0.5_DP*s3*(lli-mmi)*sds
                case(9)
                   d = (nn-0.5_DP*(ll+mm))*sdsi
                   g = (nni-0.5_DP*(lli+mmi))*sds
                case default
                   stop 'transf. not defined'
                end select
             case(2)
                select case(b)
                case(2)
                   d =                       ll *ppsi                     + (1-ll)*pppi
                   g =                       lli*pps                      + (-lli)*ppp
                case(3)
                   d =                       l*m*ppsi                -         l*m*pppi
                   g =               (li*m+l*mi)*pps                 - (li*m+l*mi)*ppp
                case(4)
                   d =                       l*n*ppsi                -         l*n*pppi
                   g =               (li*n+l*ni)*pps                 - (li*n+l*ni)*ppp
                case(5)
                   d =                   s3*ll*m*pdsi +                 m*(1-2*ll)*pdpi
                   g =          s3*(lli*m+ll*mi)*pds  + ( mi*(1-2*ll)+m*(-2*lli) )*pdp
                case(6)
                   d =                  s3*l*m*n*pdsi   -                  2*l*m*n*pdpi
                   g = s3*(li*m*n+l*mi*n+l*m*ni)*pds    - 2*(li*m*n+l*mi*n+l*m*ni)*pdp
                case(7)
                   d =                   s3*ll*n*pdsi +                 n*(1-2*ll)*pdpi
                   g =          s3*(lli*n+ll*ni)*pds  + ( ni*(1-2*ll)+n*(-2*lli) )*pdp
                case(8)
                   d =             0.5_DP*s3                 *l*(ll-mm)*pdsi +                   l*(1-ll+mm)*pdpi
                   g =             0.5_DP*s3*( li*(ll-mm)+l*(lli-mmi) )*pds  + ( li*(1-ll+mm)+l*(-lli+mmi) )*pdp
                case(9)
                   d =                            l*(nn-0.5_DP*(ll+mm))*pdsi              -          s3*l*nn*pdpi
                   g = ( li*(nn-0.5_DP*(ll+mm))+l*(nni-0.5_DP*(lli+mmi)) )*pds               - s3*(li*nn+l*nni)*pdp
                case default
                   stop 'transf. not defined'
                end select
             case(3)
                select case(b)
                case(3) !(*)
                   d =                                mm*ppsi                       + (1-mm)*pppi
                   g =                               mmi*pps                        + (-mmi)*ppp
                case(4) !(*)
                   d =                               m*n*ppsi                          - m*n*pppi
                   g =                       (mi*n+m*ni)*pps                   - (mi*n+m*ni)*ppp
                case(5) !(*)
                   d =                           s3*mm*l*pdsi                   + l*(1-2*mm)*pdpi
                   g =                  s3*(mmi*l+mm*li)*pds    + ( li*(1-2*mm)+l*(-2*mmi) )*pdp
                case(6) !(*)
                   d =                           s3*mm*n*pdsi                   + n*(1-2*mm)*pdpi
                   g =                  s3*(mmi*n+mm*ni)*pds    + ( ni*(1-2*mm)+n*(-2*mmi) )*pdp
                case(7) !(*)
                   d =                          s3*m*n*l*pdsi                      - 2*m*n*l*pdpi
                   g =         s3*(mi*n*l+m*ni*l+m*n*li)*pds      - 2*(mi*n*l+m*ni*l+m*n*li)*pdp
                case(8)
                   d =                  0.5_DP*s3*m*(ll-mm)*pdsi                  - m*(1+ll-mm)*pdpi
                   g = 0.5_DP*s3*( mi*(ll-mm)+m*(lli-mmi) )*pds  - ( mi*(1+ll-mm)+m*(lli-mmi) )*pdp
                case(9)
                   d =                m*(nn-0.5_DP*(ll+mm))*pdsi                      - s3*m*nn*pdpi
                   g = ( mi*(nn-0.5_DP*(ll+mm))+m*(nni-0.5_DP*(lli+mmi)) )*pds   -s3*(mi*nn+m*nni)*pdp
                case default
                   stop 'transf. not defined'
                end select
             case(4)
                select case(b)
                case(4) !(*)
                   d =                                nn*ppsi                       + (1-nn)*pppi
                   g =                               nni*pps                        + (-nni)*ppp
                case(5) !(*)
                   d =                          s3*l*m*n*pdsi                      - 2*m*n*l*pdpi
                   g =         s3*(li*m*n+l*mi*n+l*m*ni)*pds      - 2*(mi*n*l+m*ni*l+m*n*li)*pdp
                case(6) !(*)
                   d =                           s3*nn*m*pdsi                   + m*(1-2*nn)*pdpi
                   g =                  s3*(nni*m+nn*mi)*pds    + ( mi*(1-2*nn)+m*(-2*nni) )*pdp
                case(7) !(*)
                   d =                           s3*nn*l*pdsi                   + l*(1-2*nn)*pdpi
                   g =                  s3*(nni*l+nn*li)*pds    + ( li*(1-2*nn)+l*(-2*nni) )*pdp
                case(8)
                   d =                  0.5_DP*s3*n*(ll-mm)*pdsi                    - n*(ll-mm)*pdpi
                   g = 0.5_DP*s3*( ni*(ll-mm)+n*(lli-mmi) )*pds    - ( ni*(ll-mm)+n*(lli-mmi) )*pdp
                case(9)
                   d =                n*(nn-0.5_DP*(ll+mm))*pdsi                 + s3*n*(ll+mm)*pdpi
                   g = ( ni*(nn-0.5_DP*(ll+mm))+n*(nni-0.5_DP*(lli+mmi)) )*pds + s3*( ni*(ll+mm)+n*(lli+mmi) )*pdp
                case default
                   stop 'transf. not defined'
                end select
             case(5)
                select case(b)
                case(5)
                   d =           3*ll*mm*ddsi             + (ll+mm-4*ll*mm)*ddpi           + (nn+ll*mm)*dddi
                   g = 3*(lli*mm+ll*mmi)*dds + (lli+mmi-4*(lli*mm+ll*mmi) )*ddp + (nni+(lli*mm+ll*mmi))*ddd
                case(6)
                   d = 3*l*mm*n*ddsi + l*n*(1-4*mm)*ddpi + l*n*(mm-1)*dddi
                   g =                     3*(li*mm*n+l*mmi*n+l*mm*ni)*dds &
                        + ( li*n*(1-4*mm)+l*ni*(1-4*mm)+l*n*(-4*mmi) )*ddp &
                        +        ( li*n*(mm-1)+l*ni*(mm-1)+l*n*(mmi) )*ddd
                case(7)
                   d = 3*ll*m*n*ddsi + m*n*(1-4*ll)*ddpi + m*n*(ll-1)*dddi
                   g =                     3*(lli*m*n+ll*mi*n+ll*m*ni)*dds &
                        + ( mi*n*(1-4*ll)+m*ni*(1-4*ll)+m*n*(-4*lli) )*ddp &
                        + ( mi*n*(ll-1)+m*ni*(ll-1)+m*n*(lli) )*ddd
                case(8)
                   d = 1.5_DP*l*m*(ll-mm)*ddsi + 2*l*m*(mm-ll)*ddpi + 0.5_DP*l*m*(ll-mm)*dddi
                   g = 1.5_DP*( li*m*(ll-mm)+l*mi*(ll-mm)+l*m*(lli-mmi) )*dds &
                        + 2*( li*m*(mm-ll)+l*mi*(mm-ll)+l*m*(mmi-lli) )*ddp &
                        + 0.5_DP*( li*m*(ll-mm)+l*mi*(ll-mm)+l*m*(lli-mmi) )*ddd
                case(9)
                   d = s3*l*m*(nn-0.5_DP*(ll+mm))*ddsi - 2*s3*l*m*nn*ddpi + 0.5_DP*s3*l*m*(1+nn)*dddi
                   g = s3*( li*m*(nn-0.5_DP*(ll+mm))+l*mi*(nn-0.5_DP*(ll+mm))+l*m*(nni-0.5_DP*(lli+mmi)) )*dds &
                        -2*s3*(li*m*nn+l*mi*nn+l*m*nni)*ddp &
                        + 0.5_DP*s3*( li*m*(1+nn)+l*mi*(1+nn)+l*m*(nni) )*ddd
                case default
                   stop 'transf. not defined'
                end select
             case(6)
                select case(b)
                case(6) !(*)
                   d = 3*mm*nn*ddsi + (mm+nn-4*mm*nn)*ddpi + (ll+mm*nn)*dddi
                   g =             3*(mmi*nn+mm*nni)*dds &
                        + (mmi+nni-4*(mmi*nn+mm*nni))*ddp &
                        +         (lli+mmi*nn+mm*nni)*ddd
                case(7) !(*)
                   d = 3*m*nn*l*ddsi + m*l*(1-4*nn)*ddpi + m*l*(nn-1)*dddi
                   g =                    3*(mi*nn*l+m*nni*l+m*nn*li)*dds &
                        + ( mi*l*(1-4*nn)+m*li*(1-4*nn)+m*l*(-4*nni) )*ddp &
                        + ( mi*l*(nn-1)+m*li*(nn-1)+m*l*(nni) )*ddd
                case(8)
                   d = 1.5_DP*m*n*(ll-mm)*ddsi - m*n*(1+2*(ll-mm))*ddpi + m*n*(1+0.5_DP*(ll-mm))*dddi
                   g =                     1.5_DP*( mi*n*(ll-mm)+m*ni*(ll-mm)+m*n*(lli-mmi) )*dds &
                        - ( mi*n*(1+2*(ll-mm))+m*ni*(1+2*(ll-mm))+m*n*(2*lli-2*mmi) )*ddp &
                        + ( mi*n*(1+0.5_DP*(ll-mm))+m*ni*(1+0.5_DP*(ll-mm))+m*n*(0.5_DP*(lli-mmi)) )*ddd
                case(9)
                   d = s3*m*n*(nn-0.5_DP*(ll+mm))*ddsi + s3*m*n*(ll+mm-nn)*ddpi - 0.5_DP*s3*m*n*(ll+mm)*dddi
                   g = s3*( mi*n*(nn-0.5_DP*(ll+mm)) + m*ni*(nn-0.5_DP*(ll+mm))+m*n*(nni-0.5_DP*(lli+mmi)) ) * dds &
                        + s3*( mi*n*(ll+mm-nn)+m*ni*(ll+mm-nn)+m*n*(lli+mmi-nni) )* ddp &
                        - 0.5_DP*s3*( mi*n*(ll+mm)+m*ni*(ll+mm)+m*n*(lli+mmi) )* ddd
                case default
                   stop 'transf. not defined'
                end select
             case(7)
                select case(b)
                case(7) !(*)
                   d = 3*nn*ll*ddsi + (nn+ll-4*nn*ll)*ddpi + (mm+nn*ll)*dddi
                   g =               3*(nni*ll+nn*lli)*dds &
                        + ( nni+lli-4*(nni*ll+nn*lli) )*ddp &
                        + (mmi+nni*ll+nn*lli)*ddd
                case(8)
                   d = 1.5_DP*n*l*(ll-mm)*ddsi + n*l*(1-2*(ll-mm))*ddpi - n*l*(1-0.5_DP*(ll-mm))*dddi
                   g =                      1.5_DP*( ni*l*(ll-mm)+n*li*(ll-mm)+n*l*(lli-mmi) )*dds &
                        + ( ni*l*(1-2*(ll-mm))+n*li*(1-2*(ll-mm))+n*l*(-2*(lli-mmi)) )*ddp &
                        - ( ni*l*(1-0.5_DP*(ll-mm))+n*li*(1-0.5_DP*(ll-mm))+n*l*(-0.5_DP*(lli-mmi)) )*ddd
                case(9)
                   d =  s3*l*n*(nn-0.5_DP*(ll+mm))*ddsi + s3*l*n*(ll+mm-nn)*ddpi - 0.5_DP*s3*l*n*(ll+mm)*dddi
                   g = s3*( li*n*(nn-0.5_DP*(ll+mm))+l*ni*(nn-0.5_DP*(ll+mm))+l*n*(nni-0.5_DP*(lli+mmi)) ) *dds &
                        + s3*( li*n*(ll+mm-nn)+l*ni*(ll+mm-nn)+l*n*(lli+mmi-nni) )*ddp &
                        - 0.5_DP*s3*( li*n*(ll+mm)+l*ni*(ll+mm)+l*n*(lli+mmi) )*ddd
                case default
                   stop 'transf. not defined'
                end select
             case(8)
                select case(b)
                case(8)
                   d = 0.75_DP*(ll-mm)**2*ddsi + (ll+mm-(ll-mm)**2)*ddpi + (nn+0.25_DP*(ll-mm)**2)*dddi
                   g =         0.75_DP*2*(ll-mm)*(lli-mmi)*dds &
                        + (lli+mmi-2*(ll-mm)*(lli-mmi))*ddp &
                        + (nni+0.25_DP*2*(ll-mm)*(lli-mmi))*ddd
                case(9)
                   d = 0.5_DP*s3*(ll-mm)*(nn-0.5_DP*(ll+mm))*ddsi + s3*nn*(mm-ll)*ddpi + 0.25_DP*s3*(1+nn)*(ll-mm)*dddi
                   g = 0.5_DP*s3*( (lli-mmi)*(nn-0.5_DP*(ll+mm))+(ll-mm)*(nni-0.5_DP*(lli+mmi)) )*dds &
                        + s3*( nni*(mm-ll)+nn*(mmi-lli) )*ddp &
                        + 0.25_DP*s3*( nni*(ll-mm)+(1+nn)*(lli-mmi) )*ddd
                case default
                   stop 'transf. not defined'
                end select
             case(9)
                select case(b)
                case(9)
                   d = (nn-0.5_DP*(ll+mm))**2*ddsi + 3*nn*(ll+mm)*ddpi + 0.75_DP*(ll+mm)**2*dddi
                   g =        2*(nn-0.5_DP*(ll+mm))*(nni-0.5_DP*(lli+mmi))*dds &
                        + 3*( nni*(ll+mm)+nn*(lli+mmi) )*ddp &
                        + 0.75_DP*2*(ll+mm)*(lli+mmi)*ddd
                case default
                   stop 'transf. not defined'
                end select
             case default
                stop 'transf. not defined'
             end select

             !-------------------------------------------------------------
             !
             ! d is the part where (mel)'s are derivated and g is the part
             ! where the geometric factor is derivated (chain rule)
             !
             !
             ! The pre-factor comes from the interchange of the orbitals
             ! with respect to tabulated values
             ! [parity = (-1)^(angular momentum)]
             !
             !--------------------------------------------------------------
             diff(i,a,b) = (d+g)
             diff(i,b,a) = (d+g)*(-1)**( lo(a)+lo(b) )

          end do b_loop
       end do a_loop
    end do i_loop
!    call timer('mdiff',2)
  endsubroutine mdiff

endmodule dense_force_calc


