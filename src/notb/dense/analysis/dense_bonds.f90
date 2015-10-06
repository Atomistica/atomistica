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
!! Determine the covalent bond energies
!! See: Bornsen et al., J.Phys.: Cond. Mat. 11, L287 (1999)
!<

#include "macros.inc"

module dense_bonds
  use supplib

  use particles
  use neighbors

  use materials

  use dense_hamiltonian

contains

  !>
  !! Determine covalent energies
  !<
  subroutine bond_analysis(tb, p, nl, overlap_population, Loewdin_bond_order, &
                           e_cov)
    implicit none

    type(dense_hamiltonian_t), intent(in)  :: tb
    type(particles_t),         intent(in)  :: p
    type(neighbors_t),         intent(in)  :: nl
    real(DP),        optional, intent(out) :: overlap_population(nl%neighbors_size)
    real(DP),        optional, intent(out) :: Loewdin_bond_order(nl%neighbors_size)
    real(DP),        optional, intent(out) :: e_cov(nl%neighbors_size)

    ! ---

    integer :: i, ni, j, k, a, b, ia, jb

    real(DP) :: overlap_population_accum, Loewdin_bond_order_accum, e_cov_accum
    real(DP) :: alpha, beta

    WF_T(DP),             pointer :: rho(:, :, :), H(:, :, :), S(:, :, :)
    type(notb_element_t), pointer :: at(:)

    WF_T(DP), allocatable :: sqrt_S(:, :), Loewdin_rho(:, :, :)

    ! ---

    call timer_start("bond_analysis")

    call c_f_pointer(tb%rho, rho, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%H, H, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%S, S, [tb%norb, tb%norb, tb%nk])
    call c_f_pointer(tb%at, at, [p%nat])

    if (present(Loewdin_bond_order)) then
       allocate(sqrt_S(tb%norb, tb%norb))
       allocate(Loewdin_rho(tb%norb, tb%norb, tb%nk))
       do k = 1, tb%nk
          sqrt_S = sqrtm(S(:, :, k))
          Loewdin_rho(:, :, k) = matmul(sqrt_S, matmul(rho(:, :, k), sqrt_S))
       enddo
    endif

    i_loop: do i = 1, p%nat
       ni_loop: do ni = nl%seed(i), nl%last(i)
          j = nl%neighbors(ni)

          overlap_population_accum = 0.0_DP
          Loewdin_bond_order_accum = 0.0_DP
          e_cov_accum = 0.0_DP

          a_loop: do a = 1, at(i)%no
             ia = at(i)%o1 + a - 1
             b_loop: do b = 1, at(j)%no
                jb = at(j)%o1 + b - 1

                kpoint_loop: do k = 1, tb%nk

                   overlap_population_accum = overlap_population_accum + &
                                              rho(ia, jb, k) * S(jb, ia, k)

                   if (present(Loewdin_bond_order)) then
                      Loewdin_bond_order_accum = Loewdin_bond_order_accum + &
                                                 Loewdin_rho(ia, jb, k) + &
                                                 Loewdin_rho(jb, ia, k)
                   endif

                   ! Note: For SCC NOTB there is a contribution from phi
                   !   H(jb, ia) -= 0.5_DP*S(jb, ia)*(phi(i) + phi(j))
                   ! but this shift due to the electrostatic potential cancels
                   ! in the expression below.
                   e_cov_accum = e_cov_accum + rho(ia, jb, k) &
                        * ( H(jb, ia, k) - &
                            0.5_DP*S(jb, ia, k)*(H(ia, ia, k) + H(jb, jb, k)) &
                          )

                enddo kpoint_loop
                   
             enddo b_loop
          enddo a_loop

          if (present(overlap_population)) then
             overlap_population(ni) = overlap_population_accum
          endif
          if (present(Loewdin_bond_order)) then
             Loewdin_bond_order(ni) = 0.5_DP*Loewdin_bond_order_accum
          endif
          if (present(e_cov)) then
             e_cov(ni) = e_cov_accum
          endif
       enddo ni_loop
    enddo i_loop

    if (allocated(sqrt_S))       deallocate(sqrt_S)
    if (allocated(Loewdin_rho))  deallocate(Loewdin_rho)

    call timer_stop("bond_analysis")

  endsubroutine bond_analysis

endmodule dense_bonds
