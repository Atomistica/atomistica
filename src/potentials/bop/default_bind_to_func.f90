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
  !! Notify potential of particles, neighbors objects to use in the future
  !<
  subroutine BIND_TO_FUNC(this, p, nl, ierror)
    implicit none
 
    type(BOP_TYPE),    intent(inout) :: this
    type(particles_t), intent(inout) :: p
    type(neighbors_t), intent(inout) :: nl
    integer, optional, intent(out)   :: ierror

    ! ---

    integer          :: i, j, ii, jj, nel, npairs, Z

    real(DP)         :: x(this%db%nel*(this%db%nel+1)/2), cutoff

    ! ---

    INIT_ERROR(ierror)

    call del(this)

#ifdef SCREENING

    this%Cmin      = this%db%Cmin
    this%Cmax      = this%db%Cmax
    this%dC        = this%Cmax-this%Cmin

    !
    ! The maximum cutoff needs to be the maximum distance and atom can be away
    ! and still considered for screening.
    !
    ! This means there is a scale factor for the distance a screening neighbor
    ! can. It is given by
    !   X = (xik/xij)^2 = C^2/(4*(C-1))
    ! where xij is the bond distance and xik the distance to the screening
    ! neighbor.
    !
    ! Note that at C = 2 the maximum distance is the xik^2 = xij^2 and hence
    ! C_dr_cut = 1.0_DP below. For C < 2 we also need to consider at least
    ! xik^2 = xij^2.
    !

    this%C_dr_cut  = 1.0_DP
    where (this%Cmax > 2.0_DP)
       this%C_dr_cut = this%Cmax**2/(4*(this%Cmax-1))
    endwhere

#endif

    nel    = this%db%nel
    npairs = nel*(nel+1)/2

    this%Z2db = -1

    do i = 1, this%db%nel
       Z = atomic_number(a2s(this%db%el(:, i)))
       if (Z > 0) then
          this%Z2db(Z)  = i
       else
          RAISE_ERROR("Unknown element '" // trim(a2s(this%db%el(:, i))) // "'.", ierror)
       endif
    enddo

    do i = 1, npairs
       call init(this%cut_in(i), this%db%r1(i), this%db%r2(i))

       this%cut_in_l(i)   = this%db%r1(i)
       this%cut_in_h(i)   = this%db%r2(i)
       this%cut_in_h2(i)  = this%db%r2(i)**2

#ifdef SCREENING
       call init(this%cut_out(i), this%db%or1(i), this%db%or2(i))

       this%cut_out_l(i)  = this%db%or1(i)
       this%cut_out_h(i)  = this%db%or2(i)

       call init(this%cut_bo(i), this%db%bor1(i), this%db%bor2(i))

       this%cut_bo_l(i)   = this%db%bor1(i)
       this%cut_bo_h(i)   = this%db%bor2(i)

       this%max_cut_sq(i)   = max( &
            this%cut_in_h(i), &
            this%cut_out_h(i), &
            this%cut_bo_h(i) &
            )**2
#endif
    enddo

    !
    ! Request interaction range for each element pair
    !

#ifdef SCREENING

    x = sqrt(this%C_dr_cut(1:npairs))

#endif

    do i = 1, p%nel
       if (p%el2Z(i) > 0) then
          do j = 1, p%nel
             if (p%el2Z(j) > 0) then
                ii = this%Z2db(p%el2Z(i))
                jj = this%Z2db(p%el2Z(j))
                if (ii > 0 .and. jj > 0) then
                   nel = Z2pair(this, ii, jj)
#ifdef SCREENING
                   cutoff = x(nel)*sqrt(this%max_cut_sq(nel))
#else
                   cutoff = this%cut_in_h(nel)
#endif
                   call request_interaction_range(nl, cutoff, i, j)
#ifdef LAMMPS
                   call set_interaction_range(p, 2*cutoff, i, j)
#endif
                endif
             endif
          enddo
       endif
    enddo

  endsubroutine BIND_TO_FUNC
