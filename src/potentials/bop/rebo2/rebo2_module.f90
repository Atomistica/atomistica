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
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine INIT_FUNC(this)
    implicit none

    type(BOP_TYPE), intent(inout)   :: this

    ! ---
    
    if (this%zero_tables) then
        this%in_Fcc = 0.0_DP
        this%in_dFdi = 0.0_DP
        this%in_dFdj = 0.0_DP
        this%in_dFdk = 0.0_DP

        this%in_Fch = 0.0_DP
        this%in_Fhh = 0.0_DP

        this%in_Pcc = 0.0_DP
        this%in_Pch = 0.0_DP
        this%in_Tcc = 0.0_DP
    endif

  endsubroutine INIT_FUNC


  !>
  !! Destructor
  !<
  subroutine DEL_FUNC(this)
    implicit none

    type(BOP_TYPE), intent(inout)            :: this

    ! ---

    call rebo2_db_del(this)

    if (allocated(this%internal_el)) then
       deallocate(this%internal_el)
    endif

  endsubroutine DEL_FUNC


  subroutine BIND_TO_FUNC(this, p, nl, ierror)
    implicit none

    type(BOP_TYPE),    intent(inout) :: this
    type(particles_t), intent(inout) :: p
    type(neighbors_t), intent(inout) :: nl
    integer, optional, intent(inout) :: ierror

    ! ---

    integer   :: i, j
    real(DP)  :: c_cc, c_ch, c_hh

    ! ---

    this%els = filter_from_string(this%elements, p, ierror)
    PASS_ERROR(ierror)

    call rebo2_db_init_with_parameters( &
         this, &
         this%in_Fcc, this%in_dFdi, this%in_dFdj, this%in_dFdk, &
         this%in_Fch, this%in_Fhh, &
         this%in_Pcc, this%in_Pch, &
         this%in_Tcc)

!    call rebo2_db_init(this)

#ifdef SCREENING
    c_cc = sqrt(this%C_dr_cut)*maxval( [ this%cc_in_r2, this%cc_ar_r2, &
         this%cc_bo_r2, this%cc_nc_r2 ] )
    c_ch = this%ch_r2
    c_hh = this%hh_r2
#else
    c_cc = this%cc_in_r2
    c_ch = this%ch_r2
    c_hh = this%hh_r2
#endif

    do i = 1, p%nel
       do j = i, p%nel
          if (p%el2Z(i) == C_ .and. p%el2Z(j) == C_) then
             call request_interaction_range(nl, c_cc, i, j)
          else if ( &
               (p%el2Z(i) == C_ .and. p%el2Z(j) == H_) .or. &
               (p%el2Z(i) == H_ .and. p%el2Z(j) == C_) &
               ) then
             call request_interaction_range(nl, c_ch, i, j)
          else if (p%el2Z(i) == H_ .and. p%el2Z(j) == H_) then
             call request_interaction_range(nl, c_hh, i, j)
          endif
       enddo
    enddo

#ifdef SCREENING
    c_cc = (2+sqrt(this%C_dr_cut))*maxval( [ this%cc_in_r2, &
         this%cc_ar_r2, this%cc_bo_r2, this%cc_nc_r2 ] )
#else
    c_cc = 3*this%cc_in_r2
#endif
    call request_border(p, c_cc)

    if (allocated(this%internal_el))  deallocate(this%internal_el)

    allocate(this%internal_el(p%maxnatloc))

  endsubroutine BIND_TO_FUNC


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine COMPUTE_FUNC(this, p, nl, epot, f, wpot, epot_per_at, &
       epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)
    implicit none

    type(BOP_TYPE), intent(inout)      :: this
    type(particles_t), intent(inout)   :: p
    type(neighbors_t), intent(inout)   :: nl
    real(DP), intent(inout)            :: epot
    real(DP), intent(inout)            :: f(3, p%maxnatloc)  !< forces
    real(DP), intent(inout)            :: wpot(3, 3)
    real(DP), intent(inout), optional  :: epot_per_at(p%maxnatloc)
    real(DP), intent(inout), optional  :: epot_per_bond(nl%neighbors_size)
    real(DP), intent(inout), optional  :: f_per_bond(3, nl%neighbors_size)
#ifdef LAMMPS
    real(DP), intent(inout), optional  :: wpot_per_at(6, p%maxnatloc)
    real(DP), intent(inout), optional  :: wpot_per_bond(6, nl%neighbors_size)
#else
    real(DP), intent(inout), optional  :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), intent(inout), optional  :: wpot_per_bond(3, 3, nl%neighbors_size)
#endif
    integer, intent(inout), optional   :: ierror

    ! ---

    integer :: i, d, nebmax, nebavg

    ! ---

    call timer_start(BOP_NAME_STR // "_force")

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    if (size(this%internal_el) < p%maxnatloc) then
       deallocate(this%internal_el)
       allocate(this%internal_el(p%maxnatloc))
    endif

    this%internal_el = 0
    nebmax = 0
    nebavg = 0
    do i = 1, p%nat
       if (IS_EL(this%els, p, i)) then
          if (p%el2Z(p%el(i)) == C_) then
             this%internal_el(i) = rebo2_C_
          else if (p%el2Z(p%el(i)) == H_) then
             this%internal_el(i) = rebo2_H_
          endif
       endif
       d = nl%last(i)-nl%seed(i)+1
       nebmax = max(nebmax, d)
       nebavg = nebavg + d
    enddo
    nebavg = (nebavg+1)/max(p%nat, 1)+1

#ifdef LAMMPS
    call BOP_KERNEL( &
         this, &
         p%maxnatloc, p%natloc, p%nat, p%r_non_cyc, p%tag, this%internal_el, &
         nebmax, nebavg, nl%seed, nl%last, nl%neighbors, nl%neighbors_size, &
         epot, f, wpot, &
         epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
         ierror)
#else
    call BOP_KERNEL( &
         this, p%Abox, &
         p%maxnatloc, p%natloc, p%nat, p%r_non_cyc, this%internal_el, &
         nebmax, nebavg, nl%seed, nl%last, nl%neighbors, nl%neighbors_size, &
         nl%dc, p%shear_dx, &
         epot, f, wpot, &
         epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
         ierror)
#endif
    PASS_ERROR(ierror)

    call timer_stop(BOP_NAME_STR // "_force")

  endsubroutine COMPUTE_FUNC

