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
  subroutine INIT_FUNC(this, &
       db, &
       el, D0, r0, S, beta, gamma, c, d, h, n, m, alpha, omega, r1, r2 &
#ifdef SCREENING
       , or1, or2, bor1, bor2, Cmin, Cmax &
#endif
       , ierror)

    implicit none

    type(BOP_TYPE), intent(inout)            :: this
    type(BOP_DB_TYPE), intent(in), optional  :: db
    character(2), intent(in), optional       :: el(:)
    real(DP), intent(in),    optional        :: D0(:)
    real(DP), intent(in),    optional        :: r0(:)
    real(DP), intent(in),    optional        :: S(:)
    real(DP), intent(in),    optional        :: beta(:)
    real(DP), intent(in),    optional        :: gamma(:)
    real(DP), intent(in),    optional        :: c(:)
    real(DP), intent(in),    optional        :: d(:)
    real(DP), intent(in),    optional        :: h(:)
    real(DP), intent(in),    optional        :: n(:)
    real(DP), intent(in),    optional        :: alpha(:)
    real(DP), intent(in),    optional        :: omega(:)
    integer,  intent(in),    optional        :: m(:)
    real(DP), intent(in),    optional        :: r1(:)
    real(DP), intent(in),    optional        :: r2(:)
#ifdef SCREENING
    real(DP), intent(in),    optional        :: or1(:)
    real(DP), intent(in),    optional        :: or2(:)
    real(DP), intent(in),    optional        :: bor1(:)
    real(DP), intent(in),    optional        :: bor2(:)
    real(DP), intent(in),    optional        :: Cmin(:)
    real(DP), intent(in),    optional        :: Cmax(:)
#endif
    integer,  intent(inout), optional        :: ierror

    ! ---

    integer :: i, j

    ! ---

    call prlog("- " // BOP_NAME_STR // " -")

    if (present(db)) then

       this%db = db

       call prlog("     Using database: " // trim(this%db%ref))

    else

       ! Is the reference string set?
       ! Then search for a parameter set by string.
       if (trim(this%ref) /= "*") then
          j = -1
          do i = 1, size(BOP_DB)
             if (trim(this%ref) == BOP_DB(i)%ref(1:len_trim(this%ref))) then
                if (j > 0) then
                   RAISE_ERROR("Reference string '" // this%ref // "' not unique. Matching entries: '" // BOP_DB(j)%ref // "' and '" // BOP_DB(i)%ref // "'.", ierror)
                endif
                j = i
             endif
          enddo

          if (j > 0) then
             this%db = BOP_DB(j)
          else
             RAISE_ERROR("Could not find parameter set for reference '" // this%ref // "' in database.", ierror)
          endif

          call prlog("     Using database: " // trim(this%db%ref))
       endif

       ASSIGN_STRING_ARRAY_PROPERTY(el,    this%db%el,    this%db%nel, i)
       ASSIGN_ARRAY_PROPERTY(D0,    this%db%D0,    this%db%nD0)
       ASSIGN_ARRAY_PROPERTY(r0,    this%db%r0,    this%db%nr0)
       ASSIGN_ARRAY_PROPERTY(S,     this%db%S,     this%db%nS)
       ASSIGN_ARRAY_PROPERTY(beta,  this%db%beta,  this%db%nbeta)
       ASSIGN_ARRAY_PROPERTY(gamma, this%db%gamma, this%db%ngamma)
       ASSIGN_ARRAY_PROPERTY(c,     this%db%c,     this%db%nc)
       ASSIGN_ARRAY_PROPERTY(d,     this%db%d,     this%db%nd)
       ASSIGN_ARRAY_PROPERTY(h,     this%db%h,     this%db%nh)
       ASSIGN_ARRAY_PROPERTY(n,     this%db%n,     this%db%nn)
       ASSIGN_ARRAY_PROPERTY(alpha, this%db%alpha, this%db%nalpha)
       ASSIGN_ARRAY_PROPERTY(omega, this%db%omega, this%db%nomega)
       ASSIGN_ARRAY_PROPERTY(m,     this%db%m,     this%db%nm)
       ASSIGN_ARRAY_PROPERTY(r1,    this%db%r1,    this%db%nr1)
       ASSIGN_ARRAY_PROPERTY(r2,    this%db%r2,    this%db%nr2)
#ifdef SCREENING
       ASSIGN_ARRAY_PROPERTY(or1,   this%db%or1,   this%db%nor1)
       ASSIGN_ARRAY_PROPERTY(or2,   this%db%or2,   this%db%nor2)
       ASSIGN_ARRAY_PROPERTY(bor1,  this%db%bor1,  this%db%nbor1)
       ASSIGN_ARRAY_PROPERTY(bor2,  this%db%bor2,  this%db%nbor2)
       ASSIGN_ARRAY_PROPERTY(Cmin,  this%db%Cmin,  this%db%nCmin)
       ASSIGN_ARRAY_PROPERTY(Cmax,  this%db%Cmax,  this%db%nCmax)
#endif

    endif

#ifdef SCREENING
    this%Cmin      = this%db%Cmin
    this%Cmax      = this%db%Cmax
    this%dC        = this%Cmax-this%Cmin
    this%C_dr_cut  = this%Cmax**2/(4*(this%Cmax-1))
#endif

!    call prlog("     el     = " // this%db%el(1:this%db%nel))
    call prlog("     D0     = " // this%db%D0(1:this%db%nD0))
    call prlog("     r0     = " // this%db%r0(1:this%db%nr0))
    call prlog("     S      = " // this%db%S(1:this%db%nS))
    call prlog("     beta   = " // this%db%beta(1:this%db%nbeta))
    call prlog("     gamma  = " // this%db%gamma(1:this%db%ngamma))
    call prlog("     c      = " // this%db%c(1:this%db%nc))
    call prlog("     d      = " // this%db%d(1:this%db%nd))
    call prlog("     h      = " // this%db%h(1:this%db%nh))
    call prlog("     n      = " // this%db%n(1:this%db%nn))
    call prlog("     alpha  = " // this%db%alpha(1:this%db%nalpha))
    call prlog("     omega  = " // this%db%omega(1:this%db%nomega))
    call prlog("     m      = " // this%db%m(1:this%db%nm))
    call prlog("     r1     = " // this%db%r1(1:this%db%nr1))
    call prlog("     r2     = " // this%db%r2(1:this%db%nr2))
#ifdef SCREENING
    call prlog("     or1    = " // this%db%or1(1:this%db%nor1))
    call prlog("     or2    = " // this%db%or2(1:this%db%nor2))
    call prlog("     bor1   = " // this%db%bor1(1:this%db%nbor1))
    call prlog("     bor2   = " // this%db%bor2(1:this%db%nbor2))
    call prlog("     Cmin   = " // this%db%Cmin(1:this%db%nCmin))
    call prlog("     Cmax   = " // this%db%Cmax(1:this%db%nCMax))
#endif

  endsubroutine INIT_FUNC


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine DEL_FUNC(this)
    implicit none

    type(BOP_TYPE), intent(inout)  :: this

    ! ---

    if (this%neighbor_list_allocated) then
       deallocate(this%neb)
       deallocate(this%nbb)
#ifndef LAMMPS
       deallocate(this%dcell)
#endif
       deallocate(this%bndtyp)
       deallocate(this%bndlen)
       deallocate(this%bndnm)
       deallocate(this%cutfcnar)
       deallocate(this%cutdrvar)

#ifdef SCREENING
       deallocate(this%cutfcnbo)
       deallocate(this%cutdrvbo)
       deallocate(this%sneb_seed)
       deallocate(this%sneb_last)
       deallocate(this%sneb)
       deallocate(this%sbnd)
       deallocate(this%sfacbo)
       deallocate(this%cutdrarik)
       deallocate(this%cutdrarjk)
       deallocate(this%cutdrboik)
       deallocate(this%cutdrbojk)
#endif

       this%neighbor_list_allocated = .false.
    endif

  endsubroutine DEL_FUNC


  !>
  !! Bind to 
  !!
  !! Set-up internal database, etc.
  !<
  subroutine BIND_TO_FUNC(this, p, nl, ierror)
    implicit none

    type(BOP_TYPE), intent(inout)     :: this
    type(particles_t), intent(inout)  :: p
    type(neighbors_t), intent(inout)  :: nl
    integer, intent(inout), optional  :: ierror

    ! ---

    integer          :: i, j, npairs, Z, ii, jj, nel

#ifdef SCREENING
    real(DP)         :: x(this%db%nel**2)
#endif

    ! ---

    npairs = this%db%nel**2

    if (npairs /= this%db%nD0 .or. &
        npairs /= this%db%nr0 .or. &
        npairs /= this%db%nS  .or. &
        npairs /= this%db%nbeta .or. &
        npairs /= this%db%ngamma .or. &
        npairs /= this%db%nc .or. &
        npairs /= this%db%nd .or. &
        npairs /= this%db%nh .or. &
#ifdef SCREENING
        npairs /= this%db%nor1 .or. &
        npairs /= this%db%nor2 .or. &
        npairs /= this%db%nbor1 .or. &
        npairs /= this%db%nbor2 .or. &
        npairs /= this%db%nCmin .or. &
        npairs /= this%db%nCmax .or. &
#endif
        npairs /= this%db%nr1 .or. &
        npairs /= this%db%nr2) then

       write (*, '(A,I2)')  "nel     = ", this%db%nel
       write (*, '(A,I2)')  "nD0     = ", this%db%nD0
       write (*, '(A,I2)')  "nr0     = ", this%db%nr0
       write (*, '(A,I2)')  "nS      = ", this%db%nS
       write (*, '(A,I2)')  "nbeta   = ", this%db%nbeta
       write (*, '(A,I2)')  "ngamma  = ", this%db%ngamma
       write (*, '(A,I2)')  "nc      = ", this%db%nc
       write (*, '(A,I2)')  "nh      = ", this%db%nh
       write (*, '(A,I2)')  "nr1     = ", this%db%nr1
       write (*, '(A,I2)')  "nr2     = ", this%db%nr2
#ifdef SCREENING
       write (*, '(A,I2)')  "nor1    = ", this%db%nor1
       write (*, '(A,I2)')  "nor2    = ", this%db%nor2
       write (*, '(A,I2)')  "nbor1   = ", this%db%nbor1
       write (*, '(A,I2)')  "nbor2   = ", this%db%nbor2
       write (*, '(A,I2)')  "nCmin   = ", this%db%nCmin
       write (*, '(A,I2)')  "nCmax   = ", this%db%nCmax
#endif

       RAISE_ERROR("The number of entries must be identical for all parameters.", ierror)
    endif


    do i = 1, this%db%nel
      do j = 1, this%db%nel
        if (this%db%r0(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL)) < 0.0_DP) then
          this%db%D0(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))    = this%db%D0(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%r0(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))    = this%db%r0(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%S(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))     = this%db%S(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%beta(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))  = this%db%beta(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%gamma(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL)) = this%db%gamma(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%c(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))     = this%db%c(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%d(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))     = this%db%d(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%h(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))     = this%db%h(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%n(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))     = this%db%n(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
#ifdef SCREENING
          this%db%r1(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))    = this%db%r1(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%r2(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))    = this%db%r2(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%or1(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))   = this%db%or1(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%or2(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))   = this%db%or2(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%bor1(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))  = this%db%bor1(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%bor2(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))  = this%db%bor2(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%Cmin(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))  = this%db%Cmin(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%Cmax(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))  = this%db%Cmax(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
#else
          this%db%r1(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))    = this%db%r1(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
          this%db%r2(PAIR_INDEX_NS(i, j, JUSLIN_MAX_EL))    = this%db%r2(PAIR_INDEX_NS(j, i, JUSLIN_MAX_EL))
#endif
        endif
      enddo
    enddo


#ifdef SCREENING
    this%Cmin      = this%db%Cmin
    this%Cmax      = this%db%Cmax
    this%dC        = this%Cmax-this%Cmin
    this%C_dr_cut  = this%Cmax**2/(4*(this%Cmax-1))
#endif


    this%Z2db  = 0

    do i = 1, this%db%nel
       Z = atomic_number(a2s(this%db%el(:, i)))
       if (Z > 0 .and. Z <= MAX_Z) then
          this%Z2db(Z)  = i
       else
          RAISE_ERROR("Unknown element '" // trim(a2s(this%db%el(:, i))) // "' with database index " // i // ".", ierror)
       endif
    enddo

    do i = 1, npairs
       this%bo_exp(i)   = - 0.5_DP / this%db%n(i)
       this%bo_fac(i)   = 0.5_DP * this%bo_exp(i) * this%db%n(i)
       this%bo_exp1(i)  = this%bo_exp(i) - 1.0_DP

       this%expR(i)  = this%db%beta(i)*sqrt(2*this%db%S(i))
       this%expA(i)  = this%db%beta(i)*sqrt(2/this%db%S(i))

       this%c_sq(i)  = this%db%c(i)**2
       this%d_sq(i)  = this%db%d(i)**2
       if (this%d_sq(i) == 0.0_DP) then
          RAISE_ERROR("d = 0! This leads to problems computing c**2/d**2. Please specify d != 0.", ierror)
       endif
       this%c_d(i)   = this%c_sq(i)/this%d_sq(i)

       if (this%db%S(i) <= 1.0_DP) then
          RAISE_ERROR("S <= 1! This leads to problems computing (S-1)**(-1). Please specify S > 1.", ierror)
       endif
       this%VR_f(i)  = this%db%D0(i)/(this%db%S(i)-1)
       this%VA_f(i)  = this%db%S(i)*this%db%D0(i)/(this%db%S(i)-1)

       !
       ! cutoff constants.
       !

       this%cut_in_l(i)    = this%db%r1(i)
       this%cut_in_h(i)    = this%db%r2(i)
       this%cut_in_h2(i)   = this%db%r2(i)**2
       this%cut_in_fca(i)  = PI / ( this%db%r2(i) - this%db%r1(i) )
       this%cut_in_fc(i)   = - 0.5_DP * this%cut_in_fca(i)

#ifdef SCREENING

       !
       ! screening cutoffs
       !

       this%cut_out_l(i)    = this%db%or1(i)
       this%cut_out_h(i)    = this%db%or2(i)
       this%cut_out_fca(i)  = PI / ( this%db%or2(i) - this%db%or1(i) )
       this%cut_out_fc(i)   = - 0.5_DP * this%cut_out_fca(i)

       this%cut_bo_l(i)    = this%db%bor1(i)
       this%cut_bo_h(i)    = this%db%bor2(i)
       this%cut_bo_fca(i)  = PI / ( this%db%bor2(i) - this%db%bor1(i) )
       this%cut_bo_fc(i)   = - 0.5_DP * this%cut_bo_fca(i)

       this%max_cut_sq(i)   = max( &
            this%cut_in_h(i), &
            this%cut_out_h(i), &
            this%cut_bo_h(i) &
            )**2

#endif

    enddo

#ifdef SCREENING

    x = sqrt(this%C_dr_cut(1:npairs))

#endif

    do i = 1, p%nel
       do j = 1, p%nel
          ii = this%Z2db(p%el2Z(i))
          jj = this%Z2db(p%el2Z(j))
          nel = Z2pair(this, ii, jj)
#ifdef SCREENING
          call request_interaction_range( &
               nl, &
               x(nel)*sqrt(this%max_cut_sq(nel)), &
               i, j &
               )
#else
          call request_interaction_range( &
               nl, &
               this%cut_in_h(nel), &
               i, j &
               )
#endif
       enddo
    enddo

  endsubroutine BIND_TO_FUNC


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine COMPUTE_FUNC(this, p, nl, epot, f, wpot, mask, epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, ierror)
    implicit none

    type(BOP_TYPE),     intent(inout) :: this
    type(particles_t),  intent(inout) :: p
    type(neighbors_t),  intent(inout) :: nl
    real(DP),           intent(inout) :: epot
    real(DP),           intent(inout) :: f(3, p%maxnatloc)  !< forces
    real(DP),           intent(inout) :: wpot(3, 3)
    integer,  optional, intent(in)    :: mask(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_at(p%maxnatloc)
    real(DP), optional, intent(inout) :: epot_per_bond(nl%neighbors_size)
    real(DP), optional, intent(inout) :: f_per_bond(3, nl%neighbors_size)
#ifdef LAMMPS
    real(DP), optional, intent(inout) :: wpot_per_at(6, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot_per_bond(6, nl%neighbors_size)
#else
    real(DP), optional, intent(inout) :: wpot_per_at(3, 3, p%maxnatloc)
    real(DP), optional, intent(inout) :: wpot_per_bond(3, 3, nl%neighbors_size)
#endif
    integer,  optional, intent(inout) :: ierror

    ! ---

    integer  :: i, d, el(p%maxnatloc), nebmax, nebavg

    ! ---

    call timer_start(BOP_NAME_STR // "_force")

    call update(nl, p, ierror)
    PASS_ERROR(ierror)

    el = -1
    nebmax = 0
    nebavg = 0
    do i = 1, p%nat
       if (p%el2Z(p%el(i)) > 0) then
         el(i)  = this%Z2db(p%el2Z(p%el(i)))
       endif
       d = nl%last(i)-nl%seed(i)+1
       nebmax = max(nebmax, d)
       nebavg = nebavg + d
    enddo
    nebavg = (nebavg+1)/max(p%nat, 1)+1

#ifdef LAMMPS
    call BOP_KERNEL( &
         this, &
         p%maxnatloc, p%natloc, p%nat, p%r_non_cyc, el, &
         nebmax, nebavg, nl%seed, nl%last, nl%neighbors, nl%neighbors_size, &
         epot, f, wpot, mask, &
         epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
         ierror)
#else
    call BOP_KERNEL( &
         this, p%Abox, &
         p%maxnatloc, p%natloc, p%nat, p%r_non_cyc, el, &
         nebmax, nebavg, nl%seed, nl%last, nl%neighbors, nl%neighbors_size, &
         nl%dc, p%shear_dx, &
         epot, f, wpot, mask, &
         epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
         ierror)
#endif
    PASS_ERROR(ierror)

    call timer_stop(BOP_NAME_STR // "_force")

  endsubroutine COMPUTE_FUNC

