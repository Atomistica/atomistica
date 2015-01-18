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
       el, D0, r0, S, beta, gamma, c, d, h, mu, n, m, r1, r2 &
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
    real(DP), intent(in),    optional        :: mu(:)
    real(DP), intent(in),    optional        :: n(:)
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
       if (this%ref(1) /= "*") then
          j = -1
          do i = 1, size(BOP_DB)
             if (trim(a2s(this%ref)) == BOP_DB(i)%ref(1:len_trim(a2s(this%ref)))) then
                if (j > 0) then
                   RAISE_ERROR("Reference string '" // a2s(this%ref) // "' not unique. Matching entries: '" // BOP_DB(j)%ref // "' and '" // BOP_DB(i)%ref // "'.", ierror)
                endif
                j = i
             endif
          enddo

          if (j > 0) then
             this%db = BOP_DB(j)
          else
             RAISE_ERROR("Could not find parameter set for reference '" // a2s(this%ref) // "' in database.", ierror)
          endif

       endif

       call prlog("     Using database: " // trim(this%db%ref))

       ASSIGN_STRING_ARRAY_PROPERTY(el,    this%db%el,    this%db%nel, i)
       ASSIGN_ARRAY_PROPERTY(D0,    this%db%D0,    this%db%nD0)
       ASSIGN_ARRAY_PROPERTY(r0,    this%db%r0,    this%db%nr0)
       ASSIGN_ARRAY_PROPERTY(S,     this%db%S,     this%db%nS)
       ASSIGN_ARRAY_PROPERTY(beta,  this%db%beta,  this%db%nbeta)
       ASSIGN_ARRAY_PROPERTY(gamma, this%db%gamma, this%db%ngamma)
       ASSIGN_ARRAY_PROPERTY(c,     this%db%c,     this%db%nc)
       ASSIGN_ARRAY_PROPERTY(d,     this%db%d,     this%db%nd)
       ASSIGN_ARRAY_PROPERTY(h,     this%db%h,     this%db%nh)
       ASSIGN_ARRAY_PROPERTY(mu,    this%db%mu,    this%db%nmu)
       ASSIGN_ARRAY_PROPERTY(n,     this%db%n,     this%db%nn)
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

    do i = 1, this%db%nel
       call prlog("     el("//i//")  = " // a2s(this%db%el(:,i)))
    enddo
    call prlog("     D0     = " // this%db%D0(1:this%db%nD0))
    call prlog("     r0     = " // this%db%r0(1:this%db%nr0))
    call prlog("     S      = " // this%db%S(1:this%db%nS))
    call prlog("     beta   = " // this%db%beta(1:this%db%nbeta))
    call prlog("     gamma  = " // this%db%gamma(1:this%db%ngamma))
    call prlog("     c      = " // this%db%c(1:this%db%nc))
    call prlog("     d      = " // this%db%d(1:this%db%nd))
    call prlog("     h      = " // this%db%h(1:this%db%nh))
    call prlog("     mu     = " // this%db%mu(1:this%db%nmu))
    call prlog("     n      = " // this%db%n(1:this%db%nn))
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


#include "../default_del_func.f90"


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

#ifdef SCREENING
    real(DP)         :: x(this%db%nel*(this%db%nel+1)/2)
#endif

    ! ---

    INIT_ERROR(ierror)

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

    npairs = this%db%nel*(this%db%nel+1)/2

    if (npairs /= this%db%nD0 .or. &
        npairs /= this%db%nr0 .or. &
        npairs /= this%db%nS  .or. &
        npairs /= this%db%nbeta .or. &
        npairs /= this%db%ngamma .or. &
        npairs /= this%db%nc .or. &
        npairs /= this%db%nd .or. &
        npairs /= this%db%nh .or. &
        npairs /= this%db%nn .or. &
        npairs /= this%db%nm .or. &
        npairs /= this%db%nmu .or. &
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
       write (*, '(A,I2)')  "nn      = ", this%db%nn
       write (*, '(A,I2)')  "nm      = ", this%db%nm
       write (*, '(A,I2)')  "nmu     = ", this%db%nmu
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

    this%Z2db  = 0

    do i = 1, this%db%nel
       Z = atomic_number(a2s(this%db%el(:, i)))
       if (Z > 0) then
          this%Z2db(Z)  = i
       else
          RAISE_ERROR("Unknown element '" // trim(a2s(this%db%el(:, i))) // "'.", ierror)
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

       call init(this%cut_in(i), this%db%r1(i), this%db%r2(i))

       this%cut_in_l(i)    = this%db%r1(i)
       this%cut_in_h(i)    = this%db%r2(i)
       this%cut_in_h2(i)   = this%db%r2(i)**2

#ifdef SCREENING

       !
       ! screening cutoffs
       !

       call init(this%cut_out(i), this%db%or1(i), this%db%or2(i))

       this%cut_out_l(i)    = this%db%or1(i)
       this%cut_out_h(i)    = this%db%or2(i)

       call init(this%cut_bo(i), this%db%bor1(i), this%db%bor2(i))

       this%cut_bo_l(i)    = this%db%bor1(i)
       this%cut_bo_h(i)    = this%db%bor2(i)

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


#include "../default_compute_func.f90"
