!! ======================================================================
!! MDCORE - Interatomic potential library
!! https://github.com/pastewka/mdcore
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================
!**********************************************************************
! This files contains all the subroutines needed for initialization,
! etc. so this can be used with the dipatch module.
!**********************************************************************


  !>
  !! Constructor
  !<
  subroutine INIT_FUNC(this, db)
    implicit none

    type(BOP_TYPE),              intent(inout) :: this
    type(BOP_DB_TYPE), optional, intent(in)    :: db

    ! ---

    integer :: i

    ! ---

    call prlog("- " // BOP_NAME_STR // " -")

    if (present(db)) then
       this%db = db

       call prlog("     Using database: " // trim(this%db%ref))
    endif

    do i = 1, this%db%nel
       call prlog("     el("//i//")  = " // a2s(this%db%el(:,i)))
    enddo
    call prlog("     A       = " // this%db%A(1:this%db%nA))
    call prlog("     B       = " // this%db%B(1:this%db%nB))
    call prlog("     lambda1 = " // this%db%lambda1(1:this%db%nlambda1))
    call prlog("     lambda2 = " // this%db%lambda2(1:this%db%nlambda2))
    call prlog("     eta     = " // this%db%eta(1:this%db%neta))
    call prlog("     delta   = " // this%db%delta(1:this%db%ndelta))
    call prlog("     alpha   = " // this%db%alpha(1:this%db%nalpha))
    call prlog("     beta    = " // this%db%beta(1:this%db%nbeta))
    call prlog("     c1      = " // this%db%c1(1:this%db%nc1))
    call prlog("     c2      = " // this%db%c2(1:this%db%nc2))
    call prlog("     c3      = " // this%db%c3(1:this%db%nc3))
    call prlog("     c4      = " // this%db%c4(1:this%db%nc4))
    call prlog("     c5      = " // this%db%c5(1:this%db%nc5))
    call prlog("     h       = " // this%db%h(1:this%db%nh))
    call prlog("     r1      = " // this%db%r1(1:this%db%nr1))
    call prlog("     r2      = " // this%db%r2(1:this%db%nr2))
#ifdef SCREENING
    call prlog("     or1     = " // this%db%or1(1:this%db%nor1))
    call prlog("     or2     = " // this%db%or2(1:this%db%nor2))
    call prlog("     bor1    = " // this%db%bor1(1:this%db%nbor1))
    call prlog("     bor2    = " // this%db%bor2(1:this%db%nbor2))
    call prlog("     Cmin    = " // this%db%Cmin(1:this%db%nCmin))
    call prlog("     Cmax    = " // this%db%Cmax(1:this%db%nCMax))
#endif

  endsubroutine INIT_FUNC


#include "../default_del_func.f90"
#include "../default_bind_to_func.f90"
#include "../default_compute_func.f90"
