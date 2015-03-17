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
!    shared
! @endmeta

!**********************************************************************
! This is the kernel for the bond-order potentials of the
! Tersoff-Brenner type. Currently works with the Erhart-Albe,
! Tersoff and Brenner potential.
!**********************************************************************


#ifndef SCREENING
#define cutfcnbo  cutfcnar
#define cutdrvbo  cutdrvar
#define cutfcnnc  cutfcnar
#define cutdrvnc  cutdrvar

#define cut_ar_h  cut_in_h
#define cut_bo_h  cut_in_h
#endif

#ifndef LAMMPS
#define DCELL_INDEX(ni)  VEC(dc, ni, 3) + (2*maxdc(3)+1) * ( VEC(dc, ni, 2) + (2*maxdc(2)+1) * VEC(dc, ni, 1) )
#endif

#ifdef _OPENMP
#define NEB_TOO_SMALL(what, i, ierror)  RAISE_DELAYED_ERROR("Internal neighbor list exhausted on OpenMP thread " // omp_get_thread_num() // ", *" // what // "* too small: " // "nebtot = " // nebtot // "/" // nebsize // ", nebmax = " // nebmax // ", nebavg = " // nebavg // ", neb_last(i)-neb_seed(i)+1 = " // (neb_last(i)-neb_seed(i)+1) // ", i = " // i // "/" // natloc // "(" // nat // ")", ierror) ; nebtot = int(1 + real(nebsize, DP)*omp_get_thread_num()/omp_get_num_threads()) ; neb_last(i) = neb_seed(i)
#define SNEB_TOO_SMALL(what, i, ierror)  RAISE_DELAYED_ERROR("Internal screening neighbor list exhausted on OpenMP thread " // omp_get_thread_num() // ", *" // what // "* too small: " // "snebtot = " // snebtot // "/" // snebsize // ", nebmax = " // nebmax // ", nebavg = " // nebavg // ", this%sneb_last(i)-this%sneb_seed(i)+1 = " // (this%sneb_last(i)-this%sneb_seed(i)+1) // ", i = " // i // "/" // natloc // "(" // nat // ")", ierror) ; snebtot = int(1 + real(snebsize, DP)*omp_get_thread_num()/omp_get_num_threads()) ; this%sneb_last(i) = this%sneb_seed(i)
#else
#define NEB_TOO_SMALL(what, i, ierror)  RAISE_ERROR("Internal neighbor list exhausted, *" // what // "* too small: " // "nebtot = " // nebtot // "/" // nebsize // ", nebmax = " // nebmax // ", nebavg = " // nebavg // ", neb_last(i)-neb_seed(i)+1 = " // (neb_last(i)-neb_seed(i)+1) // ", i = " // i // "/" // natloc // " (" // nat // ")", ierror) ; nebtot = 1 ; neb_last(i) = neb_seed(i)
#define SNEB_TOO_SMALL(what, i, ierror)  RAISE_ERROR("Internal screening neighbor list exhausted, *" // what // "* too small: " // "snebtot = " // snebtot // "/" // snebsize // ", nebmax = " // nebmax // ", nebavg = " // nebavg // ", this%sneb_last(i)-this%sneb_seed(i)+1 = " // (this%sneb_last(i)-this%sneb_seed(i)+1) // ", i = " // i // "/" // natloc // " (" // nat // ")", ierror) ; snebtot = 1 ; this%sneb_last(i) = this%sneb_seed(i)
#endif

#ifdef LAMMPS
  recursive subroutine BOP_KERNEL( &
       this, &
       maxnat, natloc, nat, r, &
       el, &
       nebmax, nebavg, aptr, a2ptr, bptr, ptrmax, &
       epot, f_inout, wpot_inout, mask, &
       epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
       ierror)
#else
  recursive subroutine BOP_KERNEL( &
       this, cell, &
       maxnat, natloc, nat, r, &
       el, &
       nebmax, nebavg, aptr, a2ptr, bptr, ptrmax, dc, shear_dx, &
       epot, f_inout, wpot_inout, mask, &
       epot_per_at, epot_per_bond, f_per_bond, wpot_per_at, wpot_per_bond, &
       ierror)
#endif

    ! 
    ! copyright: keith beardmore 28/11/93.
    ! - algorithm from : phys. rev. b 42, 9458-9471(1990).
    ! - plus corrections : phys. rev. b 46, 1948(1990).
    ! - modified to use linked list and pointers
    ! - to reduce size of neb-list. 25/1/94.
    ! - pre-calculates pairwise terms. 29/1/94.
    ! 
    ! copyright: lars pastewka 2006-2009
    ! - made Fortran 90 compliant
    ! - screening functions
    !     M. I. Baskes et al., Modelling Simul. Mater. Sci. Eng. 2, 505 (1994)
    !     L. Pastewka et al., Phys. Rev. B 78, 161402(R) (2008)
    ! - optimizations. 07/2007
    ! - virial. 09/2008
    ! - OO compliant. 02/2009
    ! 
    ! 
    ! algorithm assumes atoms are el=1 (c), el=3 (h)
    ! 
    ! for REBO :
    ! 
    !       o o o       p p p
    !        \|/         \|/
    !         m m m   n n n              the energy of bond i-j is
    !          \|/     \|/               dependent on all atoms that are 
    !   o   m   k       l   n   p        first, second (or third neighbours
    !    \   \   \     /   /   /         if screening is enabled) of 
    ! o---m---k---i===j---l---n---p      i and j. the resulting forces   
    !    /   /   /     \   \   \         act upon all these atoms.
    !   o   m   k       l   n   p        
    !          /|\     /|\               in the code, the atoms and
    !         m m m   n n n              bonds are identified as shown
    !        /|\         /|\             on the left.
    !       o o o       p p p
    ! 

    use tls

#ifdef _OPENMP
    use omp_lib
#endif

    implicit none

    integer, parameter       :: typemax = 3

    ! ---

    type(BOP_TYPE),      intent(inout) :: this
#ifndef LAMMPS
    real(DP),            intent(in)    :: cell(3, 3)
#endif

    integer,             intent(in)    :: maxnat
    integer,             intent(in)    :: natloc
    integer,             intent(in)    :: nat
    real(DP),            intent(inout) :: r(3, maxnat)

    integer,             intent(in)    :: ptrmax

    real(DP),            intent(inout) :: f_inout(3, maxnat)
    real(DP),            intent(inout) :: epot
    real(DP),            intent(inout) :: wpot_inout(3, 3)
    real(DP)                           :: wpot(3, 3)

    integer,             intent(in)    :: el(maxnat)

    integer,   optional, intent(in)    :: mask(maxnat)

    real(DP),  optional, intent(inout) :: epot_per_at(nat)
    real(DP),  optional, intent(inout) :: epot_per_bond(ptrmax)

    integer,             intent(in)    :: nebmax, nebavg
    integer(NEIGHPTR_T), intent(in)    :: aptr(maxnat+1)
    integer(NEIGHPTR_T), intent(in)    :: a2ptr(maxnat+1)
    integer,             intent(in)    :: bptr(ptrmax)

#ifdef LAMMPS
    real(DP),  optional, intent(inout) :: wpot_per_at(6, maxnat)
    real(DP),  optional, intent(inout) :: wpot_per_bond(6, ptrmax)
#else
    integer,             intent(in)    :: dc(3, ptrmax)
    real(DP),            intent(in)    :: shear_dx(3)

    real(DP),  optional, intent(inout) :: wpot_per_at(3, 3, maxnat)
    real(DP),  optional, intent(inout) :: wpot_per_bond(3, 3, ptrmax)
#endif

    real(DP), optional, intent(inout)  :: f_per_bond(3, ptrmax)

    integer,  optional, intent(inout)  :: ierror

    ! ---

    ! "short" neighbor list (all neighbors which are not screened)

    integer(NEIGHPTR_T) :: jbeg,jend,jn

    real(DP)  :: rij(3)
    real(DP)  :: rlij, rlijr, rlik
    real(DP)  :: rnij(3), rnik(3)
    real(DP)  :: df(3)
    real(DP)  :: fcarij,dfcarijr,fcik,dfcikr
    real(DP)  :: VAij,dVAij_drij,VRij,dVRij_drij
    real(DP)  :: zij
    real(DP)  :: wij(3, 3), wijb(3, 3)
    real(DP)  :: dbidi(3), dbidj(3), dbidk(3, nebmax)
    real(DP)  :: disjk,dkc(3)
    real(DP)  :: costh
    real(DP)  :: g_costh,dg_dcosth
    real(DP)  :: h_Dr,dh_dDr
#ifdef SEPARATE_H_ARGUMENTS
    real(DP)  :: dh_dr2
#endif
    real(DP)  :: dzfac,dzdrij,dzdrik
    real(DP)  :: dgdi(3), dgdj(3), dgdk(3)
    real(DP)  :: dcsdij,dcsdik,dcsdjk
    real(DP)  :: dcsdi(3), dcsdj(3), dcsdk(3)
    real(DP)  :: bij
    real(DP)  :: dbij_dzij
#ifdef BO_WITH_D
    real(DP)  :: dbij_dDij, Dij, dDij_drij
#endif
    real(DP)  :: dffac

    real(DP)  :: rik(3)

    real(DP)  :: fi(3), fj(3)

#if defined(SCREENING)
    real(DP)  :: xik
#endif

#ifdef SCREENING
    real(DP)  :: rljk, dot_ij_ik, dot_ij_jk
    real(DP)  :: rjk(3)
    real(DP)  :: C, dCdrij, dCdrik, dCdrjk, Cmax_C, C_Cmin, fac
    real(DP)  :: xjk, xik_p_xjk, xik_m_xjk
#ifdef SIN_S
    real(DP)  :: csij
#endif
    real(DP)  :: sij, dsijdrij, dsijdrik, dsijdrjk
    real(DP)  :: fcboij, dfcboijr
    real(DP)  :: zfaci(nebmax)
#endif
    real(DP)  :: fcinij, dfcinijr

    integer  :: i,j,k
#ifndef LAMMPS
    integer  :: jdc,kdc
#endif
    integer  :: ij,ik

    integer  :: ikc
    integer  :: maskfac
    integer  :: eli,elj,elk
    integer  :: el2ij,ikpot
    ! seed for the "short" neighbor list
    integer  :: nebtot, neb_seed(nat), neb_last(nat), istart, ifinsh
    ! seed for the "screened" neighbor list
    integer  :: nebofi(nebmax)
#ifndef LAMMPS
    integer  :: dcofi(nebmax)
#endif

    integer  :: numnbi

    integer  :: nebsize

#ifdef SCREENING
    integer  :: seedi(nebmax)
    integer  :: lasti(nebmax)

    integer  :: snebsize

    integer  :: i1,i2
    integer(NEIGHPTR_T) :: kn

    integer  :: nijc
    integer  :: snebtot, ineb

    logical  :: screened
    logical  :: need_derivative
#endif

#ifndef LAMMPS
    integer  :: maxdc(3)
#endif

#ifdef _OPENMP
    integer  :: ierror_loc
#else
#define ierror_loc ierror
#endif

    ! ---

#ifdef DEBUG_OUTPUT

    if (.not. allocated(debug_S)) then
       allocate(debug_S(ptrmax))
       allocate(debug_fC(ptrmax))
       allocate(debug_rl(ptrmax))
    endif

    debug_S   = 0.0_DP
    debug_fC  = 0.0_DP
    debug_rl  = 0.0_DP

#endif

    this%it = this%it + 1

    ! This size should be sufficient, buffers should not overflow.
#ifdef _OPENMP
    nebsize = max(omp_get_max_threads()**2*nebmax**2, &
                  min((nat+1)*nebmax, ptrmax))+omp_get_max_threads()
#else
    nebsize = min(nat*nebavg, ptrmax)+1
#endif
#ifdef SCREENING
    ! This size can overflow. However, most bond are either screened or not
    ! screened such that number is expected to be low.
    snebsize = nebsize
#endif

#ifdef SCREENING
    if (this%neighbor_list_allocated .and. &
        (size(this%neb) < nebsize .or. size(this%sneb) < snebsize)) then
#else
    if (this%neighbor_list_allocated .and. size(this%neb) < nebsize) then
#endif

       this%neighbor_list_allocated = .false.

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

    endif

    !-------prepare brenner material
    if (.not. this%neighbor_list_allocated) then

       call prlog("- " // BOP_NAME_STR // " -")
#ifdef SCREENING
       call prlog("The " // BOP_NAME_STR // " potential has been compiled with screening functions.")
#endif
       call prlog("(Re)allocating internal neighbor list buffers.")
       call prlog("nebavg   = " // nebavg)
       call prlog("nebmax   = " // nebmax)
       call prlog("nebsize  = " // nebsize)
#ifdef SCREENING
       call prlog("snebsize = " // snebsize)
#endif

       call log_memory_start(BOP_NAME_STR)

       allocate(this%neb(nebsize))
       allocate(this%nbb(nebsize))
#ifndef LAMMPS
       allocate(this%dcell(nebsize))
#endif
       allocate(this%bndtyp(nebsize))
       allocate(this%bndlen(nebsize))
       allocate(this%bndnm(3, nebsize))
       allocate(this%cutfcnar(nebsize))
       allocate(this%cutdrvar(nebsize))

       call log_memory_estimate(this%neb)
       call log_memory_estimate(this%nbb)
#ifndef LAMMPS
       call log_memory_estimate(this%dcell)
#endif
       call log_memory_estimate(this%bndtyp)
       call log_memory_estimate(this%bndlen)
       call log_memory_estimate(this%bndnm)
       call log_memory_estimate(this%cutfcnar)
       call log_memory_estimate(this%cutdrvar)

#ifdef SCREENING
       allocate(this%cutfcnbo(nebsize))
       allocate(this%cutdrvbo(nebsize))
       allocate(this%sneb_seed(nebsize))
       allocate(this%sneb_last(nebsize))
       allocate(this%sneb(snebsize))
       allocate(this%sbnd(snebsize))
       allocate(this%sfacbo(snebsize))
       allocate(this%cutdrarik(snebsize))
       allocate(this%cutdrarjk(snebsize))
       allocate(this%cutdrboik(snebsize))
       allocate(this%cutdrbojk(snebsize))

       call log_memory_estimate(this%cutfcnbo)
       call log_memory_estimate(this%cutdrvbo)
       call log_memory_estimate(this%sneb_seed)
       call log_memory_estimate(this%sneb_last)
       call log_memory_estimate(this%sneb)
       call log_memory_estimate(this%sbnd)
       call log_memory_estimate(this%sfacbo)
       call log_memory_estimate(this%cutdrarik)
       call log_memory_estimate(this%cutdrarjk)
       call log_memory_estimate(this%cutdrboik)
       call log_memory_estimate(this%cutdrbojk)
#endif

       call log_memory_stop(BOP_NAME_STR)

       call prlog

       this%neighbor_list_allocated = .true.
       
    endif

    ! 
    ! set all p.e. and forces to zero.
    !

    wpot  = 0.0_DP

#ifndef LAMMPS
    maxdc = 0
    do i = 1, nat
       maxdc(1)  = max(maxdc(1), maxval(VEC(dc, aptr(i):a2ptr(i), 1)))
       maxdc(2)  = max(maxdc(2), maxval(VEC(dc, aptr(i):a2ptr(i), 2)))
       maxdc(3)  = max(maxdc(3), maxval(VEC(dc, aptr(i):a2ptr(i), 3)))
    enddo
#endif

!    write (*, *)  "maxdc  = ", maxdc

#ifdef SCREENING
    this%sfacbo  = 0.0_DP
#endif

    ! 
    ! calculate the main pairwise terms and store them.
    !

#ifdef _OPENMP
    ierror_loc = ERROR_NONE
#else
    INIT_ERROR(ierror_loc)
#endif

    !$omp  parallel default(none) &
    !$omp& shared(aptr, a2ptr, bptr, f_inout, el) &
#ifndef LAMMPS
    !$omp& shared(cell, dc, shear_dx) &
#endif
    !$omp& firstprivate(nat, natloc, nebmax, nebavg, nebsize) &
    !$omp& shared(neb_last, neb_seed) &
    !$omp& shared(mask, epot_per_at, epot_per_bond) &
    !$omp& shared(r, this, f_per_bond, wpot_per_at, wpot_per_bond) &
#ifdef SCREENING
    !$omp& firstprivate(snebsize) &
#endif
#ifndef LAMMPS
    !$omp& firstprivate(maxdc) &
#endif
    !$omp& private(jbeg,jend,jn) &
    !$omp& private(rij, rik) &
    !$omp& private(rlij, rlijr, rlik) &
    !$omp& private(rnij, rnik) &
    !$omp& private(maskfac,df) &
#ifdef BO_WITH_D
    !$omp& private(Dij, dDij_drij, dbij_dDij) &
#endif
    !$omp& private(fcarij,dfcarijr,fcik,dfcikr) &
    !$omp& private(VAij,dVAij_drij,VRij,dVRij_drij) &
    !$omp& private(zij) &
    !$omp& private(wij, wijb) &
    !$omp& private(dbidi, dbidj, dbidk) &
    !$omp& private(disjk,dkc) &
    !$omp& private(costh) &
    !$omp& private(g_costh,dg_dcosth) &
    !$omp& private(h_Dr,dh_dDr) &
#ifdef SEPARATE_H_ARGUMENTS
    !$omp& private(dh_dr2) &
#endif
    !$omp& private(dzfac,dzdrij,dzdrik) &
    !$omp& private(dgdi, dgdj, dgdk) &
    !$omp& private(dcsdij,dcsdik,dcsdjk) &
    !$omp& private(dcsdi, dcsdj, dcsdk) &
    !$omp& private(bij,dbij_dzij) &
    !$omp& private(dffac) &
    !$omp& private(fi, fj) &
#if defined(SCREENING)
    !$omp& private(xik) &
#endif
#ifdef SCREENING
    !$omp& private(rljk, dot_ij_ik, dot_ij_jk) &
    !$omp& private(rjk) &
    !$omp& private(fcboij, dfcboijr) &
    !$omp& private(zfaci) &
#endif
    !$omp& private(fcinij, dfcinijr) &
    !$omp& private(i,j,k) &
#ifndef LAMMPS
    !$omp& private(jdc,kdc) &
#endif
    !$omp& private(ij,ik) &
    !$omp& private(ikc) &
    !$omp& private(eli,elj,elk) &
    !$omp& private(el2ij,ikpot) &
    !$omp& private(istart, ifinsh) &
    !$omp& private(nebofi) &
#ifndef LAMMPS
    !$omp& private(dcofi) &
#endif
    !$omp& private(numnbi)&
    !$omp& private(nebtot) &
#ifdef SCREENING
    !$omp& private(i1,i2,kn) &
    !$omp& private(seedi) &
    !$omp& private(lasti) &
    !$omp& private(nijc,ineb) &
    !$omp& private(snebtot) &
    !$omp& private(sij, dsijdrij, dsijdrik, dsijdrjk) &
    !$omp& private(C, dCdrij, dCdrik, dCdrjk, Cmax_C, C_Cmin, fac) &
    !$omp& private(xjk, xik_p_xjk, xik_m_xjk) &
#ifdef SIN_S
    !$omp& private(csij) &
#endif
    !$omp& private(screened, need_derivative) &
#endif
    !$omp& reduction(+:ierror_loc) &
    !$omp& reduction(+:wpot) reduction(+:epot)

    call tls_init(nat, sca=1, vec=1)
#define pe tls_sca1
#define f  tls_vec1

#define nebmax_sq zij
    nebmax_sq = nebmax*nebmax

    ! Convert to real to avoid overflow
#ifdef _OPENMP

    ! When using OpenMP parallelization, every thread gets an equal share of the
    ! internal neighbor list buffers.
    nebtot   = int(1 + real(nebsize, DP)*omp_get_thread_num()/omp_get_num_threads())
    nebsize  = int(real(nebsize, DP)*(omp_get_thread_num()+1)/omp_get_num_threads())
#ifdef SCREENING
    snebtot  = int(1 + real(snebsize, DP)*omp_get_thread_num()/omp_get_num_threads())
    snebsize = int(real(snebsize, DP)*(omp_get_thread_num()+1)/omp_get_num_threads())
#endif ! SCREENING

#else

    nebtot   = 1
#ifdef SCREENING
    snebtot  = 1
#endif ! SCREENING

#endif ! _OPENMP

    !$omp do
    i_loop1: do i = 1, natloc
       eli = el(i)

       neb_seed(i)  = nebtot
       neb_last(i)  = nebtot-1

       i_known_el1: if (eli > 0) then

          jbeg = aptr(i)
          jend = a2ptr(i)

          jn_loop1: do jn = jbeg, jend

             !
             ! Loop over all pairs
             !

#ifdef LAMMPS
             j    = bptr(jn)+1
#else
             j    = bptr(jn)
             jdc  = DCELL_INDEX(jn)
#endif
             elj  = el(j)

             j_known_el1: if (elj > 0) then

#ifdef LAMMPS
                rij  = VEC3(r, j) - VEC3(r, i)
#else
                rij  = VEC3(r, j) - VEC3(r, i) - matmul(cell, VEC3(dc, jn)) - &
                     shear_dx*VEC(dc, jn, 3)
#endif

                rlij = dot_product(rij, rij)

                el2ij = Z2pair(this, eli, elj)

                !
                ! ...NO PARTIAL SCREENING...
                !  There are different regions that need to be handled 
                !  differently
                !  -------------- r1 ------ r2 ----------- r3 -------- r4
                !   no screening     trans.     screening      cutoff
                !       (a)           (b)          (c)           (d)
                !                 ^         ^                          ^
                !              cut_in_l  cut_in_h                   max_cut
                !
                ! ...PARTIAL SCREENING...
                !  No difference in region, but "inner cutoff" (fCin) is always
                !  on and unscreened and "outer cutoff" (fCbo, fCar) is
                !  screened. Potential needs to ensure sum equals to one.
                !
                
#ifdef DEBUG_OUTPUT

                debug_rl(jn) = sqrt(rlij)

#endif

!                write (*, *)  "A: ", sqrt(rlij), el2ij

#ifndef PARTIAL_SCREENING                   
                cutoff_region: if (rlij < this%cut_in_l(el2ij)**2) then

!                   write (*, *)  "B: ", sqrt(rlij), el2ij
                
                   !
                   ! In region (a) -> atoms are allowed to interact
                   !

                   this%cutfcnar(nebtot)  = 1.0_DP
                   this%cutdrvar(nebtot)  = 0.0_DP

#ifdef SCREENING
                   this%cutfcnbo(nebtot)  = 1.0_DP
                   this%cutdrvbo(nebtot)  = 0.0_DP
#endif

#ifdef DEBUG_OUTPUT
                   debug_S(jn)            = 1.0_DP
                   debug_fC(jn)           = 1.0_DP
#endif

                   this%neb(nebtot)       = j
                   this%nbb(nebtot)       = jn
#ifndef LAMMPS
                   this%dcell(nebtot)     = jdc
#endif

#ifdef SCREENING
                   this%sneb_seed(nebtot) = snebtot
                   this%sneb_last(nebtot) = snebtot-1
#endif

                   ! 
                   ! bond-length and direction cosines.
                   !
                   
                   rlij                     = sqrt( rlij )
                   this%bndlen(nebtot)      = rlij
                   this%bndnm(1:3, nebtot)  = rij / rlij
                   this%bndtyp(nebtot)      = el2ij

                   neb_last(i)            = nebtot
                   nebtot                 = nebtot + 1

                   if (neb_last(i)-neb_seed(i)+1 > nebmax) then
                      NEB_TOO_SMALL("nebmax", i, ierror_loc)
                   endif

                   if (nebtot > nebsize) then
                      NEB_TOO_SMALL("nebsize", i, ierror_loc)
                   endif

#endif
      
#ifdef SCREENING

#ifdef PARTIAL_SCREENING
                cutoff_region: if (rlij < this%max_cut_sq(el2ij) .and. &
                     this%cut_out_l(el2ij) < this%cut_out_h(el2ij)) then
#else
                else if (rlij < this%max_cut_sq(el2ij) .and. &
                     this%cut_out_l(el2ij) < this%cut_out_h(el2ij)) then
#endif

                   !
                   ! Compute screening function
                   !

                   screened                = .false.
                   need_derivative         = .false.
#ifdef SIN_S
                   sij       = 1.0_DP
#else
                   sij       = 0.0_DP
#endif

                   this%sneb_seed(nebtot)  = snebtot
                   this%sneb_last(nebtot)  = snebtot-1

                   !
                   ! within cutoff: compute the screening function for the
                   ! bond i-j
                   !
                      
                   ineb      = snebtot

                   dsijdrij  = 0.0_DP

                   kn = jbeg
                   do while (.not. &
                        (screened .or. &
                        sij < this%screening_threshold) .and. &
                        kn <= jend)

#ifdef LAMMPS
                      k = bptr(kn)+1
                      rik  = VEC3(r, k) - VEC3(r, i)
#else
                      k = bptr(kn)
                      rik  = VEC3(r, k) - VEC3(r, i) - &
                           matmul(cell, VEC3(dc, kn)) - shear_dx*VEC(dc, kn, 3)
#endif

                      if (dot_product(rik, rik) < this%C_dr_cut(el2ij)*rlij) &
                           then

#ifdef LAMMPS
                         k_neq_j: if (k /= j) then
#else
                         k_neq_j: if (k /= j .or. &
                              any(VEC3(dc, kn) /= VEC3(dc, jn))) then
#endif

                            dot_ij_ik = dot_product(rij, rik)

                            rlik = dot_product(rik, rik)

                            rjk  = -rij + rik

                            dot_ij_jk = dot_product(rij, rjk)

                            rljk = dot_product(rjk, rjk)

                            if (dot_ij_ik > this%dot_threshold .and. &
                                 dot_ij_jk < -this%dot_threshold) then

                               xik  = rlik/rlij
                               xjk  = rljk/rlij

                               xik_m_xjk  = xik-xjk
                               xik_p_xjk  = xik+xjk

                               fac        = 1.0_DP/(1-xik_m_xjk**2)

                               C  = (2*(xik_p_xjk)-(xik_m_xjk)**2-1)*fac

                               if (C <= this%Cmin(el2ij)) then
                                  screened = .true.
                               else if (C < this%Cmax(el2ij)) then
                                  need_derivative = .true.

                                  Cmax_C    = this%Cmax(el2ij)-C
                                  C_Cmin    = C-this%Cmin(el2ij)

#ifdef SIN_S
                                  csij      = (1 - &
                                       cos(PI*(C-Cmin)/this%dC(el2ij)))/2
                                  sij       = sij * csij
#else
                                  sij       = sij - (Cmax_C/C_Cmin)**2
#endif

                                  dCdrik    = 4*xik*fac*(1+(C-1)*xik_m_xjk)
                                  dCdrjk    = 4*xjk*fac*(1-(C-1)*xik_m_xjk)

                                  dCdrij    = -(dCdrik+dCdrjk)

#ifdef SIN_S
                                  fac       = PI/(2*this%dC(el2ij)) * &
                                       sin(PI*(C-Cmin)/this%dC(el2ij)) / csij
#else
                                  fac       = 2*Cmax_C*this%dC(el2ij)/&
                                       (C_Cmin**3)
#endif

                                  !
                                  !  the following estimates lack a factor of 
                                  !  sij
                                  !

                                  dsijdrij  = dsijdrij + fac*dCdrij
                                  dsijdrik  = fac*dCdrik
                                  dsijdrjk  = fac*dCdrjk
                                        
                                  this%sneb(snebtot)       = k
                                  this%sbnd(snebtot)       = kn

                                  this%cutdrarik(snebtot)  = dsijdrik/rlik
                                  this%cutdrarjk(snebtot)  = dsijdrjk/rljk

                                  this%sneb_last(nebtot)   = snebtot
                                  snebtot                  = snebtot + 1

                               endif

                            endif

                         endif k_neq_j

                      endif

                      kn = kn+1

                   enddo

#ifndef PARTIAL_SCREENING
                   is_fully_screened: if ((screened .or. &
                        sij < this%screening_threshold) .and. &
                        rlij > this%cut_in_h2(el2ij)) then
                         
                      !
                      ! reset our screening neighbor because the bond is
                      ! screened anyway
                      !

                      snebtot                 = ineb
                      this%sneb_last(nebtot)  = ineb - 1

                   else
#endif

                      !
                      ! not screened by another atom: add to local neighbor
                      ! list
                      !

                      this%neb(nebtot)    = j
                      this%nbb(nebtot)    = jn
#ifndef LAMMPS
                      this%dcell(nebtot)  = DCELL_INDEX(jn)
#endif

                      ! 
                      ! bond-length and direction cosines.
                      !

                      rlij                     = sqrt( rlij )
                      this%bndlen(nebtot)      = rlij
                      this%bndnm(1:3, nebtot)  = rij / rlij
                      this%bndtyp(nebtot)      = el2ij

                      is_partially_screened: if ( screened ) then

                         call fCin(this, el2ij, rlij, fcinij, dfcinijr)

                         this%cutfcnar(nebtot)  = fcinij
                         this%cutdrvar(nebtot)  = dfcinijr

                         this%cutfcnbo(nebtot) = fcinij
                         this%cutdrvbo(nebtot) = dfcinijr

#ifdef DEBUG_OUTPUT
                         debug_S(jn)   = 0.0_DP
                         debug_fC(jn)  = this%cutfcnar(nebtot)
#endif

                         snebtot                 = ineb
                         this%sneb_last(nebtot)  = ineb - 1

                      else if ( need_derivative ) then

#ifndef SIN_S
                         sij = exp( sij )
#endif

                         call fCin(this, el2ij, rlij, fcinij, dfcinijr)
                         call fCar(this, el2ij, rlij, fcarij, dfcarijr)
                         call fCbo(this, el2ij, rlij, fcboij, dfcboijr)

                         !
                         ! do also compute the derivatives with respect to the
                         ! neighbors
                         !

                         this%cutfcnar(nebtot)  = (1.0_DP-fcinij)*sij*fcarij + &
                              fcinij
                         this%cutdrvar(nebtot)  = (1.0_DP-fcinij)*sij* &
                              (dfcarijr + fcarij*dsijdrij/rlij) - &
                              dfcinijr*sij*fcarij + dfcinijr

                         this%cutfcnbo(nebtot) = (1.0_DP-fcinij)*sij*fcboij + &
                              fcinij
                         this%cutdrvbo(nebtot) = (1.0_DP-fcinij)*sij* &
                              (dfcboijr + fcboij*dsijdrij/rlij) - &
                              dfcinijr*sij*fcboij + dfcinijr

#ifdef DEBUG_OUTPUT
                         debug_S(jn)   = sij
                         debug_fC(jn)  = this%cutfcnar(nebtot)
#endif

                         !
                         ! multiply the sij and fcarij into the derivatives
                         !

                         this%cutdrboik(ineb:snebtot-1) = &
                              this%cutdrarik(ineb:snebtot-1)*sij*fcboij * &
                              (1.0_DP-fcinij)
                         this%cutdrbojk(ineb:snebtot-1) = &
                              this%cutdrarjk(ineb:snebtot-1)*sij*fcboij * &
                              (1.0_DP-fcinij)

                         this%cutdrarik(ineb:snebtot-1)  = &
                              this%cutdrarik(ineb:snebtot-1)*sij*fcarij * &
                              (1.0_DP-fcinij)
                         this%cutdrarjk(ineb:snebtot-1)  = &
                              this%cutdrarjk(ineb:snebtot-1)*sij*fcarij * &
                              (1.0_DP-fcinij)

                      else

                         !
                         ! we don't need the derivative of the screening 
                         ! function with respect to the neighbors because the
                         ! screening function is a constant (=1, locally).
                         !

                         call fCar(this, el2ij, rlij, fcarij, dfcarijr)
                         call fCbo(this, el2ij, rlij, fcboij, dfcboijr)

#ifndef PARTIAL_SCREENING
                         if (rlij < this%cut_in_h(el2ij)) then
#endif

                            call fCin(this, el2ij, rlij, fcinij, dfcinijr)

                            this%cutfcnar(nebtot)  = (1.0_DP-fcinij)*fcarij + &
                                 fcinij
                            this%cutdrvar(nebtot)  = (1.0_DP-fcinij)*dfcarijr -&
                                 dfcinijr*fcarij + dfcinijr

                            this%cutfcnbo(nebtot) = (1.0_DP-fcinij)*fcboij + &
                                 fcinij
                            this%cutdrvbo(nebtot) = (1.0_DP-fcinij)*dfcboijr - &
                                 dfcinijr*fcboij + dfcinijr
      
#ifndef PARTIAL_SCREENING                      
                         else
                            
                            this%cutfcnar(nebtot) = fcarij
                            this%cutdrvar(nebtot) = dfcarijr

                            this%cutfcnbo(nebtot) = fcboij
                            this%cutdrvbo(nebtot) = dfcboijr

                         endif
#endif

#ifdef DEBUG_OUTPUT
                         debug_S(jn)   = 1.0_DP
                         debug_fC(jn)  = this%cutfcnar(nebtot)
#endif

                      endif is_partially_screened

                      neb_last(i)  = nebtot
                      nebtot       = nebtot + 1

                      if (neb_last(i)-neb_seed(i)+1 > nebmax) then
                         NEB_TOO_SMALL("nebmax", i, ierror_loc)
                      endif

                      if (nebtot > nebsize) then
                         NEB_TOO_SMALL("nebsize", i, ierror_loc)
                      endif
                      
                      if (snebtot > snebsize) then
                         SNEB_TOO_SMALL("snebsize", i, ierror_loc)
                      endif

#ifndef PARTIAL_SCREENING
                   endif is_fully_screened
#endif

#endif ! ifdef SCREENING

                else if (rlij < this%cut_in_h2(el2ij)) then

                   !
                   ! Either we don't have screening compiled, or this bond
                   ! shouldn't be screened.
                   !

!                   write (*, *)  "Unscreened bond"

                   ! 
                   ! bond-length and direction cosines.
                   !

                   rlij                     = sqrt( rlij )
                   this%bndlen(nebtot)      = rlij
                   this%bndnm(1:3, nebtot)  = rij / rlij
                   this%bndtyp(nebtot)      = el2ij

                   !
                   ! Cut-off function
                   !

                   call fCin(this, el2ij, rlij, fcinij, dfcinijr)

                   this%cutfcnar(nebtot)  = fcinij
                   this%cutdrvar(nebtot)  = dfcinijr

#ifdef SCREENING
                   this%cutfcnbo(nebtot)  = fcinij
                   this%cutdrvbo(nebtot)  = dfcinijr
#endif

                   this%neb(nebtot)       = j
                   this%nbb(nebtot)       = jn
#ifndef LAMMPS
                   this%dcell(nebtot)     = DCELL_INDEX(jn)
#endif

#ifdef SCREENING
                   this%sneb_seed(nebtot) = snebtot
                   this%sneb_last(nebtot) = snebtot-1

                   if (this%sneb_last(nebtot)-this%sneb_seed(nebtot)+1 > &
                        nebmax_sq) then
                      SNEB_TOO_SMALL("nebmax", i, ierror_loc)
                   endif

                   if (snebtot > snebsize) then
                      SNEB_TOO_SMALL("snebsize", i, ierror_loc)
                   endif
#endif

                   neb_last(i) = nebtot
                   nebtot      = nebtot + 1

                   if (neb_last(i)-neb_seed(i)+1 > nebmax) then
                      NEB_TOO_SMALL("nebmax", i, ierror_loc)
                   endif

                   if (nebtot > nebsize) then
                      NEB_TOO_SMALL("nebsize", i, ierror_loc)
                   endif

                endif cutoff_region

             endif j_known_el1

          enddo jn_loop1

       endif i_known_el1

    enddo i_loop1

    ! 
    ! begin potential calculation.
    !

    !$omp do
    i_loop2: do i = 1, natloc

       eli  = el(i)

       i_known_el2: if (eli > 0) then

          fi  = 0.0_DP

          istart = neb_seed(i)
          ifinsh = neb_last(i)

          !
          ! have a list of all non-negligible bonds on atom i.
          ! calculate the morse terms and derivatives.
          ! i==j loop. consider all pairs of atoms i < j.
          !

          ij_loop: do ij = istart, ifinsh
             j    = this%neb(ij)

             maskfac = 2
             if (present(mask)) then
                if (mask(i) == 0 .and. mask(j) == 0) then
                   maskfac = 0
                else if (mask(i) == 0 .or. mask(j) == 0) then
                   maskfac = 1
                endif
             endif 

#ifndef LAMMPS
             jdc  = this%dcell(ij)
#endif

             el2ij    = this%bndtyp(ij)
             rlij     = this%bndlen(ij)

             ar_within_cutoff: if (maskfac > 0 .and. &
                                   rlij < this%cut_ar_h(el2ij)) then

                fj       = 0.0_DP

                elj      = el(j)
                rlijr    = 1.0_DP / rlij
                rnij     = this%bndnm(1:3, ij)

                rij      = rlij*rnij

                ! 
                ! cutoff function and derivative.
                !

                fcarij    = this%cutfcnar(ij)
                dfcarijr  = this%cutdrvar(ij)

                ! 
                ! repulsive/attractive potentials
                !

                call VA(this, el2ij, rlij, VAij, dVAij_drij)
                call VR(this, el2ij, rlij, VRij, dVRij_drij)

                VAij       = 0.5_DP*maskfac*VAij
                dVAij_drij = 0.5_DP*maskfac*dVAij_drij
                VRij       = 0.5_DP*maskfac*VRij
                dVRij_drij = 0.5_DP*maskfac*dVRij_drij

#ifdef BO_WITH_D
                !
                ! sp-splitting terms
                !

                call D(this, eli, elj, el2ij, rlij, VAij, dVAij_drij, &
                     Dij, dDij_drij)
#endif

                !
                ! reset virial contributions
                !

                wij   = 0.0_DP
                wijb  = 0.0_DP

                ! 
                ! calculate components of bond order term and derivatives
                ! with respect to bond i-j for atom i.
                !

                zij    = 0.0_DP
                dbidi  = 0.0_DP
                dbidj  = 0.0_DP

                ! 
                ! restart i-k loop; now nci, nhi and nconj have been calculated
                !

                ikc = 0
                ik_loop2: do ik = istart, ifinsh
                   ! consider all atoms bound to i, except atom j.

                   ikc          = ikc + 1

                   k            = this%neb(ik)
#ifndef LAMMPS
                   kdc          = this%dcell(ik)
#endif

                   nebofi(ikc)  = k
#ifndef LAMMPS
                   dcofi(ikc)   = kdc
#endif
#ifdef SCREENING
                   seedi(ikc)   = this%sneb_seed(ik)
                   lasti(ikc)   = this%sneb_last(ik)
#endif

                   rlik         = this%bndlen(ik)
                   rnik         = this%bndnm(1:3, ik)

                   fcik         = this%cutfcnbo(ik)

                   ik_neq_ij: if (ik /= ij) then
                      ikpot  = this%bndtyp(ik)
                      rlik   = this%bndlen(ik)

                      bo_within_cutoff1: if (rlik < this%cut_bo_h(ikpot)) then
                         k       = this%neb(ik)
                         elk     = el(k)
                         rnik    = this%bndnm(1:3, ik)
                         rik     = rlik*rnik

                         dfcikr  = this%cutdrvbo(ik)

                         ! 
                         ! calculate length dependent terms
                         ! (constant for ijk=ccc,cch,chc).
                         !

#ifdef SEPARATE_H_ARGUMENTS
                         call h(this, elj, eli, elk, el2ij, ikpot, &
                              rlij, rlik, h_Dr, dh_dDr, dh_dr2)
#else
                         call h(this, elj, eli, elk, el2ij, ikpot, &
                              rlij - rlik, h_Dr, dh_dDr)
#endif

                         ! 
                         ! calculate angle dependent terms
                         ! ( constant for j(.)-i(c)-k(.) ).
                         !

                         !
                         ! cos( thetaijk ), g( thetaijk ) & 
                         !    dg( thetaijk ) / dcos( thetaijk )
                         ! g_costh = g( thetaijk )
                         ! dg_dcosth = dg( thetaijk ) / dcos( thetaijk )
                         !

                         costh  = dot_product(rnik, rnij)
                         call g(this, elj, eli, elk, el2ij, ikpot, costh, g_costh, dg_dcosth)

                         ! 
                         ! direction cosines of rjk = ( rik - rij ) / disjk
                         !

                         dkc     = rnik * rlik - rnij * rlij
                         disjk   = sqrt( dot_product( dkc, dkc ) )
                         dkc     = dkc / disjk

                         ! 
                         ! dcos( thetaijk ) / drwh [where w=x,y,z and h =i,j,k]
                         !

                         dcsdij  =  1.0_DP / rlik - costh * rlijr
                         dcsdik  =  rlijr      - costh / rlik
                         dcsdjk  = - disjk * rlijr / rlik

                         dcsdi   = - dcsdij*rnij - dcsdik*rnik
                         dcsdj   =   dcsdij*rnij                  - dcsdjk*dkc
                         dcsdk   =                    dcsdik*rnik + dcsdjk*dkc

                         ! 
                         ! fcik * exp * dg( thetaijk ) / drwh
                         ! [where w=x,y,z and h =i,j,k]
                         !

                         dzfac        = fcik * dg_dcosth * h_Dr

                         dgdi         = dzfac * dcsdi
                         dgdj         = dzfac * dcsdj
                         dgdk         = dzfac * dcsdk

#ifdef SCREENING

                         !
                         ! save for screening function derivative
                         !

                         zfaci(ikc)  = g_costh * h_Dr

#endif

                         ! 
                         ! sum etaij
                         !

                         zij    = zij + fcik * g_costh * h_Dr

                         !
                         ! fcik * g( thetaijk ) * dexp / drij
                         !

                         dzdrij = g_costh * fcik * dh_dDr

                         !
                         ! g( thetaijk ) * 
                         !    ( dfcik / drik * exp + fcik * dexp / drik )
                         !

#ifdef SEPARATE_H_ARGUMENTS
                         dzdrik = g_costh * ( dfcikr * h_Dr + fcik * dh_dr2)
#else
                         dzdrik = g_costh * ( dfcikr * h_Dr - fcik * dh_dDr)
#endif

                         ! 
                         ! sum detaij / drwh [where w=x,y,z and h =i,j,k]
                         ! g * ( fcik * dexp/drwi + dfcik/drwi * exp ) +
                         !     fcik * exp * dg/drwi
                         !

                         dbidi  = dbidi - dzdrij*rnij - dzdrik*rnik + dgdi

                         !
                         ! g * fcik * dexp/drwj + fcik * exp * dg/drwj
                         !

                         df     = dzdrij*rnij + dgdj
                         dbidj  = dbidj + df

                         !
                         ! g * ( fcik * dexp/drwk + dfcik/drwk * exp ) + 
                         !     fcik * exp * dg/drwk
                         !

                         dbidk(1:3, ikc)  = dzdrik*rnik + dgdk

                         !
                         ! Virial
                         !

                         wijb  = wijb &
                              - outer_product(rij, df) &
                              - outer_product(rik, dbidk(1:3, ikc))

                      else

#ifdef SCREENING

                         zfaci(ikc)     = 0.0_DP

#endif

                         dbidk(1:3, ikc)  = 0.0_DP

                      endif bo_within_cutoff1

                   endif ik_neq_ij

                enddo ik_loop2
                numnbi = ikc

                ! 
                ! bij & 0.5 * fcarij * VAij * dbij / detaij 
                ! 

#ifdef BO_WITH_D
                call bo(this, eli, el2ij, zij, Dij, &
                     bij, dbij_dzij, dbij_dDij)
                dbij_dzij = dbij_dzij * VAij * fcarij
                dbij_dDij = dbij_dDij * VAij * fcarij
#else
                call bo(this, eli, el2ij, zij, fcarij, VAij, &
                     bij, dbij_dzij)
#endif
                
                ! 
                ! now calculate the potential energies and forces for i and j.
                ! vfac   = VRij + baveij * VAij
                ! hlfvij = fcarij * vfac / 2.0
                ! 

                dffac  = 0.5_DP * fcarij * ( VRij + bij * VAij )
                pe(i)  = pe(i) + dffac
                pe(j)  = pe(j) + dffac

                if (present(epot_per_bond)) then
                   epot_per_bond(this%nbb(ij))  = epot_per_bond(this%nbb(ij)) &
                        + dffac
                endif

                !
                ! dvij / drij
                ! dffac = ( dVRij_drij + baveij * dVAij_drij ) * fcarij + 
                !     dfcarijr * vfac
                !

#ifdef BO_WITH_D
                dffac = &
                     0.5_DP * ( dVRij_drij * fcarij + &
                                bij * dVAij_drij * fcarij + &
                                VRij * dfcarijr + &
                                bij * VAij * dfcarijr ) &
                     + dDij_drij * dbij_dDij
#else
                dffac = &
                     0.5_DP * ( dVRij_drij * fcarij + &
                                bij * dVAij_drij * fcarij + &
                                VRij * dfcarijr + &
                                bij * VAij * dfcarijr )
#endif

                !
                ! compute force without bond order term
                !

                df          = dffac * rnij
                fi          = fi + df
                fj          = fj - df

                wij         = wij + outer_product(rij, df) - dbij_dzij*wijb

                if (present(f_per_bond)) then
                   f_per_bond(1:3, this%nbb(ij))  = &
                        f_per_bond(1:3, this%nbb(ij)) + df
                endif

                !
                ! compute force due to bond order term
                ! - ( dvij / drwi + 
                !     0.5 * fcarij * VAij * ( dbij / drwi + dbji / drwi )
                !

                df  = - dbij_dzij*dbidi
                fi  = fi + df

                !
                ! - ( dvij / drwj + 
                !     0.5 * fcarij * VAij * ( dbij / drwj + dbji / drwj )
                !

                df  = - dbij_dzij*dbidj
                fj  = fj + df

                !
                ! calculate forces on neighbours of i.
                !

                do ikc = 1, numnbi
                   k    = nebofi(ikc)
#ifdef LAMMPS
                   if (k /= j) then
#else
                   kdc  = dcofi(ikc)
                   if (k /= j .or. kdc /= jdc) then
#endif
                      
                      !
                      ! - ( 0.5 * fcarij * VAij * dbij / drwk ).
                      !

                      df          = - dbij_dzij * dbidk(1:3, ikc)
                      VEC3(f, k)  = VEC3(f, k) + df

#ifdef SCREENING

                      i1 = seedi(ikc)
                      i2 = lasti(ikc)
                      if (i1 <= i2) then

                         !
                         ! forces due to screening
                         !

                         this%sfacbo(i1:i2) = this%sfacbo(i1:i2) + &
                              zfaci(ikc) * dbij_dzij

                      endif

#endif

                   endif

                enddo ! ikc

#ifdef SCREENING

                !
                ! calculate forces on neighbors of i and j due to screening.
                !

                dffac = 0.5_DP * ( VRij + bij * VAij )

                do nijc = this%sneb_seed(ij), this%sneb_last(ij)

                   k = this%sneb(nijc)

#ifdef LAMMPS
                   rik  = VEC3(r, k) - VEC3(r, i)
#else
                   rik  = VEC3(r, k) - VEC3(r, i) &
                        - matmul(cell, VEC3(dc, this%sbnd(nijc))) &
                        - shear_dx*VEC(dc, this%sbnd(nijc), 3)
#endif
                   rjk  = -rij + rik

                   df   = dffac * this%cutdrarik(nijc) * rik

                   fi          = fi         + df
                   VEC3(f, k)  = VEC3(f, k) + (- df)

                   wij         = wij + outer_product(rik, df)

                   df  = dffac * this%cutdrarjk(nijc) * rjk

                   fj          = fj         + df
                   VEC3(f, k)  = VEC3(f, k) + (- df)

                   wij         = wij + outer_product(rjk, df)

                enddo ! nijc

#endif

                wpot  = wpot + wij
                if (present(wpot_per_bond)) then
                   SUM_VIRIAL(wpot_per_bond, this%nbb(ij), wij)
                endif
                if (present(wpot_per_at)) then
                   wij  = wij/2
                   SUM_VIRIAL(wpot_per_at, i, wij)
                   SUM_VIRIAL(wpot_per_at, j, wij)
                endif

                VEC3(f, j)  = VEC3(f, j) + fj

             endif ar_within_cutoff
          enddo ij_loop

          VEC3(f, i)  = VEC3(f, i) + fi

       endif  i_known_el2

    enddo i_loop2

#ifdef SCREENING

    !
    ! Restart loop, now compute forces due to screening
    !

    !$omp do
    i_loop_scr: do i = 1, natloc
       i_known_el_scr: if (el(i) > 0) then

          fi  = 0.0_DP

          ij_loop_scr: do ij = neb_seed(i), neb_last(i)
             j    = this%neb(ij)

             fj   = 0.0_DP

             wij  = 0.0_DP

             rij  = this%bndlen(ij) * this%bndnm(1:3, ij)

             istart  = this%sneb_seed(ij)
             ifinsh  = this%sneb_last(ij)

             this%cutdrboik(istart:ifinsh)  = &
                  this%sfacbo(istart:ifinsh) * this%cutdrboik(istart:ifinsh)

             this%cutdrbojk(istart:ifinsh)  = &
                  this%sfacbo(istart:ifinsh) * this%cutdrbojk(istart:ifinsh)

             nijc_loop_scr: do nijc = istart, ifinsh

                k = this%sneb(nijc)

#ifdef LAMMPS
                rik  = VEC3(r, k) - VEC3(r, i)
#else
                rik  = VEC3(r, k) - VEC3(r, i) - &
                     matmul(cell, VEC3(dc, this%sbnd(nijc))) - &
                     shear_dx*VEC(dc, this%sbnd(nijc), 3)
#endif
                rjk  = -rij + rik

                df          = this%cutdrboik(nijc) * rik

                fi          = fi         + df
                VEC3(f, k)  = VEC3(f, k) + (- df)

                wij         = wij + outer_product(rik, df)

                df          = this%cutdrbojk(nijc) * rjk

                fj          = fj         + df
                VEC3(f, k)  = VEC3(f, k) + (- df)

                wij         = wij + outer_product(rjk, df)

             enddo nijc_loop_scr

             wpot = wpot + wij
             if (present(wpot_per_bond)) then
                SUM_VIRIAL(wpot_per_bond, this%nbb(ij), wij)
             endif
             if (present(wpot_per_at)) then
                wij  = wij/2
                SUM_VIRIAL(wpot_per_at, i, wij)
                SUM_VIRIAL(wpot_per_at, j, wij)
             endif

             VEC3(f, j)  = VEC3(f, j) + fj

          enddo ij_loop_scr

          VEC3(f, i)  = VEC3(f, i) + fi

       endif i_known_el_scr
    enddo i_loop_scr

#endif

    epot  = epot + 0.5_DP*sum(pe(1:nat))

    if (present(epot_per_at)) then
       tls_sca1 = 0.5_DP*tls_sca1
       call tls_reduce(nat, sca1=epot_per_at, vec1=f_inout)
    else
       call tls_reduce(nat, vec1=f_inout)
    endif

    !$omp end parallel

#ifdef _OPENMP
    INVOKE_DELAYED_ERROR(ierror_loc, ierror)
#endif

    wpot_inout  = wpot_inout + wpot

  endsubroutine BOP_KERNEL

