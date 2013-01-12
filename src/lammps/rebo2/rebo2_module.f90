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
  !>
  !! Alloc memory and return pointer to a new potential object
  !<
  subroutine CREATE_FUNC(this_cptr) bind(C)
    use, intrinsic  :: iso_c_binding

    implicit none

    type(C_PTR),    intent(out)  :: this_cptr

    ! ---

    type(BOP_TYPE), pointer  :: this

    allocate(this)
    this_cptr = c_loc(this)

  endsubroutine CREATE_FUNC


  !>
  !! Free memory
  !<
  subroutine DESTROY_FUNC(this_cptr) bind(C)
    use, intrinsic  :: iso_c_binding

    implicit none

    type(C_PTR), value  :: this_cptr

    ! ---

    type(BOP_TYPE), pointer  :: this

    call c_f_pointer(this_cptr, this)
    deallocate(this)

  endsubroutine DESTROY_FUNC


  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine INIT_FUNC(this_cptr) bind(C)
    use, intrinsic  :: iso_c_binding

    implicit none

    type(C_PTR), value ::  this_cptr

    ! ---
    
    type(BOP_TYPE), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    call rebo2_default_Fcc_table(this%in_Fcc, this%in_dFdi, this%in_dFdj, this%in_dFdk)
    call rebo2_default_Fch_table(this%in_Fch)
    call rebo2_default_Fhh_table(this%in_Fhh)
    call rebo2_default_Pcc_table(this%in_Pcc)
    call rebo2_default_Pch_table(this%in_Pch)
    call rebo2_default_Tcc_table(this%in_Tcc)

    call rebo2_db_init_with_parameters( &
         this, &
         this%in_Fcc, this%in_dFdi, this%in_dFdj, this%in_dFdk, &
         this%in_Fch, this%in_Fhh, &
         this%in_Pcc, this%in_Pch, &
         this%in_Tcc)

  endsubroutine INIT_FUNC


  !>
  !! Destructor
  !!
  !!
  !<
  subroutine DEL_FUNC(this_cptr) bind(C)
    use, intrinsic  :: iso_c_binding

    implicit none

    type(C_PTR), value  :: this_cptr

    ! ---

    type(BOP_TYPE), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    call rebo2_db_del(this)

    if (allocated(this%internal_el)) then
       deallocate(this%internal_el)
    endif

  endsubroutine DEL_FUNC


  !>
  !! Return values of cutoff
  !<
  subroutine GET_CUTOFF_FUNC(this_cptr, rcmaxsq) bind(C)
    use, intrinsic  :: iso_c_binding

    implicit none

    type(C_PTR),    value        :: this_cptr
    real(C_DOUBLE), intent(out)  :: rcmaxsq(2, 2)
    
    ! ---

    type(BOP_TYPE), pointer  :: this

    ! ---

    call c_f_pointer(this_cptr, this)

    rcmaxsq(1, 1)  = this%max_cut_sq(C_C)
    rcmaxsq(2, 2)  = this%max_cut_sq(H_H)
    rcmaxsq(2, 1)  = this%max_cut_sq(C_H)
    rcmaxsq(1, 2)  = this%max_cut_sq(C_H)

  endsubroutine GET_CUTOFF_FUNC


  !>
  !! Compute the force
  !!
  !! Compute the force
  !<
  subroutine COMPUTE_FUNC(this_cptr, natloc, nat, tag, Z, r, seed, last, nneighb, &
       neighb, epot, f, wpot, wpot_at, ierror) bind(C)
    use, intrinsic  :: iso_c_binding

    implicit none

    type(C_PTR),              value          :: this_cptr
    integer(C_INT),           value          :: natloc
    integer(C_INT),           value          :: nat
    integer(C_INT),           intent(in)     :: tag(nat)
    integer(C_INT),           intent(in)     :: Z(nat)
    real(C_DOUBLE),           intent(in)     :: r(3, nat)
    integer(C_INT),           intent(in)     :: seed(nat+1)
    integer(C_INT),           intent(in)     :: last(nat+1)
    integer(C_INT),           value          :: nneighb
    integer(C_INT),           intent(in)     :: neighb(nneighb)
    real(C_DOUBLE),           intent(inout)  :: epot
    real(C_DOUBLE),           intent(inout)  :: f(3, nat)
    real(C_DOUBLE),           intent(inout)  :: wpot(3, 3)
    real(C_DOUBLE),           intent(inout)  :: wpot_at(6, nat)
    integer(C_INT), optional, intent(inout)  :: ierror

    ! ---

    type(BOP_TYPE), pointer  :: this
    integer                  :: i

    ! ---

    call c_f_pointer(this_cptr, this)

    if (allocated(this%internal_el)) then
       if (size(this%internal_el) < nat) then
          deallocate(this%internal_el)
       endif
    endif

    if (.not. allocated(this%internal_el)) then
       allocate(this%internal_el(nat))
    endif

    this%internal_el = 0
    do i = 1, nat
       if (Z(i) == 1 .or. Z(i) == 3) then
          this%internal_el(i) = rebo2_C_
       else if (Z(i) == 2) then
          this%internal_el(i) = rebo2_H_
       endif
    enddo

    call BOP_KERNEL( &
         this, &
         nat, natloc, nat, r, &
         tag, this%internal_el, &
         seed, last, neighb, nneighb, &
         epot, f, wpot, wpot_per_at=wpot_at, &
         ierror=ierror)
    PASS_ERROR(ierror)

  endsubroutine COMPUTE_FUNC

