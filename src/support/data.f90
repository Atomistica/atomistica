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

!>
!! General dynamic data structure
!!
!! FIXME! Add Documentation
!<
module data
  use error_module
  use system_module

  use logging

  use misc

  implicit none

  private

  public :: MAX_NAME_STR
  public :: TYPE_REAL_ATTR, TYPE_REAL3_ATTR, TYPE_REAL3x3_ATTR
  public :: TYPE_INTEGER_ATTR, TYPE_INTEGER3_ATTR
  public :: TYPE_REAL, TYPE_INTEGER, TYPE_LOGICAL
  public :: TYPE_REAL3, TYPE_REAL6, TYPE_REAL3x3

  integer, parameter                  :: MAX_NAME_STR = 100

  character(MAX_NAME_STR), parameter  :: NA_STR = "N/A"

  integer, parameter                  :: TYPE_REAL_ATTR     = 1
  integer, parameter                  :: TYPE_REAL3_ATTR    = 2
  integer, parameter                  :: TYPE_REAL3x3_ATTR  = 3
  integer, parameter                  :: TYPE_REAL          = 4
  integer, parameter                  :: TYPE_INTEGER_ATTR  = 5
  integer, parameter                  :: TYPE_INTEGER3_ATTR = 6
  integer, parameter                  :: TYPE_INTEGER       = 7
  integer, parameter                  :: TYPE_LOGICAL       = 8
  integer, parameter                  :: TYPE_REAL3         = 9
  integer, parameter                  :: TYPE_REAL6         = 10
  integer, parameter                  :: TYPE_REAL3x3       = 11

  public :: data_t
  type data_t

     integer                           :: len                          !< Length of the arrays stored in this data structure

     logical                           :: allow_def                    !< Is a change in definition allowed after the structure was allocated?

     !
     ! Attributes
     !

     integer                           :: n_real_attr
     character(MAX_NAME_STR), pointer  :: name_real_attr(:)
     real(DP), pointer                 :: data_real_attr(:)
     integer, pointer                  :: tag_real_attr(:)

     integer                           :: n_real3_attr
     character(MAX_NAME_STR), pointer  :: name_real3_attr(:)
     real(DP), pointer                 :: data_real3_attr(:, :)
     integer, pointer                  :: tag_real3_attr(:)

     integer                           :: n_real3x3_attr
     character(MAX_NAME_STR), pointer  :: name_real3x3_attr(:)
     real(DP), pointer                 :: data_real3x3_attr(:, :, :)
     integer, pointer                  :: tag_real3x3_attr(:)

     integer                           :: n_integer_attr
     character(MAX_NAME_STR), pointer  :: name_integer_attr(:)
     integer, pointer                  :: data_integer_attr(:)
     integer, pointer                  :: tag_integer_attr(:)

     integer                           :: n_integer3_attr
     character(MAX_NAME_STR), pointer  :: name_integer3_attr(:)
     integer, pointer                  :: data_integer3_attr(:, :)
     integer, pointer                  :: tag_integer3_attr(:)

     !
     ! Fields
     !

     integer                           :: n_real
     character(MAX_NAME_STR), pointer  :: name_real(:)
     real(DP), pointer                 :: data_real(:, :)
     integer, pointer                  :: tag_real(:)
     character(MAX_NAME_STR), pointer  :: unit_real(:)
     real(DP), pointer                 :: conv_real(:)

     integer                           :: n_integer
     character(MAX_NAME_STR), pointer  :: name_integer(:)
     character(MAX_NAME_STR), pointer  :: alias_integer(:)
     integer, pointer                  :: data_integer(:, :)
     integer, pointer                  :: tag_integer(:)

     integer                           :: n_logical
     character(MAX_NAME_STR), pointer  :: name_logical(:)
     logical, pointer                  :: data_logical(:, :)
     integer, pointer                  :: tag_logical(:)

     integer                           :: n_real3
     character(MAX_NAME_STR), pointer  :: name_real3(:)
     real(DP), pointer                 :: data_real3(:, :, :)
     integer, pointer                  :: tag_real3(:)
     character(MAX_NAME_STR), pointer  :: unit_real3(:)
     real(DP), pointer                 :: conv_real3(:)

     integer                           :: n_real6
     character(MAX_NAME_STR), pointer  :: name_real6(:)
     real(DP), pointer                 :: data_real6(:, :, :)
     integer, pointer                  :: tag_real6(:)
     character(MAX_NAME_STR), pointer  :: unit_real6(:)
     real(DP), pointer                 :: conv_real6(:)

     integer                           :: n_real3x3
     character(MAX_NAME_STR), pointer  :: name_real3x3(:)
     real(DP), pointer                 :: data_real3x3(:, :, :, :)
     integer, pointer                  :: tag_real3x3(:)
     character(MAX_NAME_STR), pointer  :: unit_real3x3(:)
     real(DP), pointer                 :: conv_real3x3(:)

  endtype data_t

  public :: init
  interface init
     module procedure data_init, data_init_from_data
  endinterface

  public :: allocate
  interface allocate
     module procedure data_allocate
  endinterface

  public :: allocated
  interface allocated
     module procedure data_allocated
  endinterface

  public :: del
  interface del
     module procedure data_del
  endinterface

  public :: exists
  interface exists
     module procedure data_exists
  endinterface

  public :: add_real_attr
  interface add_real_attr
     module procedure data_add_real_attr
  endinterface

  public :: add_real3_attr
  interface add_real3_attr
     module procedure data_add_real3_attr
  endinterface

  public :: add_real3x3_attr
  interface add_real3x3_attr
     module procedure data_add_real3x3_attr
  endinterface

  public :: add_integer_attr
  interface add_integer_attr
     module procedure data_add_integer_attr
  endinterface

  public :: add_integer3_attr
  interface add_integer3_attr
     module procedure data_add_integer3_attr
  endinterface

  public :: add_real
  interface add_real
     module procedure data_add_real
  endinterface

  public :: add_integer
  interface add_integer
     module procedure data_add_integer
  endinterface

  public :: add_logical
  interface add_logical
     module procedure data_add_logical
  endinterface

  public :: add_real3
  interface add_real3
     module procedure data_add_real3
  endinterface

  public :: add_real6
  interface add_real6
     module procedure data_add_real6
  endinterface

  public :: add_real3x3
  interface add_real3x3
     module procedure data_add_real3x3
  endinterface

  public :: copy
  interface copy
     module procedure data_copy, data_copy_from_data, data_copy_slice_from_data
  endinterface

  public :: attr_by_name
  interface attr_by_name
     module procedure data_ptr_by_name_real_attr
     module procedure data_ptr_by_name_real3_attr
     module procedure data_ptr_by_name_real3x3_attr
     module procedure data_ptr_by_name_integer_attr
     module procedure data_ptr_by_name_integer3_attr
  endinterface

  public :: ptr_by_name
  interface ptr_by_name
     module procedure data_ptr_by_name_real
     module procedure data_ptr_by_name_integer
     module procedure data_ptr_by_name_logical
     module procedure data_ptr_by_name_realX
     module procedure data_ptr_by_name_realXxX
  endinterface

  public :: tag_by_name
  interface tag_by_name
     module procedure data_tag_by_name
  endinterface

  public :: print_to_log
  interface print_to_log
     module procedure data_print_to_log
  endinterface

  public :: set_tag_by_name
  interface set_tag_by_name
     module procedure data_set_tag_by_name
  endinterface

  public :: swap
  interface swap
     module procedure data_swap
  endinterface

#ifdef _MPI

  public :: size_by_tag
  interface size_by_tag
     module procedure data_size_by_tag
  endinterface

  public :: pack_buffer
  interface pack_buffer
     module procedure data_pack_buffer
  endinterface

  public :: unpack_buffer
  interface unpack_buffer
     module procedure data_unpack_buffer
  endinterface

#endif

  public :: index_by_name

contains

  !**********************************************************************
  ! Constructor
  !********************************************************************** 
  elemental subroutine data_init(this)
    implicit none

    type(data_t), intent(inout)  :: this

    ! ---

    this%len                = -1

    this%allow_def          = .false.

    this%n_real_attr        = 0
    this%n_real3_attr       = 0
    this%n_real3x3_attr     = 0
    this%n_integer_attr     = 0
    this%n_integer3_attr     = 0

    this%name_real_attr     => NULL()
    this%name_real3_attr    => NULL()
    this%name_real3x3_attr  => NULL()
    this%name_integer_attr  => NULL()
    this%name_integer3_attr => NULL()

    this%data_real_attr     => NULL()
    this%data_real3_attr    => NULL()
    this%data_real3x3_attr  => NULL()
    this%data_integer_attr  => NULL()
    this%data_integer3_attr => NULL()

    this%tag_real_attr      => NULL()
    this%tag_real3_attr     => NULL()
    this%tag_real3x3_attr   => NULL()
    this%tag_integer_attr   => NULL()
    this%tag_integer3_attr  => NULL()

    this%n_real             = 0
    this%n_integer          = 0
    this%n_logical          = 0
    this%n_real3            = 0
    this%n_real6            = 0
    this%n_real3x3          = 0

    this%name_real          => NULL()
    this%name_integer       => NULL()
    this%name_logical       => NULL()
    this%name_real3         => NULL()
    this%name_real6         => NULL()
    this%name_real3x3       => NULL()

    this%data_real          => NULL()
    this%data_integer       => NULL()
    this%data_logical       => NULL()
    this%data_real3         => NULL()
    this%data_real6         => NULL()
    this%data_real3x3       => NULL()

    this%tag_real           => NULL()
    this%tag_integer        => NULL()
    this%tag_logical        => NULL()
    this%tag_real3          => NULL()
    this%tag_real6          => NULL()
    this%tag_real3x3        => NULL()

    this%alias_integer      => NULL()

    this%unit_real          => NULL()
    this%unit_real3         => NULL()
    this%unit_real6         => NULL()
    this%unit_real3x3       => NULL()

    this%conv_real          => NULL()
    this%conv_real3         => NULL()
    this%conv_real6         => NULL()
    this%conv_real3x3       => NULL()

  endsubroutine data_init


  !**********************************************************************
  ! Constructor
  !********************************************************************** 
  subroutine data_init_from_data(this, from)
    implicit none

    type(data_t), intent(inout)  :: this
    type(data_t), intent(in)     :: from

    ! ---

    this%len                = -1

    this%allow_def          = .false.

    this%n_real_attr        = from%n_real_attr
    this%n_real3_attr       = from%n_real3_attr
    this%n_real3x3_attr     = from%n_real3x3_attr
    this%n_integer_attr     = from%n_integer_attr
    this%n_integer3_attr    = from%n_integer3_attr

    this%name_real_attr     => NULL()
    this%name_real3_attr    => NULL()
    this%name_real3x3_attr  => NULL()
    this%name_integer_attr  => NULL()
    this%name_integer3_attr => NULL()

    this%data_real_attr     => NULL()
    this%data_real3_attr    => NULL()
    this%data_real3x3_attr  => NULL()
    this%data_integer_attr  => NULL()
    this%data_integer3_attr => NULL()

    this%tag_real_attr     => NULL()
    this%tag_real3_attr    => NULL()
    this%tag_real3x3_attr  => NULL()
    this%tag_integer_attr  => NULL()
    this%tag_integer3_attr => NULL()

    if (this%n_real_attr > 0) then
       allocate(this%name_real_attr(this%n_real_attr))
       this%name_real_attr(:)  = from%name_real_attr(:)
       this%tag_real_attr(:)   = from%tag_real_attr(:)
    endif
    if (this%n_real3_attr > 0) then
       allocate(this%name_real3_attr(this%n_real3_attr))
       this%name_real3_attr(:)  = from%name_real3_attr(:)
       this%tag_real3_attr(:)   = from%tag_real3_attr(:)
    endif
    if (this%n_real3x3_attr > 0) then
       allocate(this%name_real3x3_attr(this%n_real3x3_attr))
       this%name_real3x3_attr(:)  = from%name_real3x3_attr(:)
       this%tag_real3x3_attr(:)   = from%tag_real3x3_attr(:)
    endif
    if (this%n_integer_attr > 0) then
       allocate(this%name_integer_attr(this%n_integer_attr))
       this%name_integer_attr(:)  = from%name_integer_attr(:)
       this%tag_integer_attr(:)   = from%tag_integer_attr(:)
    endif
    if (this%n_integer3_attr > 0) then
       allocate(this%name_integer3_attr(this%n_integer3_attr))
       this%name_integer3_attr(:)  = from%name_integer3_attr(:)
       this%tag_integer3_attr(:)   = from%tag_integer3_attr(:)
    endif


    this%n_real             = from%n_real
    this%n_integer          = from%n_integer
    this%n_logical          = from%n_logical
    this%n_real3            = from%n_real3
    this%n_real6            = from%n_real6
    this%n_real3x3          = from%n_real3x3

    this%name_real          => NULL()
    this%name_integer       => NULL()
    this%name_logical       => NULL()
    this%name_real3         => NULL()
    this%name_real6         => NULL()
    this%name_real3x3       => NULL()

    this%data_real          => NULL()
    this%data_integer       => NULL()
    this%data_logical       => NULL()
    this%data_real3         => NULL()
    this%data_real6         => NULL()
    this%data_real3x3       => NULL()

    this%tag_real           => NULL()
    this%tag_integer        => NULL()
    this%tag_logical        => NULL()
    this%tag_real3          => NULL()
    this%tag_real6          => NULL()
    this%tag_real3x3        => NULL()

    this%alias_integer      => NULL()

    this%unit_real          => NULL()
    this%unit_real3         => NULL()
    this%unit_real6         => NULL()
    this%unit_real3x3       => NULL()

    this%conv_real          => NULL()
    this%conv_real3         => NULL()
    this%conv_real6         => NULL()
    this%conv_real3x3       => NULL()

    if (this%n_real > 0) then
       allocate(this%name_real(this%n_real))
       allocate(this%tag_real(this%n_real))
       allocate(this%unit_real(this%n_real))
       allocate(this%conv_real(this%n_real))

       this%name_real(:)  = from%name_real(:)
       this%tag_real(:)   = from%tag_real(:)
       this%unit_real(:)  = from%unit_real(:)
       this%conv_real(:)  = from%conv_real(:)
    endif
    if (this%n_integer > 0) then
       allocate(this%name_integer(this%n_integer))
       allocate(this%tag_integer(this%n_integer))
       allocate(this%alias_integer(this%n_integer))

       this%name_integer(:)  = from%name_integer(:)
       this%tag_integer(:)   = from%tag_integer(:)
       this%alias_integer    = from%alias_integer
    endif
    if (this%n_logical > 0) then
       allocate(this%name_logical(this%n_logical))
       allocate(this%tag_logical(this%n_logical))

       this%name_logical(:)  = from%name_logical(:)
       this%tag_logical(:)   = from%tag_logical(:)
    endif
    if (this%n_real3 > 0) then
       allocate(this%name_real3(this%n_real3))
       allocate(this%tag_real3(this%n_real3))
       allocate(this%unit_real3(this%n_real3))
       allocate(this%conv_real3(this%n_real3))

       this%name_real3(:)  = from%name_real3(:)
       this%tag_real3(:)   = from%tag_real3(:)
       this%unit_real3(:)  = from%unit_real3(:)
       this%conv_real3(:)  = from%conv_real3(:)
    endif
    if (this%n_real6 > 0) then
       allocate(this%name_real6(this%n_real6))
       allocate(this%tag_real6(this%n_real6))
       allocate(this%unit_real6(this%n_real6))
       allocate(this%conv_real6(this%n_real6))

       this%name_real6(:)  = from%name_real6(:)
       this%tag_real6(:)   = from%tag_real6(:)
       this%unit_real6(:)  = from%unit_real6(:)
       this%conv_real6(:)  = from%conv_real6(:)
    endif
    if (this%n_real3x3 > 0) then
       allocate(this%name_real3x3(this%n_real3x3))
       allocate(this%tag_real3x3(this%n_real3x3))
       allocate(this%unit_real3x3(this%n_real3x3))
       allocate(this%conv_real3x3(this%n_real3x3))

       this%name_real3x3(:)  = from%name_real3x3(:)
       this%tag_real3x3(:)   = from%tag_real3x3(:)
       this%unit_real3x3(:)  = from%unit_real3x3(:)
       this%conv_real3x3(:)  = from%conv_real3x3(:)
    endif

  endsubroutine data_init_from_data


  !**********************************************************************
  ! Constructor
  !********************************************************************** 
  elemental subroutine data_allocate(this, len, allow_def)
    implicit none

    type(data_t), intent(inout)    :: this
    integer, intent(in)            :: len
    logical, intent(in), optional  :: allow_def

    ! ---

    this%len  = len

    if (present(allow_def)) then
       this%allow_def  = allow_def
    endif

    if (this%n_real_attr > 0) then
       allocate(this%data_real_attr(this%n_real_attr))
       this%data_real_attr(:)  = 0.0_DP
    endif

    if (this%n_real3_attr > 0) then
       allocate(this%data_real3_attr(3, this%n_real3_attr))
       this%data_real3_attr(:, :)  = 0.0_DP
    endif

    if (this%n_real3x3_attr > 0) then
       allocate(this%data_real3x3_attr(3, 3, this%n_real3x3_attr))
       this%data_real3x3_attr(:, :, :)  = 0.0_DP
    endif

    if (this%n_integer_attr > 0) then
       allocate(this%data_integer_attr(this%n_integer_attr))
       this%data_integer_attr(:)  = 0
    endif

    if (this%n_integer3_attr > 0) then
       allocate(this%data_integer3_attr(3, this%n_integer3_attr))
       this%data_integer3_attr(:, :)  = 0
    endif

    if (this%n_real > 0) then
       allocate(this%data_real(this%len, this%n_real))
       this%data_real(:, :)  = 0.0_DP
    endif

    if (this%n_integer > 0) then
       allocate(this%data_integer(this%len, this%n_integer))
       this%data_integer(:, :)  = 0
    endif

    if (this%n_logical > 0) then
       allocate(this%data_logical(this%len, this%n_logical))
       this%data_logical(:, :)  = .false.
    endif

    if (this%n_real3 > 0) then
#ifdef __SEP_XYZ__
       allocate(this%data_real3(this%len, 3, this%n_real3))
#else
       allocate(this%data_real3(3, this%len, this%n_real3))
#endif
       this%data_real3(:, :, :)  = 0.0_DP
    endif

    if (this%n_real6 > 0) then
#ifdef __SEP_XYZ__
       allocate(this%data_real6(this%len, 6, this%n_real6))
#else
       allocate(this%data_real6(6, this%len, this%n_real6))
#endif
       this%data_real6(:, :, :)  = 0.0_DP
    endif

    if (this%n_real3x3 > 0) then
       allocate(this%data_real3x3(3, 3, this%len, this%n_real3x3))
       this%data_real3x3(:, :, :, :)  = 0.0_DP
    endif

  endsubroutine data_allocate


  !**********************************************************************
  ! Check if the object has been allocated
  !********************************************************************** 
  function data_allocated(this)
    implicit none

    type(data_t), intent(in)  :: this
    logical                   :: data_allocated

    ! ---

    data_allocated  = this%len > 0

  endfunction data_allocated


  !**********************************************************************
  ! Destructor
  !********************************************************************** 
  elemental subroutine data_del(this)
    implicit none

    type(data_t), intent(inout)  :: this

    ! ---

    if (associated(this%data_real))     deallocate(this%data_real)
    if (associated(this%data_integer))  deallocate(this%data_integer)
    if (associated(this%data_logical))  deallocate(this%data_logical)
    if (associated(this%data_real3))    deallocate(this%data_real3)
    if (associated(this%data_real6))    deallocate(this%data_real6)
    if (associated(this%data_real3x3))  deallocate(this%data_real3x3)

    if (associated(this%name_real)) then
       deallocate(this%name_real)
       deallocate(this%tag_real)
       deallocate(this%unit_real)
       deallocate(this%conv_real)
    endif

    if (associated(this%name_integer)) then
       deallocate(this%name_integer)
       deallocate(this%tag_integer)
       deallocate(this%alias_integer)
    endif

    if (associated(this%name_logical)) then
       deallocate(this%name_logical)
       deallocate(this%tag_logical)
    endif

    if (associated(this%name_real3)) then
       deallocate(this%name_real3)
       deallocate(this%tag_real3)
       deallocate(this%unit_real3)
       deallocate(this%conv_real3)
    endif

    if (associated(this%name_real6)) then
       deallocate(this%name_real6)
       deallocate(this%tag_real6)
       deallocate(this%unit_real6)
       deallocate(this%conv_real6)
    endif

    if (associated(this%name_real3x3)) then
       deallocate(this%name_real3x3)
       deallocate(this%tag_real3x3)
       deallocate(this%unit_real3x3)
       deallocate(this%conv_real3x3)
    endif

  endsubroutine data_del


  !**********************************************************************
  ! Find the index for a certain name
  !********************************************************************** 
  pure function index_by_name(n, names, name)
    implicit none

    integer, intent(in)                  :: n
    character(MAX_NAME_STR), intent(in)  :: names(n)
    character(*), intent(in)             :: name

    integer                              :: index_by_name

    ! ---

    integer  :: i, j

    ! ---

!    write (*, *)  "name = ", name

    j = -1
    do i = 1, n
!       write (*, '(A,A,A,A,A)')  "'", trim(name), "' = '", trim(names(i)), "'?"

       if (equal(name, names(i))) then
!          write (*, *)  "yes"
          j = i
       endif
    enddo

    index_by_name = j

  endfunction index_by_name


  !**********************************************************************
  ! Check if a field with this name already exists
  !********************************************************************** 
  function data_exists(this, name, data_type)
    implicit none

    type(data_t), intent(in)        :: this
    character(*), intent(in)        :: name
    integer, intent(out), optional  :: data_type

    logical                         :: data_exists

    ! ---

    data_exists = .false.

    if (associated(this%name_real_attr)) then
       if ( index_by_name(this%n_real_attr, this%name_real_attr, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL_ATTR
          endif
       endif
    endif

    if (associated(this%name_real3_attr)) then
       if ( index_by_name(this%n_real3_attr, this%name_real3_attr, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL3_ATTR
          endif
       endif
    endif

    if (associated(this%name_real3x3_attr)) then
       if ( index_by_name(this%n_real3x3_attr, this%name_real3x3_attr, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL3x3_ATTR
          endif
       endif
    endif

    if (associated(this%name_integer_attr)) then
       if ( index_by_name(this%n_integer_attr, this%name_integer_attr, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_INTEGER_ATTR
          endif
       endif
    endif

    if (associated(this%name_integer3_attr)) then
       if ( index_by_name(this%n_integer3_attr, this%name_integer3_attr, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_INTEGER3_ATTR
          endif
       endif
    endif

    if (associated(this%name_real)) then
       if ( index_by_name(this%n_real, this%name_real, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL
          endif
       endif
    endif

    if (associated(this%name_integer)) then
       if ( index_by_name(this%n_integer, this%name_integer, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_INTEGER
          endif
       endif
    endif

    if (associated(this%name_logical)) then
       if ( index_by_name(this%n_logical, this%name_logical, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_LOGICAL
          endif
       endif
    endif

    if (associated(this%name_real3)) then
       if ( index_by_name(this%n_real3, this%name_real3, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL3
          endif
       endif
    endif

    if (associated(this%name_real6)) then
       if ( index_by_name(this%n_real6, this%name_real6, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL6
          endif
       endif
    endif

    if (associated(this%name_real3x3)) then
       if ( index_by_name(this%n_real3x3, this%name_real3x3, name) > 0 ) then
          data_exists = .true.
          if (present(data_type)) then
             data_type  = TYPE_REAL3x3
          endif
       endif
    endif

  endfunction data_exists


  !**********************************************************************
  ! Check if a field with this name already exists,
  ! bail out if it doesn't
  !********************************************************************** 
  subroutine data_name_check(this, name, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    integer, intent(inout), optional  :: ierror

    ! ---

    if (data_exists(this, name)) then
       RAISE_ERROR("Field '" // trim(name) // "' already exists.", ierror)
    endif

  endsubroutine data_name_check


  !**********************************************************************
  ! Add a new (real) attribute
  !********************************************************************** 
  subroutine data_add_real_attr(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in), optional     :: tag
    integer, intent(inout), optional  :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:)
    integer, pointer                  :: old_tag(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real_attr
    old_data => this%data_real_attr
    old_tag  => this%tag_real_attr

    this%n_real_attr  = this%n_real_attr + 1

    allocate(this%name_real_attr(this%n_real_attr))
    allocate(this%tag_real_attr(this%n_real_attr))

    this%name_real_attr(this%n_real_attr)    = name
    this%tag_real_attr(this%n_real_attr)     = 0
    if (present(tag)) then
       this%tag_real_attr(this%n_real_attr)  = tag
    endif

    if (associated(old_name)) then
       this%name_real_attr(1:this%n_real_attr-1)  = old_name(1:this%n_real_attr-1)
       deallocate(old_name)
    endif
    if (associated(old_tag)) then
       this%tag_real_attr(1:this%n_real_attr-1)  = old_tag(1:this%n_real_attr-1)
       deallocate(old_tag)
    endif

    if (this%len > 0) then

       allocate(this%data_real_attr(this%n_real_attr))

       this%data_real_attr(this%n_real_attr)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real_attr(1:this%n_real_attr-1)  = old_data(1:this%n_real_attr-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_real_attr


  !**********************************************************************
  ! Add a new (real3) attribute
  !********************************************************************** 
  subroutine data_add_real3_attr(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in), optional     :: tag
    integer, intent(inout), optional  :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:, :)
    integer, pointer                  :: old_tag(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real3_attr
    old_data => this%data_real3_attr
    old_tag  => this%tag_real3_attr

    this%n_real3_attr  = this%n_real3_attr + 1

    allocate(this%name_real3_attr(this%n_real3_attr))
    allocate(this%tag_real3_attr(this%n_real3_attr))

    this%name_real3_attr(this%n_real3_attr)    = name
    this%tag_real3_attr(this%n_real3_attr)     = 0
    if (present(tag)) then
       this%tag_real3_attr(this%n_real3_attr)  = tag
    endif

    if (associated(old_name)) then
       this%name_real3_attr(1:this%n_real3_attr-1)  = old_name(1:this%n_real3_attr-1)
       deallocate(old_name)
    endif
    if (associated(old_tag)) then
       this%tag_real3_attr(1:this%n_real3_attr-1)  = old_tag(1:this%n_real3_attr-1)
       deallocate(old_tag)
    endif

    if (this%len > 0) then

       allocate(this%data_real3_attr(3, this%n_real3_attr))

       this%data_real3_attr(1:3, this%n_real3_attr)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real3_attr(1:3, 1:this%n_real3_attr-1)  = old_data(1:3, 1:this%n_real3_attr-1)
          deallocate(old_data)
       endif

    endif

 endsubroutine data_add_real3_attr


  !**********************************************************************
  ! Add a new (real3) attribute
  !********************************************************************** 
  subroutine data_add_real3x3_attr(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in), optional     :: tag
    integer, intent(inout), optional  :: ierror


    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:, :, :)
    integer, pointer                  :: old_tag(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real3x3_attr
    old_data => this%data_real3x3_attr
    old_tag  => this%tag_real3x3_attr

    this%n_real3x3_attr  = this%n_real3x3_attr + 1

    allocate(this%name_real3x3_attr(this%n_real3x3_attr))
    allocate(this%tag_real3x3_attr(this%n_real3x3_attr))

    this%name_real3x3_attr(this%n_real3x3_attr)    = name
    this%tag_real3x3_attr(this%n_real3x3_attr)     = 0
    if (present(tag)) then
       this%tag_real3x3_attr(this%n_real3x3_attr)  = tag
    endif

    if (associated(old_name)) then
       this%name_real3x3_attr(1:this%n_real3x3_attr-1)  = old_name(1:this%n_real3x3_attr-1)
       deallocate(old_name)
    endif
    if (associated(old_tag)) then
       this%tag_real3x3_attr(1:this%n_real3x3_attr-1)  = old_tag(1:this%n_real3x3_attr-1)
       deallocate(old_tag)
    endif

    if (this%len > 0) then

       allocate(this%data_real3x3_attr(3, 3, this%n_real3x3_attr))

       this%data_real3x3_attr(1:3, 1:3, this%n_real3x3_attr)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real3x3_attr(1:3, 1:3, 1:this%n_real3x3_attr-1)  = old_data(1:3, 1:3, 1:this%n_real3x3_attr-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_real3x3_attr


  !**********************************************************************
  ! Add a new (integer) attribute
  !********************************************************************** 
  subroutine data_add_integer_attr(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in), optional     :: tag
    integer, intent(inout), optional  :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    integer, pointer                  :: old_data(:)
    integer, pointer                  :: old_tag(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_integer_attr
    old_data => this%data_integer_attr
    old_tag  => this%tag_integer_attr

    this%n_integer_attr  = this%n_integer_attr + 1

    allocate(this%name_integer_attr(this%n_integer_attr))
    allocate(this%tag_integer_attr(this%n_integer_attr))

    this%name_integer_attr(this%n_integer_attr)    = name
    this%tag_integer_attr(this%n_integer_attr)     = 0
    if (present(tag)) then
       this%tag_integer_attr(this%n_integer_attr)  = tag
    endif

    if (associated(old_name)) then
       this%name_integer_attr(1:this%n_integer_attr-1)  = old_name(1:this%n_integer_attr-1)
       deallocate(old_name)
    endif
    if (associated(old_tag)) then
       this%tag_integer_attr(1:this%n_integer_attr-1)  = old_tag(1:this%n_integer_attr-1)
       deallocate(old_tag)
    endif

    if (this%len > 0) then

       allocate(this%data_integer_attr(this%n_integer_attr))

       this%data_integer_attr(this%n_integer_attr)  = 0

       if (associated(old_data)) then
          this%data_integer_attr(1:this%n_integer_attr-1)  = old_data(1:this%n_integer_attr-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_integer_attr


  !**********************************************************************
  ! Add a new (integer3) attribute
  !********************************************************************** 
  subroutine data_add_integer3_attr(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in), optional     :: tag
    integer, intent(inout), optional  :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    integer, pointer                  :: old_data(:, :)
    integer, pointer                  :: old_tag(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_integer3_attr
    old_data => this%data_integer3_attr
    old_tag  => this%tag_integer3_attr

    this%n_integer3_attr  = this%n_integer3_attr + 1

    allocate(this%name_integer3_attr(this%n_integer3_attr))
    allocate(this%tag_integer3_attr(this%n_integer3_attr))

    this%name_integer3_attr(this%n_integer3_attr)    = name
    this%tag_integer3_attr(this%n_integer3_attr)     = 0
    if (present(tag)) then
       this%tag_integer3_attr(this%n_integer3_attr)  = tag
    endif

    if (associated(old_name)) then
       this%name_integer3_attr(1:this%n_integer3_attr-1)  = old_name(1:this%n_integer3_attr-1)
       deallocate(old_name)
    endif
    if (associated(old_tag)) then
       this%tag_integer3_attr(1:this%n_integer3_attr-1)  = old_tag(1:this%n_integer3_attr-1)
       deallocate(old_tag)
    endif

    if (this%len > 0) then

       allocate(this%data_integer3_attr(3, this%n_integer3_attr))

       this%data_integer3_attr(:, this%n_integer3_attr)  = 0

       if (associated(old_data)) then
          this%data_integer3_attr(:, 1:this%n_integer3_attr-1)  = old_data(:, 1:this%n_integer3_attr-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_integer3_attr


  !**********************************************************************
  ! Add a new (real) field
  !********************************************************************** 
  subroutine data_add_real(this, name, tag, unit, conv, ierror)
    implicit none

    type(data_t), intent(inout)         :: this
    character(*), intent(in)            :: name
    integer, intent(in), optional       :: tag
    character(*), intent(in), optional  :: unit
    real(DP), intent(in), optional      :: conv
    integer, intent(inout), optional    :: ierror


    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:, :)
    integer, pointer                  :: old_tag(:)
    character(MAX_NAME_STR), pointer  :: old_unit(:)
    real(DP), pointer                 :: old_conv(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real
    old_data => this%data_real
    old_tag  => this%tag_real
    old_unit => this%unit_real
    old_conv => this%conv_real

    this%n_real  = this%n_real + 1

    allocate(this%name_real(this%n_real))
    allocate(this%tag_real(this%n_real))
    allocate(this%unit_real(this%n_real))
    allocate(this%conv_real(this%n_real))

    this%name_real(this%n_real)     = name
    this%tag_real(this%n_real)      = 0
    this%unit_real(this%n_real)     = "1"
    this%conv_real(this%n_real)     = 1.0_DP
    if (present(tag)) then
       this%tag_real(this%n_real)   = tag
    endif
    if (present(unit)) then
       this%unit_real(this%n_real)  = unit
    else
       this%unit_real(this%n_real)  = NA_STR
    endif
    if (present(conv)) then
       this%conv_real(this%n_real)  = conv
    endif

    if (associated(old_name)) then
       this%name_real(1:this%n_real-1)     = old_name(1:this%n_real-1)
       this%tag_real(1:this%n_real-1)      = old_tag(1:this%n_real-1)
       this%unit_real(1:this%n_real-1)     = old_unit(1:this%n_real-1)
       this%conv_real(1:this%n_real-1)     = old_conv(1:this%n_real-1)
       deallocate(old_name)
       deallocate(old_tag)
       deallocate(old_unit)
       deallocate(old_conv)
    endif

    if (this%len > 0) then

       allocate(this%data_real(this%len, this%n_real))

       this%data_real(:, this%n_real)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real(:, 1:this%n_real-1)  = old_data(:, 1:this%n_real-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_real


  !**********************************************************************
  ! Add a new (integer) field
  !********************************************************************** 
  subroutine data_add_integer(this, name, tag, alias, ierror)
    implicit none

    type(data_t),           intent(inout) :: this
    character(*),           intent(in)    :: name
    integer,      optional, intent(in)    :: tag
    character(*), optional, intent(in)    :: alias
    integer,      optional, intent(inout) :: ierror

    ! ---

    character(MAX_NAME_STR), pointer :: old_name(:)
    integer,                 pointer :: old_data(:, :)
    integer,                 pointer :: old_tag(:)
    character(MAX_NAME_STR), pointer :: old_alias(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_integer
    old_data => this%data_integer
    old_tag => this%tag_integer
    old_alias => this%alias_integer

    this%n_integer = this%n_integer + 1

    allocate(this%name_integer(this%n_integer))
    allocate(this%tag_integer(this%n_integer))
    allocate(this%alias_integer(this%n_integer))

    this%name_integer(this%n_integer)     = name
    this%tag_integer(this%n_integer)      = 0
    if (present(tag)) then
       this%tag_integer(this%n_integer)   = tag
    endif
    if (present(alias)) then
       this%alias_integer(this%n_integer) = alias
    else
       this%alias_integer(this%n_integer) = "*"
    endif

    if (associated(old_name)) then
       this%name_integer(1:this%n_integer-1)  = old_name(1:this%n_integer-1)
       this%tag_integer(1:this%n_integer-1)   = old_tag(1:this%n_integer-1)
       this%alias_integer(1:this%n_integer-1) = old_alias(1:this%n_integer-1)
       deallocate(old_name)
       deallocate(old_tag)
       deallocate(old_alias)
    endif

    if (this%len > 0) then

       allocate(this%data_integer(this%len, this%n_integer))

       this%data_integer(:, this%n_integer)  = 0

       if (associated(old_data)) then
          this%data_integer(:, 1:this%n_integer-1)  = old_data(:, 1:this%n_integer-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_integer


  !**********************************************************************
  ! Add a new (logical) field
  !********************************************************************** 
  subroutine data_add_logical(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in), optional     :: tag
    integer, intent(inout), optional  :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    logical, pointer                  :: old_data(:, :)
    integer, pointer                  :: old_tag(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_logical
    old_data => this%data_logical
    old_tag  => this%tag_logical

    this%n_logical  = this%n_logical + 1

    allocate(this%name_logical(this%n_logical))
    allocate(this%tag_logical(this%n_logical))

    this%name_logical(this%n_logical)     = name
    this%tag_logical(this%n_logical)      = 0
    if (present(tag)) then
       this%tag_logical(this%n_logical)   = tag
    endif

    if (associated(old_name)) then
       this%name_logical(1:this%n_logical-1)     = old_name(1:this%n_logical-1)
       this%tag_logical(1:this%n_logical-1)      = old_tag(1:this%n_logical-1)
       deallocate(old_name)
       deallocate(old_tag)
    endif

    if (this%len > 0) then

       allocate(this%data_logical(this%len, this%n_logical))

       this%data_logical(:, this%n_logical)  = .false.

       if (associated(old_data)) then
          this%data_logical(:, 1:this%n_logical-1)  = old_data(:, 1:this%n_logical-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_logical


  !**********************************************************************
  ! Add a new (real3) field
  !********************************************************************** 
  subroutine data_add_real3(this, name, tag, unit, conv, ierror)
    implicit none

    type(data_t), intent(inout)         :: this
    character(*), intent(in)            :: name
    integer, intent(in), optional       :: tag
    character(*), intent(in), optional  :: unit
    real(DP), intent(in), optional      :: conv
    integer, intent(inout), optional    :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:, :, :)
    integer, pointer                  :: old_tag(:)
    character(MAX_NAME_STR), pointer  :: old_unit(:)
    real(DP), pointer                 :: old_conv(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real3
    old_data => this%data_real3
    old_tag  => this%tag_real3
    old_unit => this%unit_real3
    old_conv => this%conv_real3

    this%n_real3  = this%n_real3 + 1

    allocate(this%name_real3(this%n_real3))
    allocate(this%tag_real3(this%n_real3))
    allocate(this%unit_real3(this%n_real3))
    allocate(this%conv_real3(this%n_real3))

    this%name_real3(this%n_real3)        = name
    this%tag_real3(this%n_real3)         = 0
    this%unit_real3(this%n_real3)        = "1"
    this%conv_real3(this%n_real3)        = 1.0_DP
    if (present(tag)) then
       this%tag_real3(this%n_real3)      = tag
    endif
    if (present(unit)) then
       this%unit_real3(this%n_real3)     = unit
    else
       this%unit_real3(this%n_real3)     = NA_STR
    endif
    if (present(conv)) then
       this%conv_real3(this%n_real3)     = conv
    endif

    if (associated(old_name)) then
       this%name_real3(1:this%n_real3-1)        = old_name(1:this%n_real3-1)
       this%tag_real3(1:this%n_real3-1)         = old_tag(1:this%n_real3-1)
       this%unit_real3(1:this%n_real3-1)        = old_unit(1:this%n_real3-1)
       this%conv_real3(1:this%n_real3-1)        = old_conv(1:this%n_real3-1)
       deallocate(old_name)
       deallocate(old_tag)
       deallocate(old_unit)
       deallocate(old_conv)
    endif

    if (this%len > 0) then

#ifdef __SEP_XYZ__
       allocate(this%data_real3(this%len, 3, this%n_real3))
#else
       allocate(this%data_real3(3, this%len, this%n_real3))
#endif

       this%data_real3(:, :, this%n_real3)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real3(:, :, 1:this%n_real3-1)  = old_data(:, :, 1:this%n_real3-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_real3


  !**********************************************************************
  ! Add a new (real6) field
  !********************************************************************** 
  subroutine data_add_real6(this, name, tag, unit, conv, ierror)
    implicit none

    type(data_t), intent(inout)         :: this
    character(*), intent(in)            :: name
    integer, intent(in), optional       :: tag
    character(*), intent(in), optional  :: unit
    real(DP), intent(in), optional      :: conv
    integer, intent(inout), optional    :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:, :, :)
    integer, pointer                  :: old_tag(:)
    character(MAX_NAME_STR), pointer  :: old_unit(:)
    real(DP), pointer                 :: old_conv(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real6
    old_data => this%data_real6
    old_tag  => this%tag_real6
    old_unit => this%unit_real6
    old_conv => this%conv_real6

    this%n_real6  = this%n_real6 + 1

    allocate(this%name_real6(this%n_real6))
    allocate(this%tag_real6(this%n_real6))
    allocate(this%unit_real6(this%n_real6))
    allocate(this%conv_real6(this%n_real6))

    this%name_real6(this%n_real6)        = name
    this%tag_real6(this%n_real6)         = 0
    this%unit_real6(this%n_real6)        = "1"
    this%conv_real6(this%n_real6)        = 1.0_DP
    if (present(tag)) then
       this%tag_real6(this%n_real6)      = tag
    endif
    if (present(unit)) then
       this%unit_real6(this%n_real6)     = unit
    else
       this%unit_real6(this%n_real6)     = NA_STR
    endif
    if (present(conv)) then
       this%conv_real6(this%n_real6)     = conv
    endif

    if (associated(old_name)) then
       this%name_real6(1:this%n_real6-1)        = old_name(1:this%n_real6-1)
       this%tag_real6(1:this%n_real6-1)         = old_tag(1:this%n_real6-1)
       this%unit_real6(1:this%n_real6-1)        = old_unit(1:this%n_real6-1)
       this%conv_real6(1:this%n_real6-1)        = old_conv(1:this%n_real6-1)
       deallocate(old_name)
       deallocate(old_tag)
       deallocate(old_unit)
       deallocate(old_conv)
    endif

    if (this%len > 0) then

#ifdef __SEP_XYZ__
       allocate(this%data_real6(this%len, 6, this%n_real6))
#else
       allocate(this%data_real6(6, this%len, this%n_real6))
#endif

       this%data_real6(:, :, this%n_real6)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real6(:, :, 1:this%n_real6-1)  = old_data(:, :, 1:this%n_real6-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_real6


  !**********************************************************************
  ! Add a new (real3x3 - a tensor) field
  !********************************************************************** 
  subroutine data_add_real3x3(this, name, tag, unit, conv, ierror)
    implicit none

    type(data_t), intent(inout)         :: this
    character(*), intent(in)            :: name
    integer, intent(in), optional       :: tag
    character(*), intent(in), optional  :: unit
    real(DP), intent(in), optional      :: conv
    integer, intent(inout), optional    :: ierror

    ! ---

    character(MAX_NAME_STR), pointer  :: old_name(:)
    real(DP), pointer                 :: old_data(:, :, :, :)
    integer, pointer                  :: old_tag(:)
    character(MAX_NAME_STR), pointer  :: old_unit(:)
    real(DP), pointer                 :: old_conv(:)

    ! ---

    if (this%len > 0 .and. .not. this%allow_def) then
       RAISE_ERROR("Cannot modify the data structure after it was initialized.", ierror)
    endif

    call data_name_check(this, name)

    old_name => this%name_real3x3
    old_data => this%data_real3x3
    old_tag  => this%tag_real3x3
    old_unit => this%unit_real3x3
    old_conv => this%conv_real3x3

    this%n_real3x3  = this%n_real3x3 + 1

    allocate(this%name_real3x3(this%n_real3x3))
    allocate(this%tag_real3x3(this%n_real3x3))
    allocate(this%unit_real3x3(this%n_real3x3))
    allocate(this%conv_real3x3(this%n_real3x3))

    this%name_real3x3(this%n_real3x3)           = name
    this%tag_real3x3(this%n_real3x3)            = 0
    this%unit_real3x3(this%n_real3x3)           = "1"
    this%conv_real3x3(this%n_real3x3)           = 1.0_DP
    if (present(tag)) then
       this%tag_real3x3(this%n_real3x3)         = tag
   endif
    if (present(unit)) then
       this%unit_real3x3(this%n_real3x3)        = unit
    else
       this%unit_real3x3(this%n_real3x3)        = NA_STR
    endif
    if (present(conv)) then
       this%conv_real3x3(this%n_real3x3)        = conv
    endif

    if (associated(old_name)) then
       this%name_real3x3(1:this%n_real3x3-1)           = old_name(1:this%n_real3x3-1)
       this%tag_real3x3(1:this%n_real3x3-1)            = old_tag(1:this%n_real3x3-1)
       this%unit_real3x3(1:this%n_real3x3-1)           = old_unit(1:this%n_real3x3-1)
       this%conv_real3x3(1:this%n_real3x3-1)           = old_conv(1:this%n_real3x3-1)
       deallocate(old_name)
       deallocate(old_tag)
       deallocate(old_unit)
       deallocate(old_conv)
    endif

    if (this%len > 0) then

       allocate(this%data_real3x3(3, 3, this%len, this%n_real3x3))
    
       this%data_real3x3(:, :, :, this%n_real3x3)  = 0.0_DP

       if (associated(old_data)) then
          this%data_real3x3(:, :, :, 1:this%n_real3x3-1)  = old_data(:, :, :, 1:this%n_real3x3-1)
          deallocate(old_data)
       endif

    endif

  endsubroutine data_add_real3x3


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_real_attr(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    real(DP), pointer, intent(out)    :: ptr
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real_attr, this%name_real_attr(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown real attribute: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_real_attr(i)

  endsubroutine data_ptr_by_name_real_attr


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_real3_attr(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    real(DP), pointer, intent(out)    :: ptr(:)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real3_attr, this%name_real3_attr(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown real3 attribute: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_real3_attr(:, i)

  endsubroutine data_ptr_by_name_real3_attr


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_real3x3_attr(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    real(DP), pointer, intent(out)    :: ptr(:, :)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real3x3_attr, this%name_real3x3_attr(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown real3x3 attribute: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_real3x3_attr(:, :, i)

  endsubroutine data_ptr_by_name_real3x3_attr


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_integer_attr(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    integer, pointer, intent(out)     :: ptr
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_integer_attr, this%name_integer_attr(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown integer attribute: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_integer_attr(i)

  endsubroutine data_ptr_by_name_integer_attr


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_integer3_attr(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    integer, pointer, intent(out)     :: ptr(:)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_integer3_attr, this%name_integer3_attr(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown integer attribute: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_integer3_attr(:, i)

  endsubroutine data_ptr_by_name_integer3_attr


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_real(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    real(DP), pointer, intent(out)    :: ptr(:)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real, this%name_real(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown real field: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_real(:, i)

  endsubroutine data_ptr_by_name_real


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_integer(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    integer, pointer, intent(out)     :: ptr(:)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_integer, this%name_integer(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown integer field: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_integer(:, i)

  endsubroutine data_ptr_by_name_integer


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_logical(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    logical, pointer, intent(out)     :: ptr(:)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_logical, this%name_logical(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown logical field: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_logical(:, i)

  endsubroutine data_ptr_by_name_logical


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_realX(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    real(DP), pointer, intent(out)    :: ptr(:, :)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real3, this%name_real3, name)

    if (i < 0) then
       i = index_by_name(this%n_real6, this%name_real6, name)

       if (i < 0) then
          RAISE_ERROR("Unknown real3 field: '" // trim(name) // "'.", ierror)
       else
          ptr => this%data_real6(:, :, i)
       endif
    else
       ptr => this%data_real3(:, :, i)
    endif

  endsubroutine data_ptr_by_name_realX


  !**********************************************************************
  ! Return a pointer to the field data
  !********************************************************************** 
  subroutine data_ptr_by_name_realXxX(this, name, ptr, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    real(DP), pointer, intent(out)    :: ptr(:, :, :)
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real3x3, this%name_real3x3(:), name)

    if (i < 0) then
       RAISE_ERROR("Unknown real3x3 field: '" // trim(name) // "'.", ierror)
    endif

    ptr => this%data_real3x3(:, :, :, i)

  endsubroutine data_ptr_by_name_realXxX


  !>
  !! Set the tag
  !!
  !! Set the tag given the name of a certain field
  !<
  subroutine data_set_tag_by_name(this, name, tag, ierror)
    implicit none

    type(data_t), intent(inout)       :: this
    character(*), intent(in)          :: name
    integer, intent(in)               :: tag
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real, this%name_real(:), name)
    if (i > 0) then
       this%tag_real(i)  = tag
    else
       i = index_by_name(this%n_integer, this%name_integer(:), name)
       if (i > 0) then
          this%tag_integer(i)  = tag
       else
          i = index_by_name(this%n_logical, this%name_logical(:), name)
          if (i > 0) then
             this%tag_logical(i)  = tag
          else
             i = index_by_name(this%n_real3, this%name_real3(:), name)
             if (i > 0) then
                this%tag_real3(i)  = tag
             else
                i = index_by_name(this%n_real3x3, this%name_real3x3(:), name)
                if (i > 0) then
                   this%tag_real3x3(i)  = tag
                else
                   RAISE_ERROR("Unknown field: '" // trim(name) // "'.", ierror)
                endif
             endif
          endif
       endif
    endif

  endsubroutine data_set_tag_by_name


  !>
  !! Return the tag
  !!
  !! Return the tag given the name of a certain field
  !<
  function data_tag_by_name(this, name, ierror)
    implicit none

    type(data_t), intent(in)          :: this
    character(*), intent(in)          :: name
    integer, intent(inout), optional  :: ierror

    integer                   :: data_tag_by_name

    ! ---

    integer  :: i

    ! ---

    i = index_by_name(this%n_real, this%name_real(:), name)
    if (i > 0) then
       data_tag_by_name  = this%tag_real(i)
    else
       i = index_by_name(this%n_integer, this%name_integer(:), name)
       if (i > 0) then
          data_tag_by_name  = this%tag_integer(i)
       else
          i = index_by_name(this%n_logical, this%name_logical(:), name)
          if (i > 0) then
             data_tag_by_name  = this%tag_logical(i)
          else
             i = index_by_name(this%n_real3, this%name_real3(:), name)
             if (i > 0) then
                data_tag_by_name  = this%tag_real3(i)
             else
                i = index_by_name(this%n_real3x3, this%name_real3x3(:), name)
                if (i > 0) then
                   data_tag_by_name  = this%tag_real3x3(i)
                else
                   RAISE_ERROR("Unknown field: '" // trim(name) // "'.", ierror)
                endif
             endif
          endif
       endif
    endif

  endfunction data_tag_by_name


  !**********************************************************************
  ! Copy entry with index *s* to entry with index *t*
  !********************************************************************** 
  elemental subroutine data_copy(this, t, s)
    implicit none

    type(data_t), intent(inout)  :: this
    integer, intent(in)          :: t
    integer, intent(in)          :: s

    ! ---

    if (this%n_real > 0) then
       this%data_real(t, :)           = this%data_real(s, :)
    endif
    if (this%n_integer > 0) then
       this%data_integer(t, :)        = this%data_integer(s, :)
    endif
    if (this%n_logical > 0) then
       this%data_logical(t, :)        = this%data_logical(s, :)
    endif
    if (this%n_real3 > 0) then
#ifdef __SEP_XYZ__
       this%data_real3(t, :, :)       = this%data_real3(s, :, :)
#else
       this%data_real3(:, t, :)       = this%data_real3(:, s, :)
#endif
    endif
    if (this%n_real3x3 > 0) then
       this%data_real3x3(:, :, t, :)  = this%data_real3x3(:, :, s, :)
    endif

  endsubroutine data_copy


  !**********************************************************************
  ! Copy entry with index *s* of *from* to entry with index *t*
  !********************************************************************** 
  elemental subroutine data_copy_from_data(this, t, from, s)
    implicit none

    type(data_t), intent(inout)  :: this
    integer, intent(in)          :: t
    type(data_t), intent(in)     :: from
    integer, intent(in)          :: s

    ! ---

    integer  :: i, j

    ! ---

    if (this%n_real > 0) then
       do i = 1, this%n_real
          j = index_by_name(from%n_real, from%name_real(:), this%name_real(i))
          if (j > 0) then
             this%data_real(t, i)  = from%data_real(s, j)
          endif
       enddo
    endif
    if (this%n_integer > 0) then
       do i = 1, this%n_integer
          j = index_by_name(from%n_integer, from%name_integer(:), this%name_integer(i))
          if (j > 0) then
             this%data_integer(t, i)  = from%data_integer(s, j)
          endif
       enddo
    endif
    if (this%n_logical > 0) then
       do i = 1, this%n_logical
          j = index_by_name(from%n_logical, from%name_logical(:), this%name_logical(i))
          if (j > 0) then
             this%data_logical(t, i)  = from%data_logical(s, j)
          endif
       enddo
    endif
    if (this%n_real3 > 0) then
       do i = 1, this%n_real3
          j = index_by_name(from%n_real3, from%name_real3(:), this%name_real3(i))
          if (j > 0) then
#ifdef __SEP_XYZ__
             this%data_real3(t, :, i)  = from%data_real3(s, :, j)
#else
             this%data_real3(:, t, i)  = from%data_real3(:, s, j)
#endif
          endif
       enddo
    endif
    if (this%n_real3x3 > 0) then
       do i = 1, this%n_real3x3
          j = index_by_name(from%n_real3x3, from%name_real3x3(:), this%name_real3x3(i))
          if (j > 0) then
             this%data_real3x3(:, :, t, i)  = from%data_real3x3(:, :, s, j)
          endif
       enddo
    endif

  endsubroutine data_copy_from_data


  !**********************************************************************
  ! Copy entry with index *s* of *from* to entry with index *t*
  !********************************************************************** 
  elemental subroutine data_copy_slice_from_data(this, t1, t2, from, s1, s2)
    implicit none

    type(data_t), intent(inout)  :: this
    integer, intent(in)          :: t1
    integer, intent(in)          :: t2
    type(data_t), intent(in)     :: from
    integer, intent(in)          :: s1
    integer, intent(in)          :: s2

    ! ---

    integer  :: i, j

    ! ---

    !assert t2-t1 == s2-s1

    if (this%n_real > 0) then
       do i = 1, this%n_real
          j = index_by_name(from%n_real, from%name_real(:), this%name_real(i))
          if (j > 0) then
             this%data_real(t1:t2, i)  = from%data_real(s1:s2, j)
          endif
       enddo
    endif
    if (this%n_integer > 0) then
       do i = 1, this%n_integer
          j = index_by_name(from%n_integer, from%name_integer(:), this%name_integer(i))
          if (j > 0) then
             this%data_integer(t1:t2, i)  = from%data_integer(s1:s2, j)
          endif
       enddo
    endif
    if (this%n_logical > 0) then
       do i = 1, this%n_logical
          j = index_by_name(from%n_logical, from%name_logical(:), this%name_logical(i))
          if (j > 0) then
             this%data_logical(t1:t2, i)  = from%data_logical(s1:s2, j)
          endif
       enddo
    endif
    if (this%n_real3 > 0) then
       do i = 1, this%n_real3
          j = index_by_name(from%n_real3, from%name_real3(:), this%name_real3(i))
          if (j > 0) then
#ifdef __SEP_XYZ__
             this%data_real3(t1:t2, :, i)  = from%data_real3(s1:s2, :, j)
#else
             this%data_real3(:, t1:t2, i)  = from%data_real3(:, s1:s2, j)
#endif
          endif
       enddo
    endif
    if (this%n_real3x3 > 0) then
       do i = 1, this%n_real3x3
          j = index_by_name(from%n_real3x3, from%name_real3x3(:), this%name_real3x3(i))
          if (j > 0) then
             this%data_real3x3(:, :, t1:t2, i)  = from%data_real3x3(:, :, s1:s2, j)
          endif
       enddo
    endif

  endsubroutine data_copy_slice_from_data


  !**********************************************************************
  ! Swap entry with index *s* with entry with index *t*
  !********************************************************************** 
  elemental subroutine data_swap(this, i1, i2)
    implicit none

    type(data_t), intent(inout)  :: this
    integer, intent(in)          :: i1
    integer, intent(in)          :: i2

    ! ---

    integer  :: j

    ! ---

    if (this%n_real > 0) then
       do j = 1, this%n_real
          call swap(this%data_real(i1, j), this%data_real(i2, j))
       enddo
    endif
    if (this%n_integer > 0) then
       do j = 1, this%n_integer
          call swap(this%data_integer(i1, j), this%data_integer(i2, j))
       enddo
    endif
    if (this%n_logical > 0) then
       do j = 1, this%n_logical
          call swap(this%data_logical(i1, j), this%data_logical(i2, j))
       enddo
    endif
    if (this%n_real3 > 0) then
       do j = 1, this%n_real3
#ifdef __SEP_XYZ__
          call swap(this%data_real3(i1, 1, j), this%data_real3(i2, 1, j))
          call swap(this%data_real3(i1, 2, j), this%data_real3(i2, 2, j))
          call swap(this%data_real3(i1, 3, j), this%data_real3(i2, 3, j))
#else
          call swap(this%data_real3(1, i1, j), this%data_real3(1, i2, j))
          call swap(this%data_real3(2, i1, j), this%data_real3(2, i2, j))
          call swap(this%data_real3(3, i1, j), this%data_real3(3, i2, j))
#endif
       enddo
    endif
    if (this%n_real3x3 > 0) then
       do j = 1, this%n_real3x3
          call swap(this%data_real3x3(:, :, i1, j), this%data_real3x3(:, :, i2, j))
       enddo
    endif

  endsubroutine data_swap


  !**********************************************************************
  ! Output logging information
  !********************************************************************** 
  subroutine data_print_to_log(this)
    implicit none

    type(data_t), intent(in)  :: this

    ! ---

    integer  :: i

    ! ---

    call prlog("- data_print_to_log -")

    call log_memory_start("data_print_to_log")

    do i = 1, this%n_real
       call prlog("     real     :: "//trim(this%name_real(i)))

       call log_memory_estimate(this%data_real)
       call log_memory_estimate(this%tag_real)
       call log_memory_estimate(this%conv_real)
    enddo

    do i = 1, this%n_integer
       call prlog("     integer  :: "//trim(this%name_integer(i)))

       call log_memory_estimate(this%data_integer)
       call log_memory_estimate(this%tag_integer)
    enddo

    do i = 1, this%n_logical
       call prlog("     logical  :: "//trim(this%name_logical(i)))
     
       call log_memory_estimate(this%data_logical)
       call log_memory_estimate(this%tag_logical)
    enddo

    do i = 1, this%n_real3
       call prlog("     real3    :: "//trim(this%name_real3(i)))

       call log_memory_estimate(this%data_real3)
       call log_memory_estimate(this%tag_real3)
       call log_memory_estimate(this%conv_real3)
    enddo

    do i = 1, this%n_real3x3
       call prlog("     real3x3  :: "//trim(this%name_real3x3(i)))
     
       call log_memory_estimate(this%data_real3x3)
       call log_memory_estimate(this%tag_real3x3)
       call log_memory_estimate(this%conv_real3x3)
    enddo

    call log_memory_stop("data_print_to_log")

    call prlog

  endsubroutine data_print_to_log


#ifdef _MPI

  !**********************************************************************
  ! Size of a data block containing all fields with a certain tag
  !**********************************************************************
  subroutine data_size_by_tag(this, tag, size)
    implicit none

    type(data_t), intent(in)  :: this
    integer, intent(in)       :: tag
    integer, intent(out)      :: size

    ! ---

    integer  :: i

    ! ---

    size  = 0

    do i = 1, this%n_real
       if (iand(this%tag_real(i), tag) /= 0) then
          size  = size + 1
       endif
    enddo

    do i = 1, this%n_integer
       if (iand(this%tag_integer(i), tag) /= 0) then
          size  = size + 1
       endif
    enddo

    do i = 1, this%n_logical
       if (iand(this%tag_logical(i), tag) /= 0) then
          size  = size + 1
       endif
    enddo

    do i = 1, this%n_real3
       if (iand(this%tag_real3(i), tag) /= 0) then
          size  = size + 3
       endif
    enddo

    do i = 1, this%n_real3x3
       if (iand(this%tag_real3x3(i), tag) /= 0) then
          size  = size + 9
       endif
    enddo

  endsubroutine data_size_by_tag


  !**********************************************************************
  ! Pack buffer, copy all data entries with tag *tag* to the buffer
  !**********************************************************************
  subroutine data_pack_buffer(this, tag, data_i, buffer_i, buffer)
    implicit none

    type(data_t), intent(in)  :: this
    integer, intent(in)       :: tag
    integer, intent(in)       :: data_i
    integer, intent(inout)    :: buffer_i
    real(DP), intent(inout)   :: buffer(:)

    ! ---

    integer  :: i

    ! ---

    do i = 1, this%n_real
       if (iand(this%tag_real(i), tag) /= 0) then
          buffer_i          = buffer_i + 1
          buffer(buffer_i)  = this%data_real(data_i, i)
       endif
    enddo

    do i = 1, this%n_integer
       if (iand(this%tag_integer(i), tag) /= 0) then
          buffer_i          = buffer_i + 1
          buffer(buffer_i)  = this%data_integer(data_i, i)
       endif
    enddo

    do i = 1, this%n_logical
       if (iand(this%tag_logical(i), tag) /= 0) then
          buffer_i          = buffer_i + 1
          if (this%data_logical(data_i, i)) then
             buffer(buffer_i)  = 1.0_DP
          else
             buffer(buffer_i)  = 0.0_DP
          endif
       endif
    enddo

    do i = 1, this%n_real3
       if (iand(this%tag_real3(i), tag) /= 0) then
#ifdef __SEP_XYZ__
          buffer(buffer_i+1)  = this%data_real3(data_i, 1, i)
          buffer(buffer_i+2)  = this%data_real3(data_i, 2, i)
          buffer(buffer_i+3)  = this%data_real3(data_i, 3, i)
#else
          buffer(buffer_i+1)  = this%data_real3(1, data_i, i)
          buffer(buffer_i+2)  = this%data_real3(2, data_i, i)
          buffer(buffer_i+3)  = this%data_real3(3, data_i, i)
#endif
          buffer_i            = buffer_i + 3
       endif
    enddo

    do i = 1, this%n_real3x3
       if (iand(this%tag_real3x3(i), tag) /= 0) then
          buffer(buffer_i+1)  = this%data_real3x3(1, 1, data_i, i)
          buffer(buffer_i+2)  = this%data_real3x3(2, 1, data_i, i)
          buffer(buffer_i+3)  = this%data_real3x3(3, 1, data_i, i)
          buffer(buffer_i+4)  = this%data_real3x3(1, 2, data_i, i)
          buffer(buffer_i+5)  = this%data_real3x3(2, 2, data_i, i)
          buffer(buffer_i+6)  = this%data_real3x3(3, 2, data_i, i)
          buffer(buffer_i+7)  = this%data_real3x3(1, 3, data_i, i)
          buffer(buffer_i+8)  = this%data_real3x3(2, 3, data_i, i)
          buffer(buffer_i+9)  = this%data_real3x3(3, 3, data_i, i)

          buffer_i            = buffer_i + 9
       endif
    enddo

  endsubroutine data_pack_buffer


  !**********************************************************************
  ! Unpack buffer
  !**********************************************************************
  subroutine data_unpack_buffer(this, tag, buffer_i, buffer, data_i)
    implicit none

    type(data_t), intent(inout)  :: this
    integer, intent(in)          :: tag
    integer, intent(inout)       :: buffer_i
    real(DP), intent(in)         :: buffer(:)
    integer, intent(in)          :: data_i

    ! ---

    integer  :: i

    ! ---

    do i = 1, this%n_real
       if (iand(this%tag_real(i), tag) /= 0) then
          buffer_i                   = buffer_i + 1
          this%data_real(data_i, i)  = buffer(buffer_i)
       endif
    enddo

    do i = 1, this%n_integer
       if (iand(this%tag_integer(i), tag) /= 0) then
          buffer_i                      = buffer_i + 1
          this%data_integer(data_i, i)  = buffer(buffer_i)
       endif
    enddo

    do i = 1, this%n_logical
       if (iand(this%tag_logical(i), tag) /= 0) then
          buffer_i                      = buffer_i + 1
          this%data_logical(data_i, i)  = abs(buffer(buffer_i)) > 1e-6_DP
       endif
    enddo

    do i = 1, this%n_real3
       if (iand(this%tag_real3(i), tag) /= 0) then
#ifdef __SEP_XYZ__
          this%data_real3(data_i, 1, i)  = buffer(buffer_i+1)
          this%data_real3(data_i, 2, i)  = buffer(buffer_i+2)
          this%data_real3(data_i, 3, i)  = buffer(buffer_i+3)
#else
          this%data_real3(1, data_i, i)  = buffer(buffer_i+1)    
          this%data_real3(2, data_i, i)  = buffer(buffer_i+2)
          this%data_real3(3, data_i, i)  = buffer(buffer_i+3)
#endif
          buffer_i                       = buffer_i + 3
       endif
    enddo

    do i = 1, this%n_real3x3
       if (iand(this%tag_real3x3(i), tag) /= 0) then
          this%data_real3x3(1, 1, data_i, i)  = buffer(buffer_i+1)    
          this%data_real3x3(2, 1, data_i, i)  = buffer(buffer_i+2)
          this%data_real3x3(3, 1, data_i, i)  = buffer(buffer_i+3)
          this%data_real3x3(1, 2, data_i, i)  = buffer(buffer_i+4)    
          this%data_real3x3(2, 2, data_i, i)  = buffer(buffer_i+5)
          this%data_real3x3(3, 2, data_i, i)  = buffer(buffer_i+6)
          this%data_real3x3(1, 3, data_i, i)  = buffer(buffer_i+7)    
          this%data_real3x3(2, 3, data_i, i)  = buffer(buffer_i+8)
          this%data_real3x3(3, 3, data_i, i)  = buffer(buffer_i+9)

          buffer_i                            = buffer_i + 9
       endif
    enddo

  endsubroutine data_unpack_buffer

#endif

endmodule data
