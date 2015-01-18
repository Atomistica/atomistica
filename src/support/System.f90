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
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! H0 X
! H0 X   libAtoms+QUIP: atomistic simulation library
! H0 X
! H0 X   Portions of this code were written by
! H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
! H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
! H0 X
! H0 X   Copyright 2006-2010.
! H0 X
! H0 X   These portions of the source code are released under the GNU General
! H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
! H0 X
! H0 X   If you would like to license the source code under different terms,
! H0 X   please contact Gabor Csanyi, gabor@csanyi.net
! H0 X
! H0 X   Portions of this code were written by Noam Bernstein as part of
! H0 X   his employment for the U.S. Government, and are not subject
! H0 X   to copyright in the USA.
! H0 X
! H0 X
! H0 X   When using this software, please cite the following reference:
! H0 X
! H0 X   http://www.libatoms.org
! H0 X
! H0 X  Additional contributions by
! H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
! H0 X
! H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!X
!X  System module
!X  
!X  Basic system dependent functionality:
!X  
!X  mpi constants, default output objects, printing
!X  random number generators
!X 
!% The system module contains low-level routines for I/O, timing, random
!% number generation etc. The Inoutput type is used to abstract both
!% formatted and unformatted (i.e. binary) I/O.
!X
!X
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include "error.inc"

#ifndef __GFORTRAN__
#define isnan ieee_is_nan
#endif

module system_module
  use, intrinsic :: iso_c_binding

  use error_module

  implicit none

  private

  public :: DP, BOOL, C_NULL_CHAR
  integer, parameter :: DP = C_DOUBLE
  integer, parameter :: BOOL = C_BOOL

  public :: operator(//)
  interface operator(//)
     module procedure string_cat_logical, string_cat_int, string_cat_real
     module procedure string_cat_real_array, string_cat_complex
     module procedure string_cat_int_array, string_cat_logical_array
     module procedure string_cat_complex_array!, string_cat_string_array
!     module procedure logical_cat_string, logical_cat_logical, logical_cat_int, logical_cat_real
     module procedure int_cat_string!, int_cat_logical, int_cat_int, int_cat_real
     module procedure real_cat_string!, real_cat_logical, real_cat_int, real_cat_real
     module procedure real_array_cat_string
  end interface

  !% takes as arguments a default value and an optional argument, and 
  !% returns the optional argument value if it's present, otherwise
  !% the default value
  public :: optional_default
  interface optional_default
    module procedure optional_default_l, optional_default_i, optional_default_r
    module procedure optional_default_z
    module procedure optional_default_ia, optional_default_ra
  end interface optional_default

  public :: string_to_numerical
  interface string_to_numerical
     module procedure string_to_real_sub, string_to_integer_sub, string_to_logical_sub
     module procedure string_to_real1d, string_to_integer1d, string_to_logical1d
  end interface string_to_numerical

  public :: lower_case, upper_case, k_delta

  integer, save :: default_real_precision = 17

  character(kind=C_CHAR), save, target :: dummy_string(6) = [ "(","n","u","l","l",")" ]
  character(kind=C_CHAR), save, target :: one_string(2) = [ "?"," " ]

contains

  !% Convert an input string into an integer. If 'err' is present, it is set to true
  !% if an error occurred during the conversion.
  function string_to_int(string,err)
    character(*), intent(in)   :: string
    character(len=len(string)) :: local_string
    logical, optional, intent(out) :: err
    integer                    :: String_To_Int
    character(10)              :: format
    integer                    :: n
    integer stat

    local_string = adjustl(string)
    n = len_trim(local_string)
    write(format,'(a,i0,a)')'(i',n,')'
    string_to_int = 0
    read(local_string,format,iostat=stat) string_to_int
    if (present(err)) err = (stat /= 0)

  end function string_to_int

  !% Convert an input string into a logical. If 'err' is present, it is set to true
  !% if an error occurred during the conversion.
  function string_to_logical(string, err)
    character(*), intent(in)   :: string
    logical, optional, intent(out) :: err
    logical                    :: string_to_logical
    integer stat

    string_to_logical = .false.
    read(string,*,iostat=stat) string_to_logical

    if (present(err)) err = (stat /= 0)

  end function string_to_logical


  !% Convert an input string into a real. If 'err' is present, it is set to true
  !% if an error occurred during the conversion.
  function string_to_real(string, err)
    character(*), intent(in)   :: string
    logical, optional, intent(out) :: err
    real(dp)                   :: string_to_real
    integer stat

    if (present(err)) then
       err = .false.
       if (scan(adjustl(string), 'tfTF') == 1) then
          err = .true.
          return
       end if
    end if

    string_to_real = 0.0_dp
    read(string,*,iostat=stat) string_to_real

    if (present(err)) err = (stat /= 0)

  end function string_to_real

  subroutine string_to_real_sub(string,real_number,error)
     character(len=*), intent(in) :: string
     real(dp), intent(out) :: real_number
     integer, intent(out), optional :: error

     integer :: stat

     INIT_ERROR(error)

     real_number = 0.0_dp
     read(string,*,iostat=stat) real_number

     if(stat /= 0) then
        RAISE_ERROR("string_to_real_sub: could not convert, iostat="//stat, error)
     endif
  endsubroutine string_to_real_sub

  subroutine string_to_integer_sub(string,integer_number,error)
     character(len=*), intent(in) :: string
     integer, intent(out) :: integer_number
     integer, intent(out), optional :: error

     integer :: stat

     INIT_ERROR(error)

     integer_number = 0
     read(string,*,iostat=stat) integer_number

     if(stat /= 0) then
        RAISE_ERROR("string_to_integer_sub: could not convert, iostat="//stat, error)
     endif
  endsubroutine string_to_integer_sub

  subroutine string_to_logical_sub(string,logical_number,error)
     character(len=*), intent(in) :: string
     logical, intent(out) :: logical_number
     integer, intent(out), optional :: error

     integer :: stat

     INIT_ERROR(error)

     logical_number = .false.
     read(string,*,iostat=stat) logical_number

     if(stat /= 0) then
        RAISE_ERROR("string_to_logical_sub: could not convert, iostat="//stat, error)
     endif
  endsubroutine string_to_logical_sub

  subroutine string_to_real1d(string,real1d,error)
     character(len=*), intent(in) :: string
     real(dp), dimension(:), intent(out) :: real1d
     integer, intent(out), optional :: error

     integer :: stat

     INIT_ERROR(error)

     real1d = 0.0_dp
     read(string,*,iostat=stat) real1d

     if(stat /= 0) then
        RAISE_ERROR("string_to_real1d: could not convert, iostat="//stat, error)
     endif
  endsubroutine string_to_real1d

  subroutine string_to_integer1d(string,integer1d,error)
     character(len=*), intent(in) :: string
     integer, dimension(:), intent(out) :: integer1d
     integer, intent(out), optional :: error

     integer :: stat

     INIT_ERROR(error)

     integer1d = 0
     read(string,*,iostat=stat) integer1d

     if(stat /= 0) then
        RAISE_ERROR("string_to_integer1d: could not convert, iostat="//stat, error)
     endif
  endsubroutine string_to_integer1d

  subroutine string_to_logical1d(string,logical1d,error)
     character(len=*), intent(in) :: string
     logical, dimension(:), intent(out) :: logical1d
     integer, intent(out), optional :: error

     integer :: stat

     INIT_ERROR(error)

     logical1d = .false.
     read(string,*,iostat=stat) logical1d

     if(stat /= 0) then
        RAISE_ERROR("string_to_logical1d: could not convert, iostat="//stat, error)
     endif
  endsubroutine string_to_logical1d


  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  !X
  !% Concatenation functions.
  !% Overloadings for the // operator to make strings from various other types.
  !% In each case, we need to work out the exact length of the resultant string
  !% in order to avoid printing excess spaces.
  !X
  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  !>
  !! Return a string which is the real number 'r' rounded to 'digits' decimal
  !! digits
  !<
  function round(r,digits)

    real(dp), intent(in) :: r
    integer,  intent(in) :: digits
    ! below we work out the exact length of the resultant string
    !          space for '-' sign or not  + digits in integer part +                              space for . or not             + decimal digits
    character( int(0.5_dp-sign(0.5_dp,r)) + int(log10(max(1.0_dp,abs(r)+0.5_dp*10.0_dp**(-digits)))) + 1 + int(sign(0.5_dp,real(digits,dp)-0.5_dp)+0.5_dp) + max(0,digits)) :: round
    character(8) :: format

    if (digits > 0) then
       write(format,'(a,i0,a)')'(f0.',max(0,digits),')'
       write(round,format) r
    else
       write(round,'(i0)') int(r)
    end if

  end function round

  function string_cat_logical(string, log)
    character(*),      intent(in)  :: string
    logical,           intent(in)  :: log
    character((len(string)+1)) :: string_cat_logical
    write(string_cat_logical,'(a,l1)') string, log
  end function string_cat_logical

  function string_cat_logical_array(string, log)
    character(*),      intent(in)  :: string
    logical,           intent(in)  :: log(:)
    character((len(string)+2*size(log)-1)) :: string_cat_logical_array
    character(len=32) format

    format = '(a,'//size(log)//'(l1,1x),l1)'
    write(string_cat_logical_array,format) string, log
  end function string_cat_logical_array

  elemental function int_format_length(i) result(len)
    integer, intent(in)::i
    integer::len
    len = max(1,(-sign(1, i)+1)/2 + ceiling(log10(abs(real(i,dp))+0.01_dp)))
  end function int_format_length

  function string_cat_int(string, int)
    character(*),      intent(in)  :: string
    integer,           intent(in)  :: int
    ! below we work out the exact length of the resultant string
    character(len(string)+int_format_length(int)) :: string_cat_int

    write(string_cat_int,'(a,i0)') string, int
  end function string_cat_int

  function int_cat_string(int,string)
    character(*),      intent(in)  :: string
    integer,           intent(in)  :: int
    ! below we work out the exact length of the resultant string
    character(len(string)+int_format_length(int)) :: int_cat_string

    write(int_cat_string,'(i0,a)') int,string 
  end function int_cat_string

  function string_cat_int_array(string, values)
    character(*),      intent(in)  :: string
    integer,           intent(in)  :: values(:)
    ! below we work out the exact length of the resultant string
    character(len(string)+size(values)+sum(int_format_length(values)))::string_cat_int_array

      character(32) :: format

      if (size(values) == 1) then
         format = '(a,i0)'
         write(string_cat_int_array,format) string, values
      else if (size(values)>1) then
         format = '(a,' // (size(values)-1) //'(i0,1x),i0)'
         write(string_cat_int_array,format) string, values
      else
         write(string_cat_int_array,'(a)') string
      end if

  end function string_cat_int_array

  pure function real_sci_format_length() result(len)
    integer::len
    !  space sign 0.   fractional part                    E+00
    len = 1 + 1 + 2 + max(0,default_real_precision)+4
  end function real_sci_format_length


  function string_cat_real_array(string, values)
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character((len(string)+size(values)*real_sci_format_length())) :: string_cat_real_array
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(a,",I0,"e",I0,".",I0,")")') size(values), real_sci_format_length(), &
            default_real_precision
       write(string_cat_real_array, format) string, values 
    else
       write(string_cat_real_array, '(a)') string
    end if

  end function string_cat_real_array

  function string_cat_complex_array(string, values)
    character(*),      intent(in)  :: string
    complex(dp),          intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character((len(string)+2*size(values)*real_sci_format_length())) :: string_cat_complex_array
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(a,",I0,"e",I0,".",I0,")")') 2*size(values), real_sci_format_length(), &
            default_real_precision
       write(string_cat_complex_array, format) string, values 
    else
       write(string_cat_complex_array, '(a)') string
    end if

  end function string_cat_complex_array

  function string_cat_string_array(string, values)
    character(*),      intent(in)  :: string
    character(*),      intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character(len(string)+size(values)*len(values(1))) :: string_cat_string_array
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(a",I0,",",I0,"a",I0,")")') len(string), size(values)+1, len(values(1))
       write(string_cat_string_array, format) string, values 
    else
       write(string_cat_string_array, '(a)') string
    end if

  end function string_cat_string_array

  function real_array_cat_string(values, string)
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: values(:)
    ! we work out the exact length of the resultant string
    character((len(string)+size(values)*real_sci_format_length())) :: real_array_cat_string
    character(32) :: format

    if (size(values)>0) then
       ! replaced concatenation with write... for PGI bug, NB 22/6/2007
       write(format,'("(",I0,"e",I0,".",I0,",a)")') size(values), real_sci_format_length(), &
            default_real_precision
       write(real_array_cat_string, format) values, string
    else
       write(real_array_cat_string, '(a)') string
    end if

  end function real_array_cat_string

  pure function real_format_length(r) result(len)
#ifndef __GFORTRAN__
    use ieee_arithmetic
#endif
    real(dp), intent(in)::r
    integer::len

    if(isnan(r)) then
       len = 3
    else       !         sign                           int part         space?          decimal point                                                        fractional part
       len = int(0.5_dp-sign(0.5_dp,r)) + int(log10(max(1.0_dp,abs(r)))) + 1 + & 
           & int(sign(0.5_dp,real(default_real_precision,dp)-0.5_dp)+0.5_dp) &
           & + max(0,default_real_precision)

#ifdef GFORTRAN_ZERO_HACK
       !gfortran hack - 0.0000... is printed as .00000000
       if (r == 0.0) len = len - 1
#endif

    end if
  end function real_format_length

  pure function complex_format_length(c) result(len)
    complex(dp), intent(in)::c
    integer::len

    len = real_format_length(real(c))+1+real_format_length(imag(c))
  end function complex_format_length

  function real_cat_string(r, string)
#ifndef __GFORTRAN__
    use ieee_arithmetic
#endif
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: r
    ! we work out the exact length of the resultant string
    character( len(string)+real_format_length(r)) :: real_cat_string
    character(12) :: format

    if (default_real_precision > 0) then
       write(format,'(a,i0,a)')'(f0.',max(0,default_real_precision),',a)'
       if (isnan(r)) then
          write(real_cat_string,'(a,a)') "NaN", string
       else
          write(real_cat_string,format) r, string
       endif
    else
       write(real_cat_string,'(i0,a)') int(r), string
    end if
  end function real_cat_string

  function string_cat_real(string, r)
#ifndef __GFORTRAN__
    use ieee_arithmetic
#endif
    character(*),      intent(in)  :: string
    real(dp),          intent(in)  :: r
    ! we work out the exact length of the resultant string
    character( len(string)+real_format_length(r)) :: string_cat_real
    character(12) :: format

    if (default_real_precision > 0) then
       if (isnan(r)) then
	 write(string_cat_real,'(a,a)') string,"NaN"
       else
	 write(format,'(a,i0,a)')'(a,f0.',max(0,default_real_precision),')'
	 write(string_cat_real,format) string, r
       endif
    else
       write(string_cat_real,'(a,i0)') string, int(r)
    end if
  end function string_cat_real

  function string_cat_complex(string, c)
    character(*),      intent(in)  :: string
    complex(dp),          intent(in)  :: c
    ! we work out the exact length of the resultant string
    character( len(string)+complex_format_length(c)) :: string_cat_complex
    character(24) :: format

    if (default_real_precision > 0) then
       write(format,'(a,i0,a,i0,a)')'(a,f0.',max(0,default_real_precision),'," ",f0.', &
 	                                     max(0,default_real_precision),')'
       write(string_cat_complex,format) string, c
    else
       write(string_cat_complex,'(i0," ",i0)') string, int(real(c)), int(imag(c))
    end if
  end function string_cat_complex


  !% Return the mpi size and rank for the communicator 'comm'.
  !% this routine aborts of _MPI is not defined
  subroutine get_mpi_size_rank(comm, nproc, rank)
    
    integer, intent(in)  :: comm  !% MPI communicator
    integer, intent(out) :: nproc  !% Total number of processes
    integer, intent(out) :: rank  !% Rank of this process

#ifdef _MPI
    include 'mpif.h'
#endif

#ifdef _MPI
    
    integer::error_code

    call MPI_COMM_SIZE(comm, nproc, error_code)
    if (error_code .ne. MPI_SUCCESS) then
       rank=-1
       nproc=-1
       return
    endif
    call MPI_COMM_RANK(comm, rank, error_code)
    if (error_code .ne. MPI_SUCCESS) then
       rank=-1
       nproc=-1
       return
    endif
#else
    rank = 0
    nproc = 1
#endif
  end subroutine get_mpi_size_rank


  !>
  !! Take the values from 'date_and_time' and make a nice string
  !<
  function date_and_time_string(values)
    character(21)       :: date_and_time_string
    integer, intent(in) :: values(8)
    character(2)        :: time(7)
    character(4)        :: year
    integer             :: i

    write(year,'(i0)') values(1)
    do i = 2, 7
       if (i==4) cycle ! We don't use the local adjustment to UTC
       write(time(i),'(i0.2)') values(i)
    end do
    write(date_and_time_string,'(11a)') time(3),'/',time(2),'/',year,'   ',time(5),':',time(6),':',time(7)

  end function date_and_time_string


  !>
  !! Return the correct ordinal ending (st,nd,rd,th) for the given integer
  !<
  elemental function th(n)
    integer, intent(in) :: n
    character(2)        :: th
    integer             :: l,m

    l = mod(n,100)
    m = mod(n,10)

    if (l > 10 .and. l < 20) then
       th = 'th'
    else 
       select case(m)
       case(1)
          th = 'st'
       case(2)
          th = 'nd'
       case(3)
          th = 'rd'
       case default
          th = 'th'
       end select
    end if         

  end function th


  pure function optional_default_l(def, opt_val)
    logical, intent(in) :: def
    logical, intent(in), optional :: opt_val
    logical :: optional_default_l

    if (present(opt_val)) then
      optional_default_l = opt_val
    else
      optional_default_l = def
    endif

  end function optional_default_l


  pure function optional_default_i(def, opt_val)
    integer, intent(in) :: def
    integer, intent(in), optional :: opt_val
    integer :: optional_default_i

    if (present(opt_val)) then
      optional_default_i = opt_val
    else
      optional_default_i = def
    endif

  end function optional_default_i


  pure function optional_default_ia(def, opt_val)
    integer, intent(in) :: def(:)
    integer, intent(in), optional :: opt_val(size(def))
    integer :: optional_default_ia(size(def))

    if (present(opt_val)) then
      optional_default_ia = opt_val
    else
      optional_default_ia = def
    endif

  end function optional_default_ia


  pure function optional_default_r(def, opt_val)
    real(dp), intent(in) :: def
    real(dp), intent(in), optional :: opt_val
    real(dp) :: optional_default_r

    if (present(opt_val)) then
      optional_default_r = opt_val
    else
      optional_default_r = def
    endif

  end function optional_default_r


  pure function optional_default_ra(def, opt_val)
    real(dp), intent(in) :: def(:)
    real(dp), intent(in), optional :: opt_val(size(def))
    real(dp) :: optional_default_ra(size(def))

    if (present(opt_val)) then
      optional_default_ra = opt_val
    else
      optional_default_ra = def
    endif

  end function optional_default_ra


  pure function optional_default_z(def, opt_val)
    complex(dp), intent(in) :: def
    complex(dp), intent(in), optional :: opt_val
    complex(dp) :: optional_default_z

    if (present(opt_val)) then
      optional_default_z = opt_val
    else
      optional_default_z = def
    endif

  end function optional_default_z


  !% String to padded character array of length l
  function pad(s,l) result(a)
    character(len=*), intent(in) :: s
    integer, intent(in) :: l
    character(len=1), dimension(l) :: a

    integer i

    a = ' '
    do i=1,min(len(s),size(a))
       a(i) = s(i:i)
    end do
  end function pad


  !% Convert a word to upper case
  function upper_case(word)
    character(*) , intent(in) :: word
    character(len(word)) :: upper_case

    integer :: i,ic,nlen

    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= 97 .and. ic <= 122) then
          upper_case(i:i) = char(ic-32)
       else
          upper_case(i:i) = char(ic)
       end if
    end do
  end function upper_case

  !% Convert a word to lower case
  function lower_case(word)
    character(*) , intent(in) :: word
    character(len(word)) :: lower_case

    integer :: i,ic,nlen

    nlen = len(word)
    do i=1,nlen
       ic = ichar(word(i:i))
       if (ic >= 65 .and. ic <= 90) then
          lower_case(i:i) = char(ic+32)
       else
          lower_case(i:i) = char(ic)
       end if
    end do
  end function lower_case


  !>
  !! 
  !!              function k_delta
  !!
  !! returns the Kronecker delta of its arguments
  !!
  !<
  function k_delta(a,b) result( res )
    implicit none
    integer, intent(in) :: a,b
    integer :: res
    if( a==b ) then
       res = 1
    else
       res = 0
    endif
  endfunction k_delta

endmodule system_module
