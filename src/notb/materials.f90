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
!   shared
! @endmeta

!>
!! Materials database for the tight-binding module.
!!
!! Materials database for the tight-binding module. Contains routines to read
!! Pekka's database format and Frauenheims/Elstners .skf format.
!! (See http://www.dftb.org/ for the latter.)
!<

#include "macros.inc"

module materials
  use, intrinsic  :: iso_c_binding

  use supplib

  use io
  use logging
  use nonuniform_spline
  use misc

  use particles

  implicit none

  private

  public :: MAX_NORB
  integer, parameter  :: MAX_NORB = 10

  ! Notation for the orbital-integrals:
  !  dds ddp ddd pds pdp pps ppp sds sps sss
  !   1   2   3   4   5   6   7   8   9   10
  public :: O_dds, O_ddp, O_ddd, O_pds, O_pdp, O_pps, O_ppp, O_sds, O_sps, O_sss
  integer, parameter :: O_dds = 1
  integer, parameter :: O_ddp = 2
  integer, parameter :: O_ddd = 3
  integer, parameter :: O_pds = 4
  integer, parameter :: O_pdp = 5
  integer, parameter :: O_pps = 6
  integer, parameter :: O_ppp = 7
  integer, parameter :: O_sds = 8
  integer, parameter :: O_sps = 9
  integer, parameter :: O_sss = 10

#define HTAB 1:MAX_NORB
#define STAB MAX_NORB+1:2*MAX_NORB

  character(3), parameter :: electronic_configuration(9) = &
     ["s  ", "---", " p ", "sp ", "  d", "s d", "---", " pd", "spd"]

  ! IF YOU MODIFY THIS STRUCTURE, *ALWAYS* ALSO MODIFY THE CORRESPONDING
  ! STRUCTURE IN materials.h
  public :: notb_element_t
  type, bind(C) :: notb_element_t

     logical(C_BOOL)         :: exists    = .false.

     character(kind=C_CHAR)  :: name(2)   = ["X","X"]  ! name of element
     character(kind=C_CHAR)  :: cname(10) = ["n","o","n","a","m","e"," "," "," "," "]  ! common name of element
     integer(C_INT)          :: elem      = 10000      ! number of element (official)
     integer(C_INT)          :: no        = 10000      ! number of orbitals
     integer(C_INT)          :: l(9)      = [0,1,1,1,2,2,2,2,2] !angular momenta of orbitals
     integer(C_INT)          :: lmax      = 1000       ! maximum angular momentum
     real(C_DOUBLE)          :: e(9)      = 1E30       ! on-site energies [ e(1:no) ]
     real(C_DOUBLE)          :: el_max    = 0          ! max number of valence electrons on an atom
     real(C_DOUBLE)          :: U         = 1E30       ! Hubbard U
     real(C_DOUBLE)          :: q0        = 1E30       ! charge (nr of electrons in neutral)
  
     ! variables for HOTBIT
     integer(C_INT)          :: o1        = 1E5        ! index of the first orbital
     integer(C_INT)          :: enr       = 1E5        ! element number in the internal book-keeping

     ! spin-related variables
     logical(C_BOOL)         :: spin      = .false.    ! spin-parameters set?
     real(C_DOUBLE)          :: W(0:2,0:2)             ! W parameter values, 0,1,2 = s,p,d, W(0,0) = Wss, W(0,1) = Wsp etc.

  endtype notb_element_t

  public :: materials_t
  type materials_t

     character(1000)               :: folder

     integer                       :: nel             ! number of elements in materials database

     type(notb_element_t), pointer :: e(:)            ! elements in the material database

     real(DP), pointer             :: cut(:, :)       ! cut-off for the Slater-Koster tables

     type(spline_t), pointer       :: HS(:, :)        ! the Hamiltonian and overlap matrix
     type(spline_t), pointer       :: R(:, :)         ! repulsive potential

  endtype materials_t


  public :: read_database
  interface read_database
     module procedure materials_read_database
  endinterface

  public :: write_tables
  interface write_tables
     module procedure materials_write_tables
  endinterface

  public :: element_by_symbol
  interface element_by_symbol
     module procedure materials_element_by_symbol
  endinterface

  public :: element_by_Z
  interface element_by_Z
     module procedure materials_element_by_Z
  endinterface

  public :: get_orbital
  interface get_orbital
     module procedure materials_get_orbital
  endinterface

  integer, parameter :: valence_orbitals(116) =  &
     [ 1,-1, &
       4, 4, 4, 4, 4, 4, 4,-1, &
       4, 4, 4, 4, 4, 4, 4,-1, &
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
      -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,     &
      -1,-1,-1,-1]

  !>
  !! Temporary storage for data read from bonds.bx file
  !<
  type bopfox_table_t
     integer               :: n = 0
     real(DP), allocatable :: x(:)
     real(DP), allocatable :: HS(:, :)
  endtype bopfox_table_t

contains

  !>
  !! Convert condensed orbital index into absolute orbital index
  !<
  function materials_get_orbital(no, a0) result(a)
    implicit none

    integer, intent(in) :: no, a0
    integer             :: a

    ! ---
    
    a = a0
    ! if no == 1, no == 4 or no == 9 we are all set
    !     ^s       ^sp        ^spd
    ! if no == 5, this element has just d orbitals defined
    if (no == 5)  a = a+4
    ! if no == 5, this element has p and d orbitals defined
    if (no == 8)  a = a+1
    ! if no == 6, this element has just s and d orbitals defined
    if (no == 6 .and. a > 1)  a = a+3

  endfunction materials_get_orbital


  !>
  !! Convert string to all lower case
  !<
  elemental subroutine lowercase(s)
    implicit none

    character(*), intent(inout)  :: s

    ! ---

    integer  :: i

    ! ---

    do i = 1, len(s)
       if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') then
          s(i:i) = char(ichar(s(i:i))+32)
       endif
    enddo

  endsubroutine lowercase


  !>
  !! reads info-lines from the beginning of a datafile
  !! that is supposed to read. The info (or comment) lines
  !! start with the letter '#', compatible with GNUPLOT and many
  !! other comment lines.
  !<
  subroutine filestart(un)
    implicit none

    integer, intent(in) :: un

    ! ---

    integer         :: i,k,io
    character(1000) :: line

    ! ---

    i=0
    do
       read (un, '(1000a)',iostat=io) line
       if (io/=0) exit
       if (line(1:1)/="#") exit
       i = i+1
       if (i>1000) stop "filestart: too many comment lines?"
    end do
    rewind(un)
    do k = 1, i
       read (un, '(1000a)') line
    enddo

  endsubroutine filestart


  !>
  !! Clean the string from not-nice ascii characters (like carriage returns)
  !<
  subroutine clean_string(str)
    implicit none

    character(*), intent(inout) :: str

    ! ---

    integer :: i,j,asc

    ! ---

    do i = 1, len(str)
       asc = ichar(str(i:i))
       if (asc<32 .or. asc==127) exit
    enddo
    do j = i, len(str)
       str(i:i)=' '
    enddo
  endsubroutine clean_string


  !>
  !! from opened file unit 'un', find key-value pairs in the format
  !! key = value
  !!
  !! Example:    (code)                           (data.in)
  !!           open(10,file='data.in')            ...(some data)...
  !!           call find_key(10,'mass',re=m)      mass = 1.234
  !!           close(10)                          ...(some data)...
  !>
  subroutine find_value(un, key, re, in, ch, lg, str, test, ignore, error)
    implicit none

    integer,                intent(in)    :: un
    character(*),           intent(in)    :: key
    real(8),      optional, intent(inout) :: re
    integer,      optional, intent(inout) :: in
    logical,      optional, intent(inout) :: lg
    character(*), optional, intent(inout) :: ch
    character(*), optional, intent(inout) :: str
    logical,      optional, intent(out)   :: test
    logical,      optional, intent(in)    :: ignore 
    integer,      optional, intent(out)   :: error

    ! ---

    character(500) :: line
    integer :: io,i

    ! ---

    INIT_ERROR(error)

    rewind(un)
    do
       ! try to find the key, if not found, exit
       read (un,'(500a)',iostat=io) line
       call clean_string(line)
       if (io/=0) then
          inquire(unit=un,name=line)
          ! call flog('Key '//trim(keyd)//' not found in file '//trim(line))

          if (present(test)) then
             test=.false.
             return
          else if(present(ignore)) then
             if(ignore) then
                return
             else
                RAISE_ERROR("End of file "//trim(line)//", key"//trim(key)//"not found.", error)
             endif
          else
             RAISE_ERROR("End of file, key not found.", error)
          endif
       endif
       line = adjustl(line)
       if (line(1:1)=='#') cycle !don't consider comments
       i = scan(line,'=')
       if (i==0) cycle
       if (trim(line(1:i-1))/=trim(adjustl(key))) cycle

       ! -----------------------------
       ! key was found, read the value
       ! -----------------------------
       if (present(re)) then
          read(line(i+1:),*) re
       else if (present(in)) then
          read(line(i+1:),*) in
       else if (present(ch)) then
          read(line(i+1:),*) ch
          ch=trim(ch)
       else if (present(test)) then
          test=.true.
       else if (present(str)) then
          str=trim(line(i+1:))
       else if (present(lg) ) then
          if (trim(line(i+1:))=='T' .or. trim(line(i+1:))=='TRUE' .or. &
               trim(line(i+1:))=='yes' .or. trim(line(i+1:))=='y' ) then
             lg=.true.
          else if (trim(line(i+1:))=='F' .or. trim(line(i+1:))=='FALSE' .or. &
               trim(line(i+1:))=='no' .or. trim(line(i+1:))=='n') then
             lg=.false.
          else
             RAISE_ERROR("Not a valid logical value:"//line, error)
          endif
       endif

       !--------------------------------
       ! if re,in,ch not present,
       ! leave the cursor in this place
       ! for reading of more complicated
       ! data structure
       !--------------------------------
       return
    enddo

  endsubroutine find_value


  !>
  !! Returns the internal element number for given symbol
  !<
  logical function materials_element_by_symbol(this, sym, el, enr) result(r)
    implicit none
    
    type(materials_t), intent(in)                :: this
    character(2), intent(in)                     :: sym
    type(notb_element_t), intent(out), optional  :: el
    integer, intent(out), optional               :: enr
    

    ! ---

    integer  :: i

    ! ---

    r  = .false.
    do i = 1, this%nel
       if (this%e(i)%exists) then
          if (trim(a2s(this%e(i)%name)) == trim(sym)) then
             if (present(el)) then
                el  = this%e(i)
             endif
             if (present(enr)) then
                enr  = i
             endif
             r   = .true.
          endif
       endif
    enddo

  endfunction materials_element_by_symbol


  !>
  !! Returns the internal element number for given symbol
  !<
  logical function materials_element_by_Z(this, Z, el, enr) result(r)
    implicit none
    
    type(materials_t), intent(in)                :: this
    integer, intent(in)                          :: Z
    type(notb_element_t), intent(out), optional  :: el
    integer, intent(out), optional               :: enr

    ! ---

    integer  :: i

    ! ---

    r  = .false.
    do i = 1, this%nel
       if (this%e(i)%exists) then
          if (this%e(i)%elem == Z) then
             if (present(el)) then
                el  = this%e(i)
             endif
             if (present(enr)) then
                enr  = i
             endif
             r   = .true.
          endif
       endif
    enddo

  endfunction materials_element_by_Z


  !>
  !! Load the Slater-Koster tables (HOTBIT format)
  !<
  subroutine materials_read_sltab_hotbit(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout)  :: db
    real(DP), intent(in)              :: econv, lconv
    integer, intent(inout), optional  :: error

    ! ---

    integer          :: un, i1, i2
    character(2)     :: e1, e2
    character(1000)  :: fil, fil2
    logical          :: ex, ex2, vex

    real(DP)         :: conv(2*MAX_NORB)

    ! ---

    call prlog("- materials_read_sltab_hotbit -")

    conv(HTAB) = econv
    conv(STAB) = 1.0

    do i1=1,db%nel
       if (db%e(i1)%exists) then
          do i2=i1,db%nel
             if (db%e(i2)%exists) then
                e1=trim(a2s(db%e(i1)%name))
                e2=trim(a2s(db%e(i2)%name))

                fil =trim(db%folder)//'/'//trim(e1)//'_'//trim(e2)//'.par'
                fil2=trim(db%folder)//'/'//trim(e2)//'_'//trim(e1)//'.par'
                inquire(file=fil ,exist=ex)
                inquire(file=fil2,exist=ex2)
                if( ex ) then
                   un = fopen(fil)
                   call prlog("HOTBIT tables for "//trim(e1)//"-"//trim(e2)//" found.")
                else if( ex2 ) then
                   un = fopen(fil2)
                   call prlog("HOTBIT tables for "//trim(e1)//"-"//trim(e2)//" found.")
                else
!                   write (ilog, '(5X,A,A,A,A,A)')  "WARNING: Skipping parametrizations for "//e1//" and "//e2//". Could not find '", trim(fil), "' or '", trim(fil2), "' file."
                   if( i1/=i2 )then
                      db%cut(i1,i2) = -1
                      db%cut(i2,i1) = -1
                   else
                      db%cut(i1,i2) = 10
                   end if
                   cycle
                end if

                ! 
                ! Read repulsive potentials.
                !
                call find_value(un,'repulsion',test=vex,error=error) 
                PASS_ERROR(error)
                call read2(db%R(i1, i2), un, 2, lconv, (/ econv /), error)
                PASS_ERROR(error)
                if (i1 /= i2) then
                   call associate(db%R(i2, i1), db%R(i1, i2))
                endif

#ifdef DEBUG
                call write(db%R(i1, i2), "rep_"//trim(e1)//'_'//trim(e2)//".out")
#endif

                !
                ! read matrix elements for pairs of atom species, first
                ! <1|...|2> and then <2|...|1>
                ! ordering: dds ddp ddd pds pdp pps ppp sds sps sss 
                ! (first for H, then for S)
                !
                call find_value(un,trim(e1)//'_'//trim(e2)//'_table')
!                write (*, '(A,X,A,X,A,X,I5,I5)')  "table 1", e1, e2, i1, i2
                call read2(db%HS(i1, i2), un, 2*MAX_NORB+1, lconv, conv, error)
                PASS_ERROR(error)

#ifdef DEBUG
                call write(db%HS(i1, i2), "HS_"//trim(e1)//'_'//trim(e2)//".out")
#endif

                if( i1/=i2 ) then
                   call find_value(un,trim(e2)//'_'//trim(e1)//'_table')
!                   write (*, *)  "table 2", e2, e1
                   call read2(db%HS(i2, i1), un, 2*MAX_NORB+1, lconv, conv, error)
                   PASS_ERROR(error)

#ifdef DEBUG
                   call write(db%HS(i2, i1), "HS_"//trim(e2)//'_'//trim(e1)//".out")
#endif
                end if

                call fclose(un)

             endif
          enddo
       endif
    enddo

    call prlog

  endsubroutine materials_read_sltab_hotbit


  !>
  !! Load the Slater-Koster tables (DFTB format)
  !<
  subroutine materials_read_sltab_dftb(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout)  :: db
    real(DP), intent(in)              :: econv, lconv
    integer, intent(inout), optional  :: error

    ! ---

    real(DP), parameter  :: REP_DX    = 0.005_DP
    real(DP), parameter  :: REP_X0    = 0.0_DP

    integer, parameter   :: MAX_DATA  = 10000

    ! ---

    integer                :: un, i, j, i1, i2, n, io
    character(2)           :: e1, e2, f1, f2
    character(1000)        :: fil, fil3
    logical                :: ex, ex3

    real(DP)               :: conv(2*MAX_NORB)

    real(DP)               :: eself(3), espin, u(3)
    real(DP)               :: cx, dx, c1, c2, c3, x1, x2, splc(6), cutoff
    real(DP)               :: q(3)
    character(200)         :: line

    real(DP), allocatable  :: x(:), y(:)

    ! ---

    call prlog("- materials_read_sltab_dftb -")

    allocate(x(MAX_DATA), y(MAX_DATA))

    conv(HTAB) = econv
    conv(STAB) = 1.0

    do i1 = 1, db%nel
       if (db%e(i1)%no > 0) then
          do i2 = 1, db%nel
             if (db%e(i2)%no > 0) then
                e1    = trim(a2s(db%e(i1)%name))
                e2    = trim(a2s(db%e(i2)%name))

                f1    = e1
                call lowercase(f1)
                f1 = uppercase(f1(1:1)) // f1(2:2)
                f2    = e2
                call lowercase(f2)
                f2 = uppercase(f2(1:1)) // f2(2:2)

                fil   = trim(db%folder)//'/'//trim(f1)//trim(f2)//'.spl'
                !fil3  = trim(db%folder)//'/'//trim(uppercase(f1))//'-'//trim(uppercase(f2))//'.skf'
                fil3  = trim(db%folder)//'/'//trim(f1)//'-'//trim(f2)//'.skf'
                inquire(file=fil, exist=ex)
                inquire(file=fil3, exist=ex3)

                un = -1

                if (ex) then
                   un = fopen(fil)
                   call prlog("DFTB tables for "//trim(e1)//"-"//trim(e2)//" found.")
                else if (ex3) then
                   un = fopen(fil3)
                   call prlog("DFTB tables for "//trim(e1)//"-"//trim(e2)//" found.")
                endif

                file_exists: if (un > 0) then

                   if (db%HS(i1, i2)%n > 0 .or. db%R(i1, i2)%n > 0) then
                      call prlog("WARNING: "//e1//"-"//e2//" tables already read.")
                   endif

                   ! cutoff, number of grid points in Hamiltonian/Overlap
                   read (un, *) dx, n
                   n = n-1 ! Sometimes, there seems to be a line thats missing

                   do i = 1, n
                      x(i) = (i-1)*dx
                   enddo

                   if (i1 == i2) then
                      ! We were able to get fundamental data on this element
                      db%e(i1)%exists = .true.

                      ! self energies, spin polarization energy(?), Hubbard U's, number of electrons
                      read (un, *) eself(:), espin, u(:), q(:)

                      eself  = eself * econv

                      if (i1 == i2) then
                         call prlog("WARNING: Overriding self-energies and Hubbard-U from 'elements.dat'.")

                         db%e(i1)%e  = (/ &
                              eself(3), &
                              eself(2), eself(2), eself(2), &
                              eself(1), eself(1), eself(1), eself(1), eself(1) &
                              /)

                         db%e(i1)%U  = u(3) * econv
                      endif

                      db%e(i1)%exists = .true.

                   endif

                   ! read spline data
                   call read(db%HS(i1, i2), un, 2*MAX_NORB+1, lconv, conv, n, x, error)
                   PASS_ERROR(error)

#ifdef DEBUG
                   call write(db%HS(i1, i2), "HS_"//trim(e1)//'_'//trim(e2)//".out")
#endif

                   ! Look for 'Spline'
                   read (un, *, iostat=io) line
                   do while (trim(line) /= "Spline" .and. io == 0)
                      read (un, *, iostat=io) line
                   enddo

                   if (io /= 0) then
                      RAISE_ERROR("End of file reached while looking for keyword 'Spline'.", error)
                   endif

                   !
                   ! Reading the repulsive part. We will construct a tabulated version of the repulsion
                   ! and then use our spline module.
                   !

                   read (un, *) n, cutoff
                   read (un, *) c1, c2, c3

                   read (un, *) x1, x2, splc(1:4)

                   !
                   ! The tail of the repulsive function is given by an exponential
                   !

                   cx  = REP_X0
                   i   = 1
                   do while (cx < x1)
                      x(i)  = cx
                      y(i)  = c3 + exp(c2-c1*cx)

                      cx  = cx + REP_DX
                      i   = i + 1
                   enddo

                   n = n - 1

                   !
                   ! The rest is spline coefficients
                   !

                   do while (n > 1)

                      do while (cx < x2)
                         x(i)  = cx
                         y(i)  = splc(1)
                         dx    = cx-x1
                         do j = 2, 4
                            y(i)  = y(i) + splc(j)*dx
                            dx    = dx*(cx-x1)
                         enddo

                         cx  = cx + REP_DX
                         i  = i + 1
                      enddo

                      read (un, *) x1, x2, splc(1:4)

                      n  = n - 1

                   enddo

                   do while (cx < x2)
                      x(i)  = cx
                      y(i)  = splc(1)
                      dx    = cx-x1
                      do j = 2, 4
                         y(i)  = y(i) + splc(j)*dx
                         dx    = dx*(cx-x1)
                      enddo

                      cx  = cx + REP_DX
                      i  = i + 1
                   enddo

                   !
                   ! The last one is an fifth order polynomial
                   !

                   read (un, *) x1, x2, splc(1:6)

                   do while (cx < x2)
                      x(i)  = cx
                      y(i)  = splc(1)
                      dx    = cx-x1
                      do j = 2, 6
                         y(i)  = y(i) + splc(j)*dx
                         dx    = dx*(cx-x1)
                      enddo

                      cx  = cx + REP_DX
                      i   = i + 1
                   enddo

                   if (i > MAX_DATA) then
                      RAISE_ERROR("i > MAX_DATA", error)
                   endif

                   !
                   ! Construct spline
                   !

                   x  = x * lconv
                   y  = y * econv
                   call nonuniform_spline_init(db%R(i1, i2), MAX_DATA, i-1, x, 1, (/ y /))

#ifdef DEBUG
                   call write(db%R(i1, i2), "rep_"//trim(e1)//'_'//trim(e2)//".out")
#endif

                   call fclose(un)

                endif file_exists

             endif
          enddo
       endif
    enddo

    deallocate(x)
    deallocate(y)

    call prlog

  endsubroutine materials_read_sltab_dftb


  !>
  !! Load the Slater-Koster tables (DFTB format)
  !<
  subroutine materials_read_sltab_bopfox(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout) :: db
    real(DP),          intent(in)    :: econv, lconv
    integer, optional, intent(inout) :: error

    ! ---

    integer               :: i, k, n, eli, elj, noi, noj, un, io
    real(DP)              :: scaling
    real(DP), allocatable :: d(:, :)
    character(1024)       :: fn, line, values, key
    character(2)          :: symi, symj
    character(3)          :: vali, valj
    logical               :: ex
 
    type(bopfox_table_t)  :: tab(db%nel, db%nel)

    ! ---

    call prlog("- materials_read_sltab_bopfox -")

    fn = trim(db%folder) // "/bonds.bx"
    inquire(file=fn, exist=ex)

    if (ex) then
       un = fopen(trim(fn))
   
       eli = -1
       elj = -1
       do
          read (un, '(200a)', iostat=io)  line
          if (io /= 0)  exit ! EOF
          if (line(1:2) == '/')  cycle ! Comment
          k = scan(line, '=')
          if (k /= 0) then
             key    = lower_case(adjustl(line(1:k-1)))
             values = line(k+1:)
          else
             ! Skip all lines that are different from "key = values"
             cycle
          endif
      
          select case(trim(key))
          case("bond")   ! starts the set for new element
             read (values, *)  symi, symj
             if (.not. element_by_symbol(db, symi, enr=eli)) then
                RAISE_ERROR("Could no find element '"//trim(symi)//"' in NOTB database.", error)
             endif
             if (.not. element_by_symbol(db, symj, enr=elj)) then
                RAISE_ERROR("Could no find element '"//trim(symj)//"' in NOTB database.", error)
             endif
          case("valence")
             read (values, *)  vali, valj
             noi = 0
             noj = 0
             do i = 1, 3
                if (vali(i:i) == 's')  noi = noi + 1
                if (vali(i:i) == 'p')  noi = noi + 3
                if (vali(i:i) == 'd')  noi = noi + 5
                if (valj(i:i) == 's')  noj = noj + 1
                if (valj(i:i) == 'p')  noj = noj + 3
                if (valj(i:i) == 'd')  noj = noj + 5
             enddo
             if (noi /= db%e(eli)%no) then
                RAISE_ERROR("'bonds.bx' reports '"//trim(vali)//"' valence for element '"//trim(symi)//"' ("//noi//" orbitals), but 'atoms.bx' reports "//db%e(eli)%no//" orbitals.", error)
             endif
             if (noj /= db%e(elj)%no) then
                RAISE_ERROR("'bonds.bx' reports '"//trim(valj)//"' valence for element '"//trim(symj)//"' ("//noj//" orbitals), but 'atoms.bx' reports "//db%e(elj)%no//" orbitals.", error)
             endif
          case("scaling")
             read (values, *)  scaling
             if (abs(scaling - 1.0) > 1e-9) then
                RAISE_ERROR("'bonds.bx' reports scaling != 1 for "//trim(symi)//"-"//trim(symj)//" bond integrals. Don't know how to handle this.", error)
             endif
          case("bondtable")
             call prlog("Hamiltonian for "//trim(symi)//"-"//trim(symj)//" found.")

             read (values, *)  n
             allocate(d(15, n))
             read (un, *)  d

             if (tab(eli, elj)%n > 0) then
                if (tab(eli, elj)%n /= n) then
                   RAISE_ERROR("Mismatch in number of grid points for Hamiltonian and overlap tables for "//trim(symi)//"-"//trim(symj)//" bond integrals.", error)
                endif
                if (any(abs(tab(eli, elj)%x - d(1, :)) > 1d-9)) then
                   RAISE_ERROR("Mismatch in number of grid positions for Hamiltonian and overlap tables for "//trim(symi)//"-"//trim(symj)//" bond integrals.", error)
                endif
             else
                tab(eli, elj)%n = n
                allocate(tab(eli, elj)%x(n))
                allocate(tab(eli, elj)%HS(n, 20))
                tab(eli, elj)%x = d(1, :)

                if (eli /= elj) then
                   tab(elj, eli)%n = n
                   allocate(tab(elj, eli)%x(n))
                   allocate(tab(elj, eli)%HS(n, 20))
                   tab(elj, eli)%x = d(1, :)
                endif
             endif

             ! BOPFOX is a bit more clever. It stores the 14 independent
             ! bond-integrals for a pair of elements. We store them separately
             ! for pairs i-j and j-i, which makes ten bond-integrals each.
             ! We need to spread BOPFOXs data out.

             tab(eli, elj)%HS(:, O_sss) = d(2, :) ! sss
             tab(elj, eli)%HS(:, O_sss) = d(2, :) ! sss
             tab(eli, elj)%HS(:, O_sps) = d(3, :) ! sps
             tab(elj, eli)%HS(:, O_sps) = d(4, :) ! pss
             tab(eli, elj)%HS(:, O_pps) = d(5, :) ! pps
             tab(elj, eli)%HS(:, O_pps) = d(5, :) ! pps
             tab(eli, elj)%HS(:, O_ppp) = d(6, :) ! ppp
             tab(elj, eli)%HS(:, O_ppp) = d(6, :) ! ppp
             tab(eli, elj)%HS(:, O_sds) = d(7, :) ! sds
             tab(elj, eli)%HS(:, O_sds) = d(8, :) ! dss
             tab(eli, elj)%HS(:, O_pds) = d(9, :) ! pds
             tab(elj, eli)%HS(:, O_pds) = d(10, :) ! dps
             tab(eli, elj)%HS(:, O_pdp) = d(11, :) ! pdp
             tab(elj, eli)%HS(:, O_pdp) = d(12, :) ! dpp
             tab(eli, elj)%HS(:, O_dds) = d(13, :) ! dda
             tab(elj, eli)%HS(:, O_dds) = d(13, :) ! dda
             tab(eli, elj)%HS(:, O_ddp) = d(14, :) ! ddp
             tab(elj, eli)%HS(:, O_ddp) = d(14, :) ! ddp
             tab(eli, elj)%HS(:, O_ddd) = d(15, :) ! ddd
             tab(elj, eli)%HS(:, O_ddd) = d(15, :) ! ddd

             deallocate(d)
          case("overtable")
             call prlog("Overlap matrix for "//trim(symi)//"-"//trim(symj)//" found.")

             read (values, *)  n
             allocate(d(15, n))
             read (un, *)  d

             if (tab(eli, elj)%n > 0) then
                if (tab(eli, elj)%n /= n) then
                   RAISE_ERROR("Mismatch in number of grid points for Hamiltonian and overlap tables for "//trim(symi)//"-"//trim(symj)//" bond integrals.", error)
                endif
                if (any(abs(tab(eli, elj)%x - d(1, :)) > 1d-9)) then
                   RAISE_ERROR("Mismatch in number of grid positions for Hamiltonian and overlap tables for "//trim(symi)//"-"//trim(symj)//" bond integrals.", error)
                endif
             else
                tab(eli, elj)%n = n
                allocate(tab(eli, elj)%x(n))
                allocate(tab(eli, elj)%HS(n, 20))
                tab(eli, elj)%x = d(1, :)

                if (eli /= elj) then
                   tab(elj, eli)%n = n
                   allocate(tab(elj, eli)%x(n))
                   allocate(tab(elj, eli)%HS(n, 20))
                   tab(elj, eli)%x = d(1, :)
                endif
             endif

             ! BOPFOX is a bit more clever. It stores the 14 independent
             ! bond-integrals for a pair of elements. We store them separately
             ! for pairs i-j and j-i, which makes ten bond-integrals each.
             ! We need to spread BOPFOXs data out.

             tab(eli, elj)%HS(:, 10+O_sss) = d(2, :) ! sss
             tab(elj, eli)%HS(:, 10+O_sss) = d(2, :) ! sss
             tab(eli, elj)%HS(:, 10+O_sps) = d(3, :) ! sps
             tab(elj, eli)%HS(:, 10+O_sps) = d(4, :) ! pss
             tab(eli, elj)%HS(:, 10+O_pps) = d(5, :) ! pps
             tab(elj, eli)%HS(:, 10+O_pps) = d(5, :) ! pps
             tab(eli, elj)%HS(:, 10+O_ppp) = d(6, :) ! ppp
             tab(elj, eli)%HS(:, 10+O_ppp) = d(6, :) ! ppp
             tab(eli, elj)%HS(:, 10+O_sds) = d(7, :) ! sds
             tab(elj, eli)%HS(:, 10+O_sds) = d(8, :) ! dss
             tab(eli, elj)%HS(:, 10+O_pds) = d(9, :) ! pds
             tab(elj, eli)%HS(:, 10+O_pds) = d(10, :) ! dps
             tab(eli, elj)%HS(:, 10+O_pdp) = d(11, :) ! pdp
             tab(elj, eli)%HS(:, 10+O_pdp) = d(12, :) ! dpp
             tab(eli, elj)%HS(:, 10+O_dds) = d(13, :) ! dda
             tab(elj, eli)%HS(:, 10+O_dds) = d(13, :) ! dda
             tab(eli, elj)%HS(:, 10+O_ddp) = d(14, :) ! ddp
             tab(elj, eli)%HS(:, 10+O_ddp) = d(14, :) ! ddp
             tab(eli, elj)%HS(:, 10+O_ddd) = d(15, :) ! ddd
             tab(elj, eli)%HS(:, 10+O_ddd) = d(15, :) ! ddd

             deallocate(d)
          case("reptable")
             call prlog("Repulsion for "//trim(symi)//"-"//trim(symj)//" found.")

             read (values, *)  n
             allocate(d(2, n))
             read (un, *)  d

             call init(db%R(eli, elj), n, n, d(1, :), 1, d(2:2, :))
             call init(db%R(elj, eli), n, n, d(1, :), 1, d(2:2, :))

             deallocate(d)
          case default
             cycle
          end select
       end do

       call fclose(un)

       ! Initialize splines
       do eli = 1, db%nel
          do elj = 1, db%nel
             if (tab(eli, elj)%n > 0) then
                call init(db%HS(eli, elj), tab(eli, elj)%n, tab(eli, elj)%n, &
                          tab(eli, elj)%x, 20, tab(eli, elj)%HS)

                deallocate(tab(eli, elj)%x)
                deallocate(tab(eli, elj)%HS)
             endif
          enddo
       enddo
    else
       call prlog("No 'bonds.bx' file found in database directory. Skipping.")
    endif

    call prlog

  endsubroutine materials_read_sltab_bopfox


  !>
  !! Reads element data from HOTBITs elements.dat
  !<
  subroutine materials_read_elements_hotbit(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout) :: db
    real(DP),          intent(in)    :: econv, lconv
    integer, optional, intent(out)   :: error

    ! ---

    integer               :: i, j, k, un, io, lmx(9)
    character(200)        :: line, dat, key, fn
    logical               :: ex

    type(notb_element_t)  :: hlp(MAX_Z)

    ! ---

    INIT_ERROR(error)

    call prlog("- materials_read_elements_hotbit -")

    lmx = 1000
    lmx(1)=0; lmx(4)=1; lmx(9)=2 !lmax = lmx(no)

    fn = trim(db%folder) // "/elements.dat"

    inquire(file=fn, exist=ex)

    if (ex) then
       un = fopen(trim(fn))
       call filestart(un)
   
       j=0
       do
          read(un,'(200a)',iostat=io) line  
          if( io/=0 ) exit !EOF
          k = scan(line,'=')
          if( k/=0 ) then
             key = adjustl(line(1:k-1))
             dat = line(k+1:)
          else
             cycle 
          end if
      
          select case(trim(key))
          case("element")   ! starts the set for new element
             j=j+1
             hlp(j)%name = ' '
             hlp(j)%name(1:min(2,len_trim(dat))) = s2a(trim(dat))
          case("Z");     read(dat,*) hlp(j)%elem
          case("common");
             hlp(j)%cname = ' '
             hlp(j)%cname(1:min(10,len_trim(dat))) = s2a(trim(dat))
          case("q0");    read(dat,*) hlp(j)%q0
          case("no")
             read(dat,*) hlp(j)%no
             hlp(j)%lmax = lmx(hlp(j)%no)
             read(un ,*) hlp(j)%e(1:hlp(j)%no)
             hlp(j)%e(1:hlp(j)%no) = hlp(j)%e(1:hlp(j)%no) * econv
          case("U")
             read(dat,*) hlp(j)%U
             hlp(j)%U = hlp(j)%U * econv
          case default
             cycle
          end select
       end do
   
       call fclose(un)
   
       ! set dependent variables and sort according to Z
       do i=1,j
   
          if (hlp(i)%elem > 0 .and. hlp(i)%elem <= MAX_Z) then
             if (trim(a2s(hlp(i)%name)) /= trim(ElementName(hlp(i)%elem))) then
                call prlog("WARNING: Name '"//a2s(hlp(i)%name)//"' in 'elements.dat' not equal common element name '"//ElementName(hlp(i)%elem)//"'.")
          endif

          hlp(i)%el_max = 0d0
          do k=1,hlp(j)%no
             hlp(i)%el_max = hlp(i)%el_max + (2d0*hlp(i)%l(k)+1d0)
          end do
          db%e(hlp(i)%elem) = hlp(i)
          db%e(hlp(i)%elem)%exists = .true.
          db%e(hlp(i)%elem)%enr = hlp(i)%elem
       else
          call prlog("WARNING: Unknown element found in 'elements.dat' (Z = "//hlp(i)%elem//").")
       endif
   
       end do
   
       call prlog(""//j//" elements found in 'elements.dat'.")
    else
       call prlog("No 'elements.dat' found in database directory. Skipping.")
    endif
    call prlog

  endsubroutine materials_read_elements_hotbit


  !>
  !! Reads element data from HOTBITs elements.dat
  !<
  subroutine materials_read_elements_bopfox(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout) :: db
    real(DP),          intent(in)    :: econv, lconv
    integer, optional, intent(out)   :: error

    ! ---

    ! Translation table for population of on-site energies read from BOPFOX
    ! file. Number is index of energy in BOPFOX file.
    integer, parameter :: l(9, 9) = &
       reshape([ 1, 0, 0, 0, 0, 0, 0, 0, 0,   & ! s
                 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                 0, 1, 1, 1, 1, 0, 0, 0, 0,   & ! p
                 1, 2, 2, 2, 0, 0, 0, 0, 0,   & ! sp
                 0, 0, 0, 0, 1, 1, 1, 1, 1,   & ! d
                 1, 0, 0, 0, 2, 2, 2, 2, 2,   & ! sd
                 0, 0, 0, 0, 0, 0, 0, 0, 0,   &
                 0, 1, 1, 1, 2, 2, 2, 2, 2,   & ! pd
                 1, 2, 2, 2, 3, 3, 3, 3, 3 ], & ! spd
                 [9, 9])

    integer              :: i, j, k, un, io
    real(DP)             :: onsitelevels(3)
    character(1024)      :: fn, line, values, key
    logical              :: ex

    type(notb_element_t) :: e(MAX_Z)

    ! ---

    INIT_ERROR(error)

    call prlog("- materials_read_elements_bopfox -")

    fn = trim(db%folder) // "/atoms.bx"

    inquire(file=fn, exist=ex)

    if (ex) then
       un = fopen(trim(fn))
       call filestart(un)
   
       j = 0
       do
          read(un, '(1024a)', iostat=io)  line
          if (io /= 0) exit !EOF
          if (line(1:2) == '/')  cycle ! Comment
          k = scan(line, '=')
          if (k /= 0) then
             key = lower_case(adjustl(line(1:k-1)))
             values = adjustl(line(k+1:))
          else
             ! Skip all lines that are different from "key = value"
             cycle
          endif

          select case(trim(key))
          case("atom")   ! starts the set for new element
             j = j + 1
             e(j)%name = '  '
             e(j)%name(1:min(2, len_trim(values))) = s2a(trim(values))
             e(j)%elem = atomic_number_from_symbol(a2s(e(j)%name))
          case("valenceorbitals")
             read(values, *)  e(j)%no
          case("valenceelectrons")
             read(values,*)  e(j)%q0
          case("onsitelevels")
             k = maxval(l(:, e(j)%no))
             read(values, *)  onsitelevels(1:k)
             do i = 1, 9
                if (l(i, e(j)%no) > 0) then
                   e(j)%e(i) = onsitelevels(l(i, e(j)%no))
                endif
             enddo
          case("jii")
             read(values, *)  e(j)%U
          case default
             cycle
          end select
       end do
   
       call fclose(un)
   
       ! set dependent variables and sort according to Z
       do i = 1, j  
          e(i)%el_max = 0d0
          do k = 1, e(j)%no
             e(i)%el_max = e(i)%el_max + (2d0*e(i)%l(k)+1d0)
          enddo
          db%e(e(i)%elem) = e(i)
          db%e(e(i)%elem)%exists = .true.
          db%e(e(i)%elem)%enr = e(i)%elem
       enddo
   
       call prlog(""//j//" elements found in 'atoms.bx'.")
    else
       call prlog("No 'atoms.bx' found in database directory. Skipping.")
    endif
    call prlog

  endsubroutine materials_read_elements_bopfox


  !>
  !! Reads element data from elements.dat
  !<
  subroutine materials_init_elements(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout) :: db
    real(DP),          intent(in)    :: econv, lconv
    integer, optional, intent(out)   :: error

    ! ---

    integer :: i

    ! ---

    INIT_ERROR(error)

    db%nel = MAX_Z
    allocate(db%e(db%nel))
    ! Initialize default valency. This can be overridden by the Slater-Koster
    ! tables.
    do i = 1, MAX_Z
       db%e(i)%name = s2a(ElementName(i))
       db%e(i)%elem = i
       db%e(i)%enr  = i
       db%e(i)%q0   = valence_orbitals(i)
       db%e(i)%no   = valence_orbitals(i)
    enddo

    call materials_read_elements_hotbit(db, econv, lconv, error)
    PASS_ERROR(error)
    call materials_read_elements_bopfox(db, econv, lconv, error)
    PASS_ERROR(error)

  endsubroutine materials_init_elements


  !>
  !! Read spin-dependent parameters (W)
  !!
  !! Read spin-dependent parameters (W)
  !!
  !! The expected format is one comment line, then one line per
  !! element, each with element name and W parameters:
  !!   el Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd
  !<
  subroutine materials_read_spin_params(db, econv, error)
    implicit none

    type(materials_t), intent(inout)  :: db                   !< Materials database
    real(DP), intent(in)              :: econv                !< Energy unit conversion
    integer, intent(inout), optional  :: error               !< Errors

    ! ---

    integer                :: i                               ! loops
    character(2)           :: el                              ! current element name

    character(1000)        :: file = "spin_parameters.dat"    ! input file name
    logical                :: exists                          ! the file exists?
    integer                :: un                              ! unit for the file
    integer                :: io                              ! for checking reading status
    character(1000)        :: dummy                           ! for reading comment line

    real(DP)               :: W(9)                            ! for reading in values from file
    logical                :: inserted                        ! currently read parameters to database succesfully?

    ! ---

    call prlog("- materials_read_spin_params -")

    ! reset
    do i = 1, db%nel
       db%e(i)%W = 0.0_DP
    end do

    ! open file and return if none found
    file = trim(db%folder)//'/'//trim(file)
    inquire(file=file, exist=exists)
    if(exists) then
       un = fopen(file)
       call prlog("Found input file for spin parameters: " // trim(file))
       call prlog("Expecting on each line: el Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd")
    else
       call prlog("No spin parameters found, expecting non-spin-polarized calculation.")
       return
    end if

    ! comment line
    read (un, *, iostat=io) dummy

    ! read data
    do
       ! read line and corrent units
       read (un, *, iostat=io) el, W(:)
       W(:) = W(:) * econv

       ! exit if end of file
       if(io /=0) exit

       ! find element
       inserted = .false.
       do i = 1, db%nel
          if (db%e(i)%exists) then
             if(trim(a2s(db%e(i)%name)) == el) then
                db%e(i)%spin = .true.
                db%e(i)%W(0,0) = W(1)
                db%e(i)%W(0,1) = W(2)
                db%e(i)%W(0,2) = W(3)
                db%e(i)%W(1,0) = W(4)
                db%e(i)%W(1,1) = W(5)
                db%e(i)%W(1,2) = W(6)
                db%e(i)%W(2,0) = W(7)
                db%e(i)%W(2,1) = W(8)
                db%e(i)%W(2,2) = W(9)
                inserted = .true.
                exit
             end if
          end if
       end do

       ! report
       if(inserted) then
          call prlog("Read spin parameters for "//trim(el)//" -> "//W(1:9))
       else
          call prlog("For element "//trim(el)//" in file, could not find corresponding database entry.")
       end if

    end do

    call prlog

  end subroutine materials_read_spin_params


  !>
  !! Load the materials database
  !!
  !! Loads the materials database including Slater-Koster tables
  !! for tight-binding. If the directory param/ does not exist,
  !! looks for the tables from the directory pointed to by the
  !! enviroment variable TBPARAM.
  !<
  subroutine materials_read_database(db, econv, lconv, folder, error)
    implicit none

    type(materials_t), intent(inout)    :: db
    real(DP), intent(in)                :: econv, lconv
    character(*), intent(in), optional  :: folder
    integer, intent(inout), optional    :: error

    ! ---

    integer  :: i1, i2, i
    real(DP) :: onsite(9)
    logical  :: params_exists

    ! ---

    call prlog("- materials_read_database -")

    if (present(folder)) then
       db%folder = folder
    else
       call get_environment_variable("TBPARAM", value=db%folder, status=i)
       if (i > 0) then
          db%folder = '.'
       endif
    endif

    call prlog("Looking for tables in directory '"//trim(db%folder)//"'.")
    call prlog

    call materials_init_elements(db, econv, lconv, error)
    PASS_ERROR(error)

    allocate(db%cut(db%nel, db%nel))
    allocate(db%HS(db%nel, db%nel))
    allocate(db%R(db%nel, db%nel))

    call materials_read_sltab_hotbit(db, econv, lconv, error)
    PASS_ERROR(error)
    call materials_read_sltab_dftb(db, econv, lconv, error)
    PASS_ERROR(error)
    call materials_read_sltab_bopfox(db, econv, lconv, error)
    PASS_ERROR(error)

    do i1 = 1, db%nel
       do i2 = 1, db%nel
          if (db%e(i1)%exists .and. db%e(i2)%exists) then
             db%cut(i1, i2) = db%R(i1, i2)%cut
             do i = 1, 2*MAX_NORB
                db%cut(i1, i2) = max(db%cut(i1, i2), db%HS(i1, i2)%cut)
             enddo
          endif
       enddo
    enddo

    call materials_read_spin_params(db, econv, error)

    call prlog("element number charge orbitals Hubbard-U on-site levels")
    call prlog("======= ====== ====== ======== ========= ==============")
    do i1 = 1, db%nel
       if (db%e(i1)%exists) then
          if (db%e(i1)%no > 0 .and. db%e(i1)%no < 10) then
             do i2 = 1, db%e(i1)%no
                onsite(i2) = db%e(i1)%e(get_orbital(db%e(i1)%no, i2))
             enddo
             if (ilog /= -1) then
                write (ilog, '(5X,A7,I7,F7.3,A9,F10.3,9F10.3)') a2s(db%e(i1)%name), db%e(i1)%elem, db%e(i1)%q0, electronic_configuration(db%e(i1)%no), db%e(i1)%U, onsite(1:db%e(i1)%no)
             endif
          else
             if (ilog /= -1) then
                write (ilog, '(5X,A7,I7,F7.3,I9)') a2s(db%e(i1)%name), db%e(i1)%elem, db%e(i1)%q0, db%e(i1)%no
             endif
          endif
       endif
    enddo

    call prlog

  endsubroutine materials_read_database


  !>
  !! Write tables to a file
  !<
  subroutine materials_write_tables(this)
    implicit none

    type(materials_t), intent(in)  :: this

    ! ---

    integer  :: i, j

    ! ---

    do i = 1, this%nel
       do j = 1, this%nel
          if (this%e(i)%exists .and. this%e(j)%exists) then
             if (this%HS(i, j)%n > 0) then
                call write(this%HS(i, j), trim(a2s(this%e(i)%name)) // "-" // trim(a2s(this%e(j)%name)) // "_HS.out")
             endif
             if (this%R(i, j)%n > 0) then
                call write(this%R(i, j), trim(a2s(this%e(i)%name)) // "-" // trim(a2s(this%e(j)%name)) // "_rep.out")
             endif
          endif
       enddo
    enddo

  endsubroutine materials_write_tables
  
endmodule materials
