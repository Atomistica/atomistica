!! ======================================================================
!! Atomistica - Interatomic potential library
!! https://github.com/pastewka/atomistica
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others.
!! See the AUTHORS file in the top-level Atomistica directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
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

  use libAtoms_module

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

  ! IF YOU MODIFY THIS STRUCTURE, *ALWAYS* ALSO MODIFY THE CORRESPONDING
  ! STRUCTURE IN materials.h
  public :: notb_element_t
  type, bind(C) :: notb_element_t

     logical(C_BOOL)         :: exists = .false.

     character(kind=C_CHAR)  :: name(2)=["X","X"]    ! name of element
     character(kind=C_CHAR)  :: cname(10)=["n","o","n","a","m","e"," "," "," "," "]  ! common name of element
     integer(C_INT)          :: elem=10000      ! number of element (official)
     integer(C_INT)          :: no=10000        ! number of orbitals
     integer(C_INT)          :: l(9)=(/0,1,1,1,2,2,2,2,2/) !angular momenta of orbitals
     integer(C_INT)          :: lmax=1000       ! maximum angular momentum
     real(C_DOUBLE)          :: m=1E10          ! mass
     real(C_DOUBLE)          :: e(9)=1E30       ! orbital energies [ e(1:no) ]
     real(C_DOUBLE)          :: occ(9)=1E10     ! occupations in neutral atom
     real(C_DOUBLE)          :: el_max=0        ! max number of valence electrons on an atom
     real(C_DOUBLE)          :: U=1E30          ! Hubbard U
     real(C_DOUBLE)          :: q0=1E30         ! charge (nr of electrons in neutral)
     real(C_DOUBLE)          :: FWHM=1E30       ! .. of Gaussian charge distribution
     real(C_DOUBLE)          :: vib=0d0         ! vibrational frequency of dimer [cm^-1]
     real(C_DOUBLE)          :: Dnn=0d0         ! nearest neighbor distance in bulk 
  
     ! variables for HOTBIT
     real(C_DOUBLE)          :: guess_dq=0d0    ! guessed initial excess charge
     real(C_DOUBLE)          :: fixed_dq=1E10   ! initial charges for TDTB propagation
     integer(C_INT)          :: o1=1E5          ! index of the first orbital
     integer(C_INT)          :: enr=1E5         ! element number in the internal book-keeping

     ! spin-related variables
     logical(C_BOOL)         :: spin = .false.  ! spin-parameters set?
     real(C_DOUBLE)          :: W(0:2,0:2)      ! W parameter values, 0,1,2 = s,p,d, W(0,0) = Wss, W(0,1) = Wsp etc.

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

contains

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

    write (ilog, '(A)')  "- materials_read_sltab_hotbit -"

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
                   write (ilog, '(5X,A)')  "HOTBIT tables for "//trim(e1)//"-"//trim(e2)//" found."
                else if( ex2 ) then
                   un = fopen(fil2)
                   write (ilog, '(5X,A)')  "HOTBIT tables for "//trim(e1)//"-"//trim(e2)//" found."
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

    write (ilog, *)

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

    real(DP), parameter  :: REP_DX    = 0.005
    real(DP), parameter  :: REP_X0    = 1.00

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

    write (ilog, '(A)')  "- materials_read_sltab_dftb -"

    allocate(x(MAX_DATA), y(MAX_DATA))

    conv(HTAB) = econv
    conv(STAB) = 1.0

    do i1 = 1, db%nel
       if (db%e(i1)%exists) then
          do i2 = 1, db%nel
             if (db%e(i2)%exists) then
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
                   write (ilog, '(5X,A)')  "DFTB tables for "//trim(e1)//"-"//trim(e2)//" found."
                else if (ex3) then
                   un = fopen(fil3)
                   write (ilog, '(5X,A)')  "DFTB tables for "//trim(e1)//"-"//trim(e2)//" found."
                endif

                file_exists: if (un > 0) then

                   if (db%HS(i1, i2)%n > 0 .or. db%R(i1, i2)%n > 0) then
                      write (ilog, '(5X,A)')  "WARNING: "//e1//"-"//e2//" tables already read."
                   endif

                   ! cutoff, number of grid points in Hamiltonian/Overlap
                   read (un, *) dx, n
                   n = n-1 ! Sometimes, there seems to be a line thats missing

                   do i = 1, n
                      x(i) = (i-1)*dx
                   enddo

                   if (i1 == i2) then

                      ! self energies, spin polarization energy(?), Hubbard U's, number of electrons
                      read (un, *) eself(:), espin, u(:), q(:)

                      eself  = eself * econv

                      if (i1 == i2) then
                         write (ilog, '(5X,A)')  "WARNING: Overriding self-energies and Hubbard-U from 'elements.dat'."

                         db%e(i1)%e  = (/ &
                              eself(3), &
                              eself(2), eself(2), eself(2), &
                              eself(1), eself(1), eself(1), eself(1), eself(1) &
                              /)

                         db%e(i1)%U  = u(3) * econv
                      endif

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

    write (ilog, *)

  endsubroutine materials_read_sltab_dftb


  !>
  !! Reads element data from elements.dat
  !<
  subroutine materials_read_elements(db, econv, lconv, error)
    implicit none

    type(materials_t), intent(inout) :: db
    real(DP),          intent(in)    :: econv, lconv
    integer, optional, intent(out)   :: error

    ! ---

    integer               :: i,j,k,un,io,lmx(9)
    character(200)        :: line,dat,key,fn
    logical               :: ex

    type(notb_element_t)  :: hlp(MAX_Z)

    ! ---

    INIT_ERROR(error)

    write (ilog, '(A)')  "- materials_read_elements -"

    lmx = 1000
    lmx(1)=0; lmx(4)=1; lmx(9)=2 !lmax = lmx(no)

    db%nel = MAX_Z
    allocate(db%e(db%nel))

    fn = trim(db%folder) // "/elements.dat"

    inquire(file=fn, exist=ex)

    if (.not. ex) then
       RAISE_ERROR("ERROR: Could not open '" // trim(fn) // "'.", error)
    endif

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
          read(dat,*) hlp(j)%name
       case("Z");     read(dat,*) hlp(j)%elem
       case("common");read(dat,*) hlp(j)%cname
       case("m");     read(dat,*) hlp(j)%m
       case("q0");    read(dat,*) hlp(j)%q0
       case("no")
          read(dat,*) hlp(j)%no
          hlp(j)%lmax = lmx(hlp(j)%no)
          read(un ,*) hlp(j)%e(1:hlp(j)%no)
          hlp(j)%e(1:hlp(j)%no) = hlp(j)%e(1:hlp(j)%no) * econv
       case("U")
          read(dat,*) hlp(j)%U
          hlp(j)%U = hlp(j)%U * econv
       case("FWHM");  read(dat,*) hlp(j)%FWHM
       case("vib");   read(dat,*) hlp(j)%vib
       case("Dnn");   read(dat,*) hlp(j)%Dnn
       case("occ");   read(dat,*) hlp(j)%occ(1:hlp(j)%no)
       case default
          cycle
       end select
    end do

    call fclose(un)

    ! set dependent variables and sort according to Z
    do i=1,j

       if (hlp(i)%elem > 0 .and. hlp(i)%elem <= MAX_Z) then
          if (trim(a2s(hlp(i)%name)) /= trim(ElementName(hlp(i)%elem))) then
             write (ilog, '(5X,5A)')  "WARNING: Name '", hlp(i)%name, "' in 'elements.dat' not equal common element name '", hlp(i)%elem, "'."
          endif

          hlp(i)%el_max = 0d0
          do k=1,hlp(j)%no
             hlp(i)%el_max = hlp(i)%el_max + (2d0*hlp(i)%l(k)+1d0)
          end do
          db%e(hlp(i)%elem) = hlp(i)
          db%e(hlp(i)%elem)%exists = .true.
          db%e(hlp(i)%elem)%enr = hlp(i)%elem
       else
          write (ilog, '(5X,A,I5,A)')  "WARNING: Unknown element found in 'elements.dat' (Z = ", hlp(i)%elem, ")."
       endif

    end do

    write (ilog, '(5X,I5,A)')  j, " elements found in 'elements.dat'."
    write (ilog, *)

  endsubroutine materials_read_elements


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

    write (ilog, '(A)')  "- materials_read_spin_params -"

    ! reset
    do i = 1, db%nel
       db%e(i)%W(:,:) = 0.0_DP
    end do

    ! open file and return if none found
    file = trim(db%folder)//'/'//trim(file)
    inquire(file=file, exist=exists)
    if(exists) then
       un = fopen(file)
       write (ilog, '(5X,A)')  "Found input file for spin parameters: " // trim(file)
       write (ilog, '(5X,A)')  "Expecting on each line: el Wss Wsp Wsd Wps Wpp Wpd Wds Wdp Wdd"
    else
       write (ilog, '(5X,A)')  "No spin parameters found, expecting non-spin-polarized calculation."
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
          write (ilog, '(5X,A,9F10.3)')  "Read spin parameters for " // trim(el) // " -> ", W(1:9)
       else
          write (ilog, '(5X,A)')  "For element " // trim(el) // " in file, could not find corresponding database entry."
       end if

    end do

    write (ilog, *)

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
    logical  :: params_exists

    ! ---

    write (ilog, '(A)')  "- materials_read_database -"

    if (present(folder)) then
       db%folder = folder
    else
       ! Get input directory
       inquire(file='param/elements.dat', exist=params_exists)
       if(params_exists) then
          db%folder = "param"
       else
          call getenv("TBPARAM", db%folder)
          db%folder = trim(db%folder)
       endif
    endif

    write (ilog, '(5X,A)')  "Looking for tables in directory '", trim(db%folder), "'."
    write (ilog, *)         ""

    call materials_read_elements(db, econv, lconv, error)
    PASS_ERROR(error)

    allocate(db%cut(db%nel, db%nel))
    allocate(db%HS(db%nel, db%nel))
    allocate(db%R(db%nel, db%nel))

    call materials_read_sltab_hotbit(db, econv, lconv, error)
    PASS_ERROR(error)
    call materials_read_sltab_dftb(db, econv, lconv, error)
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

    write (ilog, *)

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
