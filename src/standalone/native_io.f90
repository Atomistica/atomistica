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
!! I/O, native MDCORE format
!!
!! I/O, native MDCORE format
!<

#include "macros.inc"

module native_io
#ifdef _MP
  use mpi
#endif

  use libAtoms_module

  use io
  use logging
  use misc

  use timer

  use data
  use particles
  use cyclic
  use molecules

  private

  character(*), parameter  :: MODULE_STR  = "NativeIO"

  character(MAX_NAME_STR), parameter  :: T_STR            = "temperatures"
  character(MAX_NAME_STR), parameter  :: DISSIPATION_STR  = "langevin_dissipation"

  public :: read_atoms
  interface read_atoms
     module procedure native_io_read_atoms
  endinterface

  public :: write_atoms 
  interface write_atoms
     module procedure native_io_write_atoms
  endinterface

  public :: read_Z_and_groups_from_atoms
  public :: read_cell_from_atoms
  
contains


  !>
  !! Read the particle positions, etc. from an atoms.dat file
  !!
  !! Read the particle positions, etc. from an atoms.dat file
  !<
  subroutine native_io_read_atoms(p, fn, mol, skip_cell, allow_def, error)
    implicit none

    type(particles_t), intent(inout)            :: p
    character(*), intent(in)                    :: fn
    type(molecules_t), intent(inout), optional  :: mol
    logical, intent(in), optional               :: skip_cell
    logical, intent(in), optional               :: allow_def
    integer, intent(out), optional              :: error

    ! ---

    integer               :: un, i, j, k, l, nat, stat, wc, data_type, findex
    integer               :: next, nat_not_fixed, totnat, Z
    real(DP)              :: r3(3), cur_diss, cur_T

    character(1000)       :: line

    logical, allocatable  :: found_real(:)
    logical, allocatable  :: found_integer(:)
    logical, allocatable  :: found_real3(:)
    logical, allocatable  :: found_real3x3(:)

#ifdef _MP
    real(DP)           :: r, r3x3(3, 3)
    type(MPI_context)  :: mpi
#endif

    real(DP), pointer   :: v(:, :)
    real(DP), pointer   :: f(:, :)
    real(DP), pointer   :: T(:)
    real(DP), pointer   :: dissipation(:)

    ! ---

    INIT_ERROR(error)
    DEBUG_WRITE("===> read_atoms " // fn)

    call prlog("- read_atoms -")
    if(.not. present(mol)) then
       call prlog("     No molecules object specified, ignoring molecule information (next)")
    end if

    if (.not. initialized(p)) then
       call init(p)
    endif

    v            => NULL()
    f            => NULL()
    T            => NULL()
    dissipation  => NULL()

    allocate(found_real(p%data%n_real))
    allocate(found_integer(p%data%n_integer))
    allocate(found_real3(p%data%n_real3))
    allocate(found_real3x3(p%data%n_real3x3))

    found_real     = .false.
    found_integer  = .false.
    found_real3    = .false.
    found_real3x3  = .false.

    un = fopen(fn, F_READ, error=error)
    PASS_ERROR(error)

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    read(un, *, iostat=stat)  nat
    if (stat /= 0) then
       RAISE_ERROR("Could not read number of atoms.", error)
    endif

    if (.not. allocated(p)) then
       if (present(mol)) then
          call register_data(mol, p)
       endif

       call allocate(p, nat, allow_def=allow_def, error=error)
       PASS_ERROR(error)
    else
       call set_total_nat(p, nat)
    endif

    if (exists(p%data, V_STR)) then
       call ptr_by_name(p%data, V_STR, v)
    endif
    if (exists(p%data, F_STR)) then
       call ptr_by_name(p%data, F_STR, f)
    endif
    if (exists(p%data, T_STR)) then
       call ptr_by_name(p%data, T_STR, T)
    endif
    if (exists(p%data, DISSIPATION_STR)) then
       call ptr_by_name(p%data, DISSIPATION_STR, dissipation)
    endif

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("Could not read occupation information.", error)
    endif

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    nat_not_fixed  = 0

    DEBUG_WRITE("Start reading positions")

    j   = 1
    wc  = 0
    do i = 1, nat
       if (j > p%maxnatloc) then
           RAISE_ERROR("Particles object too small, maxnatloc exceeded (j = " // j // ", maxnatloc = " // p%maxnatloc // ").", error)
       endif

       read (un, '(A)', iostat=stat)  line
       if (stat /= 0) then
          RAISE_ERROR("Error reading from file (Premature end-of-file?).", error)
       endif

       read(line, *, iostat=stat)  p%sym(j), p%m(j), r3, p%g(j), cur_diss, cur_T, findex, next
       if (stat /= 0) then
          findex = i
          read(line, *, iostat=stat)  p%sym(j), p%m(j), r3, p%g(j), cur_diss, cur_T, next

          if (stat /= 0) then
             next = 0
             read(line, *, iostat=stat)  p%sym(j), p%m(j), r3, p%g(j), cur_diss, cur_T

             if (stat /= 0) then
                RAISE_ERROR("Error parsing element info/positions information (i = " // i // ", j = " // j // ", maxnatloc = " // p%maxnatloc // ", line = '" // trim(line) // "')", error)
             endif
          endif
       endif

       if (associated(dissipation)) then
          dissipation(j)  = cur_diss
       endif
       if (associated(T)) then
          T(j)  = cur_T
       endif

       if (findex /= i) then
          write (*, '(A,I10)')  "findex  = ", findex
          write (*, '(A,I10)')  "i       = ", i
          RAISE_ERROR("Wrong index in input file (findex /= i, i = " // i // ", findex = " // findex // ".", error)
       endif

!       r3  = cyclic_in_cell(p, r3)

       ! next array only taken into account if molecules object present
       if (present(mol)) then
          call molecules_verify(mol, p)
          mol%next(j)   = next
       endif
       p%index(j)  = findex

#ifndef IMPLICIT_R
       POS3(p, j)         = r3
#endif
       PNC3(p, j)         = r3
       VEC3(p%r_cont, j)  = r3

       Z = atomic_number(p%sym(j))
       if (Z > 0 .and. Z <= MAX_Z) then

          p%m(j) = ElementMass(Z)
          p%Z(j) = Z

       else
          call prlog("     m    = "//p%m(j))
          call prlog("     sym  = "//p%sym(j))
          RAISE_ERROR("Mass negative or equal to zero and atom symbol '" // p%sym(j) // "' unknown.", error)
       endif

       ! next array only taken into account if molecules object present
       if(present(mol)) then
          if (next <= 0) then

             ! -------------------------
             ! Autodetect water
             !
             if (p%sym(j) == "O") then
                wc = wc+1
             else if (p%sym(j) == "H" .and. (wc == 1 .or. wc == 2)) then
                wc = wc+1
             else
                wc = 0
             endif

             if (wc == 3) then
                ! This is an O-H-H, a water

                if (p%global2local(i-2) > 0) then
                   mol%next(p%global2local(i-2)) = i-1
                endif

                if (p%global2local(i-1) > 0) then
                   mol%next(p%global2local(i-1)) = i
                endif

                wc = 0
             endif
             !
             ! -------------------------

          endif
       endif

#ifdef _MP
       if (all(r3 >= p%lower) .and. all(r3 < p%upper)) then
#endif

       if (i > p%totnat) then
           RAISE_ERROR("global2local index too small, totnat exceeded (i = " // i // ", totnat = " // p%totnat // ").", error)
       endif

       p%global2local(i) = j

       if (p%g(j) <= 0) then
          if (associated(v)) &
               VEC3(v, j) = 0.0_DP
          if (associated(f)) &
               VEC3(f, j) = 0.0_DP
       else
          nat_not_fixed = nat_not_fixed+1
       endif

       j = j+1

#ifdef _MP
       endif
#endif

    enddo

    DEBUG_WRITE("Done reading positions")

    p%nat     = j-1
    p%natloc  = j-1

    totnat    = j-1

#ifdef _MP
    DEBUG_WRITE("dmp_sum")
    call initialise(mpi)
    call sum_in_place(mpi, totnat)
    call finalise(mpi)
    DEBUG_WRITE("p%totnat = " // p%totnat)
#endif

    p%dof = 3*p%totnat-3

    call prlog("     nat        = "//p%nat)
    call prlog("     maxnatloc  = "//p%maxnatloc)
    call prlog("     totnat     = "//totnat)

    if (totnat /= nat) then
       RAISE_ERROR("Something wrong: Different number of particles (totnat = " // totnat // ") loaded than defined in the input (nat = " // nat // "). Are some particles outside of the simulation cell?", error)
    endif

    read (un, '(A)', iostat=stat)  line
    do while (stat == 0 .and. len_trim(line) == 0)
       read (un, '(A)', iostat=stat)  line
    enddo
    do while (stat == 0)
       do while ((line(1:1) == '<' .or. line(1:1) == '-' .or. line(1:1) == '#' .or. line(1:1) == ' ') .and. len_trim(line) > 0)
          line  = line(2:)
       enddo

       DEBUG_WRITE(line)

       i  = index(line, '=')
       if (i == 0) then

          !
          ! This is a field or an attribute
          !

          ! FIXME!!! Dirty. Skip cell information because this has already been read.
          if (.not. (present(skip_cell) .and. equal(line, "cell")) .and. exists(p%data, line, data_type)) then

             select case (data_type)

             case (TYPE_REAL_ATTR)
                i = index_by_name(p%data%n_real_attr, p%data%name_real_attr, trim(line))

                read (un, *)  p%data%data_real_attr(i)

             case (TYPE_REAL3_ATTR)
                i = index_by_name(p%data%n_real3_attr, p%data%name_real3_attr, trim(line))

                read (un, *)  p%data%data_real3_attr(:, i)

             case (TYPE_REAL3x3_ATTR)
                i = index_by_name(p%data%n_real3x3_attr, p%data%name_real3x3_attr, trim(line))

                read (un, *)  ( p%data%data_real3x3_attr(:, j, i), j = 1, 3 )

             case (TYPE_REAL)
                i = index_by_name(p%data%n_real, p%data%name_real, trim(line))

#ifdef _MP
                do j = 1, p%totnat
                   read (un, *)  r
                   if (p%global2local(j) > 0) then
                      p%data%data_real(p%global2local(j), i)  = r
                   endif
                enddo
#else
                read (un, *)  ( p%data%data_real(j, i), j = 1, nat )
#endif

                found_real(i)  = .true.

             case (TYPE_INTEGER)
                i = index_by_name(p%data%n_integer, p%data%name_integer, trim(line))

#ifdef _MP
                do j = 1, p%totnat
                   read (un, *)  r
                   if (p%global2local(j) > 0) then
                      p%data%data_integer(p%global2local(j), i)  = r
                   endif
                enddo
#else
                read (un, *)  ( p%data%data_integer(j, i), j = 1, nat )
#endif

                found_integer(i)  = .true.

             case (TYPE_REAL3)
                i = index_by_name(p%data%n_real3, p%data%name_real3, trim(line))

#ifdef _MP
                do j = 1, p%totnat
                   read (un, *)  r3
                   if (p%global2local(j) > 0) then
                      p%data%data_real3(:, p%global2local(j), i)  = r3
                   endif
                enddo
#else

                read (un, *)  ( ( p%data%data_real3(k, j, i), k = 1, 3 ), j = 1, nat )

#endif

                found_real3(i)  = .true.

             case (TYPE_REAL3x3)
                i = index_by_name(p%data%n_real3x3, p%data%name_real3x3, trim(line))

#ifdef _MP
                do j = 1, p%totnat
                   read (un, *)  ( ( r3x3(k, l), k = 1, 3 ), l = 1, 3 )
                   if (p%global2local(j) > 0) then
                      p%data%data_real3x3(:, :, p%global2local(j), i)  = r3x3(:, :)
                   endif
                enddo
#else
                read (un, *)  ( ( ( p%data%data_real3x3(k, l, j, i), k = 1, 3 ), l = 1, 3 ), j = 1, nat )
#endif

                found_real3x3(i)  = .true.

             case default
                RAISE_ERROR("Don't know how to read data type of field/attribute '" // trim(line) // "." , error)
             endselect

             read (un, '(A)', iostat=stat)  line
             do while (stat == 0 .and. len_trim(line) == 0)
                read (un, '(A)', iostat=stat)  line
             enddo

          else

             if (.not. equal(line, "cell")) then

#ifdef _MP
                if (mpi_id() == ROOT) then
#endif
                   WARN("Undefined field/attribute '" // trim(line) // "' found in input file.")
#ifdef _MP
                endif
#endif
                WARN("Undefined field/attribute '" // trim(line) // "' found in input file.")

             endif

             read (un, '(A)', iostat=stat)  line
             line = adjustl(line)
             do while (stat == 0 .and. line(1:1) /= '<' .and. line(1:1) /= '#')
                read (un, '(A)', iostat=stat)  line
                line = adjustl(line)
             enddo

          endif

       else

          !
          ! This is an attribute
          !

          if (exists(p%data, trim(line(1:i-1)), data_type)) then

             select case (data_type)

             case (TYPE_REAL_ATTR)
                j = index_by_name(p%data%n_real_attr, p%data%name_real_attr, trim(line(1:i-1)))

                read (line(i+1:), *)  p%data%data_real_attr(j)

             case (TYPE_REAL3_ATTR)
                j = index_by_name(p%data%n_real3_attr, p%data%name_real3_attr, trim(line(1:i-1)))

                read (line(i+1:), *)  p%data%data_real3_attr(:, j)

             case (TYPE_REAL3x3_ATTR)
                j = index_by_name(p%data%n_real3x3_attr, p%data%name_real3x3_attr, trim(line(1:i-1)))

                read (line(i+1:), *) &
                     p%data%data_real3x3_attr(:, 1, j), &
                     p%data%data_real3x3_attr(:, 2, j), &
                     p%data%data_real3x3_attr(:, 3, j)

             case default
                RAISE_ERROR("Don't know how to read data type of field '" // trim(line) // "." , error)
             endselect

             read (un, '(A)', iostat=stat)  line
             do while (stat == 0 .and. len_trim(line) == 0)
                read (un, '(A)', iostat=stat)  line
             enddo

          else

             WARN("Warning: Undefined attribute '" // trim(line(1:i-1)) // "' found in input file.")

             read (un, '(A)', iostat=stat)  line
             line = adjustl(line)
             do while (stat == 0 .and. line(1:1) /= '<' .and. line(1:1) /= '#')
                read (un, '(A)', iostat=stat)  line
                line = adjustl(line)
             enddo

          endif
      
       endif

    enddo

    call fclose(un)

    if (.not. present(skip_cell)) then
       ! Initialize reciprocal lattice vectors
       call set_cell(p, p%Abox, error=error)
       PASS_ERROR(error)
    endif

    call update_elements(p)

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_RESTART) /= 0 .and. .not. found_real(i)) then
          WARN("Field '" // trim(p%data%name_real(i)) // "' is requested, however this field was not found in the input file.")
       endif
    enddo

    do i = 1, p%data%n_integer
       if (iand(p%data%tag_integer(i), F_RESTART) /= 0 .and. .not. found_integer(i)) then
          WARN("Field '" // trim(p%data%name_integer(i)) // "' is requested, however this field was not found in the input file.")
       endif
    enddo

    do i = 1, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_RESTART) /= 0 .and. .not. found_real3(i)) then
          WARN("Field '" // trim(p%data%name_real3(i)) // "' is requested, however this field was not found in the input file.")
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_RESTART) /= 0 .and. .not. found_real3x3(i)) then
          WARN("Field '" // trim(p%data%name_real3x3(i)) // "' is requested, however this field was not found in the input file.")
       endif
    enddo

    deallocate(found_real)
    deallocate(found_integer)
    deallocate(found_real3)
    deallocate(found_real3x3)

    if (any(p%shear_dx /= 0.0_DP)) then
       call set_lees_edwards(p, p%shear_dx, error=error)
       PASS_ERROR(error)
       call prlog("     shear_dx   = "//p%shear_dx)
    endif

    call prlog

    DEBUG_WRITE("<=== read_atoms")

  endsubroutine native_io_read_atoms


  !**********************************************************************
  ! Read the particle positions, etc. from an atoms.dat file
  !**********************************************************************    
  subroutine read_cell_from_atoms(p, fn, allow_def, error)
    implicit none

    type(particles_t), intent(inout)  :: p
    character(*), intent(in)          :: fn
    logical, intent(in), optional     :: allow_def
    integer, intent(inout), optional  :: error

    ! ---

    integer            :: un, j, nat, stat

    character(1000)    :: line

    ! ---

    call prlog("- read_cell_from_atoms -")

    un = fopen(fn, F_READ, error=error)
    PASS_ERROR(error)

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    read(un, *, iostat=stat)  nat
    if (stat /= 0) then
       RAISE_ERROR("Could not read number of atoms.", error)
    endif

    if (.not. allocated(p)) then
       call allocate(p, nat, allow_def=allow_def, error=error)
       PASS_ERROR(error)
    endif

    read (un, '(A)', iostat=stat)  line
    do while (stat == 0 .and. len_trim(line) == 0)
       read (un, '(A)', iostat=stat)  line
    enddo
    do while (stat == 0)
       do while ((line(1:1) == '<' .or. line(1:1) == '-' .or. line(1:1) == '#' .or. line(1:1) == ' ') .and. len_trim(line) > 0)
          line  = line(2:)
       enddo

!       write (*, *)  trim(line)

       if (equal(line, "cell")) then
          read (un, *)  ( p%Abox(:, j), j = 1, 3 )

           read (un, '(A)', iostat=stat)  line
          do while (stat == 0 .and. len_trim(line) == 0)
             read (un, '(A)', iostat=stat)  line
          enddo
       else
          read (un, '(A)', iostat=stat)  line
          line = adjustl(line)
          do while (stat == 0 .and. line(1:1) /= '<' .and. line(1:1) /= '#')
             read (un, '(A)', iostat=stat)  line
             line = adjustl(line)
          enddo
       endif
    enddo

    call fclose(un)

    call set_cell(p, p%Abox, error=error)
    PASS_ERROR(error)

    call prlog

  endsubroutine read_cell_from_atoms


  !**********************************************************************
  ! Read only Z, groups from atoms.dat
  !**********************************************************************    
  subroutine read_Z_and_groups_from_atoms(p, fn, mol, error)
    implicit none

    type(particles_t), intent(inout)            :: p
    character(*), intent(in)                    :: fn
    type(molecules_t), intent(inout), optional  :: mol
    integer, intent(out), optional              :: error

    ! ---

    integer            :: un, i, j, nat, stat, wc, Z
    real(DP)           :: r(3), cur_diss, cur_T

    real(DP), pointer  :: T(:)
    real(DP), pointer  :: dissipation(:)

    ! ---

    INIT_ERROR(error)

    call prlog("- read_Z_and_groups_from_atoms -")

    T            => NULL()
    dissipation  => NULL()
    if (exists(p%data, T_STR)) then
       call ptr_by_name(p%data, T_STR, T)
    endif
    if (exists(p%data, DISSIPATION_STR)) then
       call ptr_by_name(p%data, DISSIPATION_STR, dissipation)
    endif

    un = fopen(fn, F_READ, error=error)
    PASS_ERROR(error)

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    read(un, *, iostat=stat)  nat
    if (stat /= 0) then
       RAISE_ERROR("Could not read number of atoms.", error)
    endif

    if (p%nat /= nat) then
       RAISE_ERROR("Number of atoms from '" // trim(fn) // "' does not match current number of atoms.", error)
    endif

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("Could not read occupation information.", error)
    endif

    read(un, *, iostat=stat)
    if (stat /= 0) then
       RAISE_ERROR("End-of-file reached while reading '" // fn // "'.", error)
    endif

    j   = 1
    wc  = 0
    do i = 1, nat
       read(un, *, iostat=stat)  p%sym(j), p%m(j), r, p%g(j), cur_diss, cur_T
       if (stat /= 0) then
          RAISE_ERROR("Error reading positions.", error)
       endif

       if (associated(dissipation)) then
          dissipation(j)  = cur_diss
       endif
       if (associated(T)) then
          T(j)  = cur_T
       endif

       p%index(j)  = i

       Z = atomic_number(p%sym(j))
       if (Z > 0 .and. Z <= MAX_Z) then

          p%m(j) = ElementMass(Z)
          p%Z(j) = Z

       else
          call prlog("     m    = "//p%m(j))
          call prlog("     sym  = "//p%sym(j))
          RAISE_ERROR("Mass negative or equal to zero and atom symbol unknown.", error)
       endif

       ! -------------------------
       ! Autodetect water
       !
       if (p%sym(j) == "O") then
          wc = wc+1
       else if (p%sym(j) == "H" .and. (wc == 1 .or. wc == 2)) then
          wc = wc+1
       else
          wc = 0
       endif

       if (present(mol) .and. wc == 3) then
          ! This is an O-H-H, a water

          if (p%global2local(i-2) > 0) then
             mol%next(p%global2local(i-2)) = i-1
          endif

          if (p%global2local(i-1) > 0) then
             mol%next(p%global2local(i-1)) = i
          endif

          wc = 0
       endif
       !
       ! -------------------------

       j  = j + 1

    enddo

    call fclose(un)

    call update_elements(p)

    call prlog

  endsubroutine read_Z_and_groups_from_atoms

#ifdef _MP

  !>
  !! Write particles to native file format (parallel version)
  !!
  !! Write the particles to an atoms.out file (parallel version)
  !<
  subroutine native_io_write_atoms(this, fn, mol, error)
    implicit none

    type(particles_t), intent(in)     :: this
    character(*), intent(in)          :: fn
    type(molecules_t), intent(in), optional :: mol
    integer, intent(out), optional    :: error

    ! ---

    call internal_native_io_write_atoms(this, fn, mpi_n_procs(), mol, error)
    PASS_ERROR(error)

  endsubroutine native_io_write_atoms

  !>
  !! Write particles to native file format (parallel version)
  !!
  !! Write the particles to an atoms.out file (parallel version)
  !<
  subroutine internal_native_io_write_atoms(this, fn, mpi_n_procs, mol, error)
    implicit none

    type(particles_t), intent(in)     :: this
    character(*), intent(in)          :: fn
    integer, intent(in)               :: mpi_n_procs
    type(molecules_t), intent(in), optional :: mol
    integer, intent(out), optional    :: error

    ! ---

    integer, parameter  :: MAX_BUFFER_SIZE = 1000
    integer, parameter  :: NEXT_SECTION = -1

    integer             :: un, i, j, k, g, p, buffer_size, curi, ierr

    integer             :: status(MPI_STATUS_SIZE)

    integer             :: buffer_pos(mpi_n_procs-1)
    integer             :: n(mpi_n_procs-1)
    integer             :: indbuf(MAX_BUFFER_SIZE, mpi_n_procs-1)
    integer             :: ibuf(MAX_BUFFER_SIZE, mpi_n_procs-1)
    real(DP)            :: xbuf(MAX_BUFFER_SIZE, mpi_n_procs-1)
    real(DP)            :: ybuf(MAX_BUFFER_SIZE, mpi_n_procs-1)
    real(DP)            :: zbuf(MAX_BUFFER_SIZE, mpi_n_procs-1)

    ! Specialized buffer
    real(DP)            :: mass(MAX_BUFFER_SIZE, mpi_n_procs-1)
    integer             :: group(MAX_BUFFER_SIZE, mpi_n_procs-1)
    real(DP)            :: dissipation(MAX_BUFFER_SIZE, mpi_n_procs-1)
    real(DP)            :: T(MAX_BUFFER_SIZE, mpi_n_procs-1)
    integer             :: next(MAX_BUFFER_SIZE, mpi_n_procs-1)

    integer             :: cur_Z, cur_g, cur_n
    real(DP)            :: cur_m, cur_x, cur_r(3), cur_d, cur_T

    logical             :: found

    real(DP), pointer   :: this_v(:, :)
    real(DP), pointer   :: this_f(:, :)
    real(DP), pointer   :: this_T(:)
    real(DP), pointer   :: this_dissipation(:)

    ! ---

    INIT_ERROR(error)

    call timer_start("write_atoms")

    if (ROOT /= 0) then
       stop "Assuming ROOT = 0"
    endif

    this_v            => NULL()
    this_f            => NULL()
    this_T            => NULL()
    this_dissipation  => NULL()
    if (exists(this%data, V_STR)) then
       call ptr_by_name(this%data, V_STR, this_v)
    endif
    if (exists(this%data, F_STR)) then
       call ptr_by_name(this%data, F_STR, this_f)
    endif
    if (exists(this%data, T_STR)) then
       call ptr_by_name(this%data, T_STR, this_T)
    endif
    if (exists(this%data, DISSIPATION_STR)) then
       call ptr_by_name(this%data, DISSIPATION_STR, this_dissipation)
    endif

    if (mpi_id() == ROOT) then
       un = fopen(fn, F_WRITE, error=error)
       PASS_ERROR(error)

       write (un, '(A)')  "<--- Total number of atoms"

       write (un, '(I20)')  this%totnat

       write (un, '(A)')  "<--- *** The following line is ignored ***"
       write (un, *)  
       write (un, '(A)')  "<--- Element, atomic mass, coordinates, group, dissipation, temperature, (next)"

       n           = MAX_BUFFER_SIZE
       buffer_pos  = n + 1

       do g = 1, this%totnat

          do p = 1, mpi_n_procs-1
             if (n(p) == MAX_BUFFER_SIZE .and. buffer_pos(p) > n(p)) then

                call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
                PASS_MPI_ERROR(ierr, error)

                call mpi_recv(indbuf(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)

                call mpi_get_count(status, MPI_INTEGER, n(p), ierr)
                PASS_MPI_ERROR(ierr, error)

                call mpi_recv(ibuf(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)
                call mpi_recv(mass(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)
                call mpi_recv(xbuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)
                call mpi_recv(ybuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)
                call mpi_recv(zbuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)
                call mpi_recv(group(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)

                dissipation(:, p)  = 0.0_DP
                if (associated(this_dissipation)) then
                   call mpi_recv(dissipation(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
                   PASS_MPI_ERROR(ierr, error)
                endif
                T(:, p)  = 0.0_DP
                if (associated(this_T)) then
                   call mpi_recv(T(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
                   PASS_MPI_ERROR(ierr, error)
                endif
                ! Tommi: XXX
                call mpi_recv(next(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
                PASS_MPI_ERROR(ierr, error)

                buffer_pos(p)  = 1

             endif
          enddo
       
          i  = this%global2local(g)

          if (i > 0 .and. i <= this%natloc) then
             ! Okay, this particle is on this processor

             cur_Z  = this%Z(i)
             cur_m  = this%m(i)
             cur_r  = POS3(this, i)
             cur_g  = this%g(i)
             cur_d  = 0.0_DP
             if (associated(this_dissipation)) then
                cur_d  = this_dissipation(i)
             endif
             cur_T  = 0.0_DP
             if (associated(this_T)) then
                cur_T  = this_T(i)
             endif
             if(present(mol)) then
                cur_n  = mol%next(i)
             else
                cur_n = 0
             end if

          else
             ! Particle must be somewhere else
             
             found = .false.

             do p = 1, mpi_n_procs-1
                if (n(p) > 0) then

                   j  = buffer_pos(p)

                   if (.not. (j < 1 .or. j > n(p) .or. j > MAX_BUFFER_SIZE)) then

                      if (j <= n(p) .and. indbuf(j, p) == g) then
                         if (found) then
                            RAISE_ERROR("Particle was found twice.", error)
                         endif

                         found = .true.

                         cur_Z  = ibuf(j, p)
                         cur_m  = mass(j, p)
                         cur_r  = (/ xbuf(j, p), ybuf(j, p), zbuf(j, p) /)
                         cur_g  = group(j, p)
                         cur_d  = dissipation(j, p)
                         cur_T  = T(j, p)
                         cur_n  = next(j, p)
                         
                         buffer_pos(p)  = j + 1

                      endif
                   endif

                endif
             enddo

             if (.not. found) then
                RAISE_ERROR("Particle not found.", error)
             endif

          endif

          if (cur_Z > 0 .and. cur_Z <= MAX_Z) then
             if (cur_n > 0) then
                write(un, '(1X,A4,4ES20.10,I5,2ES20.10,I10)') &
                     ElementName(cur_Z), cur_m, cur_r, &
                     cur_g, cur_d, cur_T, cur_n
             else
                write(un, '(1X,A4,4ES20.10,I5,2ES20.10)') &
                     ElementName(cur_Z), cur_m, cur_r, &
                     cur_g, cur_d, cur_T
             endif
          else
             RAISE_ERROR("Unknown atomic number encountered.", error)
          endif
       enddo

       call mpi_bcast(NEXT_SECTION, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
       PASS_MPI_ERROR(ierr, error)

       if (associated(this_v)) then
          write (un, '(A)')  "<--- Velocities"

          n           = MAX_BUFFER_SIZE
          buffer_pos  = n + 1

          do g = 1, this%totnat
             
             call recv_real3

             i  = this%global2local(g)

             if (i > 0 .and. i <= this%natloc) then
                ! Okay, this particle is on this processor

                cur_r  = VEC3(this_v, i)

             else
                ! Particle must be somewhere else
             
                call find_real3(g, cur_r)

             endif

             write(un, '(1X,3ES20.10)')  cur_r
          enddo

          call mpi_bcast(NEXT_SECTION, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
          PASS_MPI_ERROR(ierr, error)
       endif

       write (un, '(A)')  "<--- Forces"

       n           = MAX_BUFFER_SIZE
       buffer_pos  = n + 1

       do g = 1, this%totnat

          call recv_real3

          i  = this%global2local(g)

          if (i > 0 .and. i <= this%natloc) then
             ! Okay, this particle is on this processor

             cur_r  = 0.0_DP
             if (associated(this_f)) &
                  cur_r  = VEC3(this_f, i)

          else
             ! Particle must be somewhere else
             
             call find_real3(g, cur_r)

          endif

          write(un, '(1X,3ES20.10)')  cur_r
       enddo

       call mpi_bcast(NEXT_SECTION, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
       PASS_MPI_ERROR(ierr, error)

       !
       ! Everything else ... dynamic!
       !

       do i = 1, this%data%n_real_attr
          write (un, '(1X,A,A)')      "<--- ", trim(this%data%name_real_attr(i))
          write (un, '(1X,ES20.10)')  this%data%data_real_attr(i)
       enddo

       do i = 1, this%data%n_real3_attr
          write (un, '(1X,A,A)')       "<--- ", trim(this%data%name_real3_attr(i))
          write (un, '(1X,3ES20.10)')  this%data%data_real3_attr(:, i)
       enddo

       do i = 1, this%data%n_real3x3_attr
          write (un, '(1X,A,A)')       "<--- ", trim(this%data%name_real3x3_attr(i))
          write (un, '(1X,3ES20.10)')  ( this%data%data_real3x3_attr(:, j, i), j = 1, 3 )
       enddo

       !
       ! Arrays
       !

       do k = 1, this%data%n_real
          if (iand(this%data%tag_real(k), F_RESTART) /= 0) then
             write (un, '(1X,A,A)')  "<--- ", trim(this%data%name_real(k))

             n           = MAX_BUFFER_SIZE
             buffer_pos  = n + 1

             do g = 1, this%totnat

                call recv_real

                i  = this%global2local(g)

                if (i > 0 .and. i <= this%natloc) then
                   ! Okay, this particle is on this processor

                   cur_x  = this%data%data_real(i, k)

                else
                   ! Particle must be somewhere else
             
                   call find_real(g, cur_x)

                endif

                write(un, '(1X,ES20.10)')  cur_x
             enddo

             call mpi_bcast(NEXT_SECTION, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

          endif
       enddo

       do k = 1, this%data%n_integer
          if (iand(this%data%tag_integer(k), F_RESTART) /= 0) then
             write (un, '(1X,A,A)')  "<--- ", trim(this%data%name_integer(k))

             n           = MAX_BUFFER_SIZE
             buffer_pos  = n + 1

             do g = 1, this%totnat

                call recv_integer

                i  = this%global2local(g)

                if (i > 0 .and. i <= this%natloc) then
                   ! Okay, this particle is on this processor

                   cur_Z  = this%data%data_integer(i, k)

                else
                   ! Particle must be somewhere else
             
                   call find_integer(g, cur_Z)

                endif

                write(un, '(1X,I10)')  cur_Z
             enddo

             call mpi_bcast(NEXT_SECTION, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

          endif
       enddo

       do k = 1, this%data%n_real3
          if (iand(this%data%tag_real3(k), F_RESTART) /= 0 .and. &
              uppercase(trim(this%data%name_real3(k))) /= "VELOCITIES" .and. &
              uppercase(trim(this%data%name_real3(k))) /= "FORCES") then
             write (un, '(1X,A,A)')  "<--- ", trim(this%data%name_real3(k))

             n           = MAX_BUFFER_SIZE
             buffer_pos  = n + 1

             do g = 1, this%totnat

                call recv_real3

                i  = this%global2local(g)

                if (i > 0 .and. i <= this%natloc) then
                   ! Okay, this particle is on this processor

                   cur_r  = this%data%data_real3(1:3, i, k)

                else
                   ! Particle must be somewhere else
             
                   call find_real3(g, cur_r)

                endif

                write(un, '(1X,3ES20.10)')  cur_r
             enddo

             call mpi_bcast(NEXT_SECTION, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

          endif
       enddo

       call fclose(un)

    else

       !
       ! Send requested information to the root
       !

       call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
       PASS_MPI_ERROR(ierr, error)

       curi  = 1
       do while (p /= NEXT_SECTION)

          if (p == mpi_id()) then

             buffer_size  = 0

             do while (curi <= this%totnat .and. buffer_size < MAX_BUFFER_SIZE)
                j  = this%global2local(curi)
                if (j > 0 .and. j <= this%natloc) then
                   buffer_size  = buffer_size + 1
                   indbuf(buffer_size, 1)       = this%index(j)
                   ibuf(buffer_size, 1)         = this%Z(j)
                   mass(buffer_size, 1)         = this%m(j)
                   xbuf(buffer_size, 1)         = POS(this, j, 1)
                   ybuf(buffer_size, 1)         = POS(this, j, 2)
                   zbuf(buffer_size, 1)         = POS(this, j, 3)
                   group(buffer_size, 1)        = this%g(j)
                   if (associated(this_dissipation)) then
                      dissipation(buffer_size, 1)  = this_dissipation(j)
                   endif
                   if (associated(this_T)) then
                      T(buffer_size, 1)            = this_T(j)
                   endif
                   if(present(mol)) then
                      next(buffer_size, 1)         = mol%next(j)
                   else
                      next(buffer_size, 1)         = 0
                   end if
                endif

                curi  = curi + 1

             enddo

             call mpi_send(indbuf(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

             call mpi_send(ibuf(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)
             call mpi_send(mass(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)
             call mpi_send(xbuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)
             call mpi_send(ybuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)
             call mpi_send(zbuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)
             call mpi_send(group(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)
             if (associated(this_dissipation)) then
                call mpi_send(dissipation(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
                PASS_MPI_ERROR(ierr, error)
             endif
             if (associated(this_T)) then
                call mpi_send(T(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
                PASS_MPI_ERROR(ierr, error)
             endif
             call mpi_send(next(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

          endif

          call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
          PASS_MPI_ERROR(ierr, error)

       enddo

       if (associated(this_v)) then
       
          !
          ! Velocities
          !

          call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
          PASS_MPI_ERROR(ierr, error)

          curi = 1
          do while (p /= NEXT_SECTION)

             if (p == mpi_id()) then
                call send_real3(curi, this_v)
             endif

             call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

          enddo

       endif

       if (associated(this_f)) then

          !
          ! Forces
          !

          call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
          PASS_MPI_ERROR(ierr, error)

          curi = 1
          do while (p /= NEXT_SECTION)

             if (p == mpi_id()) then
                call send_real3(curi, this_f)
             endif

             call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

          enddo

       endif

       !
       ! Real arrays
       !

       do k = 1, this%data%n_real
          if (iand(this%data%tag_real(k), F_RESTART) /= 0) then

             call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

             curi = 1
             do while (p /= NEXT_SECTION)

                if (p == mpi_id()) then
                   call send_real(curi, this%data%data_real(:, k))
                endif

                call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
                PASS_MPI_ERROR(ierr, error)

             enddo

          endif
       enddo

       !
       ! Integer arrays
       !

       do k = 1, this%data%n_integer
          if (iand(this%data%tag_integer(k), F_RESTART) /= 0) then

             call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

             curi = 1
             do while (p /= NEXT_SECTION)

                if (p == mpi_id()) then
                   call send_integer(curi, this%data%data_integer(:, k))
                endif

                call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
                PASS_MPI_ERROR(ierr, error)

             enddo

          endif
       enddo

       !
       ! Real3 arrays
       !

       do k = 1, this%data%n_real3
          if (iand(this%data%tag_real3(k), F_RESTART) /= 0 .and. &
              uppercase(trim(this%data%name_real3(k))) /= "VELOCITIES" .and. &
              uppercase(trim(this%data%name_real3(k))) /= "FORCES") then

             call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
             PASS_MPI_ERROR(ierr, error)

             curi = 1
             do while (p /= NEXT_SECTION)

                if (p == mpi_id()) then
                   call send_real3(curi, this%data%data_real3(:, :, k))
                endif

                call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
                PASS_MPI_ERROR(ierr, error)

             enddo

          endif
       enddo

    endif

    call timer_stop("write_atoms")

  contains

    subroutine recv_real
      implicit none

      integer  :: p, ierr

      ! ---

      do p = 1, mpi_n_procs-1
         if (n(p) == MAX_BUFFER_SIZE .and. buffer_pos(p) > n(p)) then

            call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_recv(indbuf(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_get_count(status, MPI_INTEGER, n(p), ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_recv(xbuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)

            buffer_pos(p)  = 1

         endif
      enddo

    endsubroutine recv_real


    subroutine send_real(curi, r)
      implicit none

      integer, intent(inout)  :: curi
      real(DP), intent(in)    :: r(this%maxnatloc)

      ! ---

      integer  :: buffer_size, j

      ! ---

      buffer_size  = 0

      do while (curi <= this%totnat .and. buffer_size < MAX_BUFFER_SIZE)
         j  = this%global2local(curi)
         if (j > 0 .and. j <= this%natloc) then
            buffer_size  = buffer_size + 1
            indbuf(buffer_size, 1)  = this%index(j)
            xbuf(buffer_size, 1)    = r(j)
         endif

         curi  = curi + 1
         
      enddo

      call mpi_send(indbuf(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)

      call mpi_send(xbuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)

    endsubroutine send_real


    subroutine recv_integer
      implicit none

      integer  :: p, ierr

      ! ---

      do p = 1, mpi_n_procs-1
         if (n(p) == MAX_BUFFER_SIZE .and. buffer_pos(p) > n(p)) then

            call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_recv(indbuf(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_get_count(status, MPI_INTEGER, n(p), ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_recv(ibuf(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)

            buffer_pos(p)  = 1

         endif
      enddo

    endsubroutine recv_integer


    subroutine send_integer(curi, r)
      implicit none

      integer, intent(inout)  :: curi
      integer, intent(in)     :: r(this%maxnatloc)

      ! ---

      integer  :: buffer_size, j

      ! ---

      buffer_size  = 0

      do while (curi <= this%totnat .and. buffer_size < MAX_BUFFER_SIZE)
         j  = this%global2local(curi)
         if (j > 0 .and. j <= this%natloc) then
            buffer_size  = buffer_size + 1
            indbuf(buffer_size, 1)  = this%index(j)
            ibuf(buffer_size, 1)    = r(j)
         endif

         curi  = curi + 1
         
      enddo

      call mpi_send(indbuf(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)

      call mpi_send(ibuf(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)

    endsubroutine send_integer


    subroutine recv_real3
      implicit none

      integer  :: p, ierr

      ! ---

      do p = 1, mpi_n_procs-1
         if (n(p) == MAX_BUFFER_SIZE .and. buffer_pos(p) > n(p)) then

            call mpi_bcast(p, 1, MPI_INTEGER, ROOT, MPI_COMM_WORLD, ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_recv(indbuf(:, p), MAX_BUFFER_SIZE, MPI_INTEGER, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_get_count(status, MPI_INTEGER, n(p), ierr)
            PASS_MPI_ERROR(ierr, error)

            call mpi_recv(xbuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)
            call mpi_recv(ybuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)
            call mpi_recv(zbuf(:, p), MAX_BUFFER_SIZE, MPI_DOUBLE_PRECISION, p, 0, MPI_COMM_WORLD, status, ierr)
            PASS_MPI_ERROR(ierr, error)

            buffer_pos(p)  = 1

         endif
      enddo

    endsubroutine recv_real3


    subroutine send_real3(curi, r)
      implicit none

      integer, intent(inout)  :: curi
      real(DP), intent(in)    :: r(3, this%maxnatloc)

      ! ---

      integer  :: buffer_size, j

      ! ---

      buffer_size  = 0

      do while (curi <= this%totnat .and. buffer_size < MAX_BUFFER_SIZE)
         j  = this%global2local(curi)
         if (j > 0 .and. j <= this%natloc) then
            buffer_size  = buffer_size + 1
            ASSERT(curi == this%index(j), "curi == this%index(j)", error)
            indbuf(buffer_size, 1)  = this%index(j)
            xbuf(buffer_size, 1)    = r(1, j)
            ybuf(buffer_size, 1)    = r(2, j)
            zbuf(buffer_size, 1)    = r(3, j)
         endif

         curi  = curi + 1
         
      enddo

      call mpi_send(indbuf(:, 1), buffer_size, MPI_INTEGER, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)

      call mpi_send(xbuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)
      call mpi_send(ybuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)
      call mpi_send(zbuf(:, 1), buffer_size, MPI_DOUBLE_PRECISION, ROOT, 0, MPI_COMM_WORLD, ierr)
      PASS_MPI_ERROR(ierr, error)

    endsubroutine send_real3


    subroutine find_real(g, x_out, error)
      implicit none

      integer, intent(in)               :: g
      real(DP), intent(out)             :: x_out
      integer, intent(inout), optional  :: error

      ! ---

      logical  :: found

      integer  :: p, j

      ! ---

      found = .false.

      do p = 1, mpi_n_procs-1
         if (n(p) > 0) then
            j  = buffer_pos(p)
         
            if (.not. (j < 1 .or. j > n(p) .or. j > MAX_BUFFER_SIZE)) then

               if (j <= n(p) .and. indbuf(j, p) == g) then
                  if (found) then
                     RAISE_ERROR("Particle was found twice.", error)
                  endif

                  found = .true.

                  x_out  = xbuf(j, p)

                  buffer_pos(p)  = j + 1
               endif
            endif

         endif
      enddo

      if (.not. found) then
         RAISE_ERROR("Particle not found.", error)
      endif
      
    endsubroutine find_real


    subroutine find_integer(g, i_out, error)
      implicit none

      integer, intent(in)               :: g
      integer, intent(out)              :: i_out
      integer, intent(inout), optional  :: error

      ! ---

      logical  :: found

      integer  :: p, j

      ! ---

      found = .false.

      do p = 1, mpi_n_procs-1
         if (n(p) > 0) then
            j  = buffer_pos(p)

            if (.not. (j < 1 .or. j > n(p) .or. j > MAX_BUFFER_SIZE)) then

               if (j <= n(p) .and. indbuf(j, p) == g) then
                  if (found) then
                     RAISE_ERROR("Particle was found twice.", error)
                  endif

                  found = .true.

                  i_out  = ibuf(j, p)

                  buffer_pos(p)  = j + 1
               endif
            endif

         endif
      enddo

      if (.not. found) then
         RAISE_ERROR("Particle not found.", error)
      endif
      
    endsubroutine find_integer


    subroutine find_real3(g, r_out, error)
      implicit none

      integer, intent(in)               :: g
      real(DP), intent(out)             :: r_out(3)
      integer, intent(inout), optional  :: error

      ! ---

      logical  :: found

      integer  :: p, j

      ! ---

      found = .false.

      do p = 1, mpi_n_procs-1
         if (n(p) > 0) then
            j  = buffer_pos(p)

            if (.not. (j < 1 .or. j > n(p) .or. j > MAX_BUFFER_SIZE)) then

               if (j <= n(p) .and. indbuf(j, p) == g) then
                  if (found) then
                     RAISE_ERROR("Particle was found twice.", error)
                  endif

                  found  = .true.

                  r_out  = (/ xbuf(j, p), ybuf(j, p), zbuf(j, p) /)
                  
                  buffer_pos(p)  = j + 1
               endif
            endif

         endif
      enddo

      if (.not. found) then
         RAISE_ERROR("Particle not found.", error)
      endif
      
    endsubroutine find_real3
    
  endsubroutine internal_native_io_write_atoms

#else

  !**********************************************************************
  ! Write the particles to an atoms.out file (serial version)
  !**********************************************************************    
  subroutine native_io_write_atoms(p, fn, mol, error)
    implicit none

    type(particles_t), intent(in)     :: p
    character(*), intent(in)          :: fn
    type(molecules_t), intent(in), optional :: mol
    integer, intent(inout), optional  :: error

    ! ---

    integer            :: un, i, j, k, l, g
    real(DP)           :: cur_T, cur_diss

    real(DP), pointer  :: v(:, :)
    real(DP), pointer  :: f(:, :)
    real(DP), pointer  :: T(:)
    real(DP), pointer  :: dissipation(:)

    ! ---

    call timer_start("write_atoms")

    v            => NULL()
    f            => NULL()
    T            => NULL()
    dissipation  => NULL()
    if (exists(p%data, V_STR)) then
       call ptr_by_name(p%data, V_STR, v)
    endif
    if (exists(p%data, F_STR)) then
       call ptr_by_name(p%data, F_STR, f)
    endif
    if (exists(p%data, T_STR)) then
       call ptr_by_name(p%data, T_STR, T)
    endif
    if (exists(p%data, DISSIPATION_STR)) then
       call ptr_by_name(p%data, DISSIPATION_STR, dissipation)
    endif

    un = fopen(fn, F_WRITE, error=error)
    PASS_ERROR(error)

    write (un, '(A)')  "<--- Total number of atoms"

    write (un, '(I20)')  p%natloc

    write (un, '(A)')  "<--- *** The following line is ignored ***"
    write (un, *)
    write (un, '(A)')  "<--- Element, atomic mass, coordinates, group, dissipation, temperature, (next)"

    do g = 1, p%natloc
       i  = p%global2local(g)

       if (p%Z(i) > 0 .and. p%Z(i) <= MAX_Z) then
          cur_T  = 0.0_DP
          if (associated(T)) then
             cur_T  = T(i)
          endif
          cur_diss  = 0.0_DP
          if (associated(dissipation)) then
             cur_diss  = dissipation(i)
          endif

          if(present(mol)) then
             if (mol%next(i) > 0) then
                write(un, '(1X,A4,4ES20.10,I5,2ES20.10,I10)') &
                     ElementName(p%Z(i)), ElementMass(p%Z(i)), &
                     POS(p, i, 1), POS(p, i, 2), POS(p, i, 3), &
                     p%g(i), cur_diss, cur_T, mol%next(i)
             else
                write(un, '(1X,A4,4ES20.10,I5,2ES20.10)') &
                     ElementName(p%Z(i)), ElementMass(p%Z(i)), &
                     POS(p, i, 1), POS(p, i, 2), POS(p, i, 3), &
                     p%g(i), cur_diss, cur_T
             endif
          else
             write(un, '(1X,A4,4ES20.10,I5,2ES20.10)') &
                  ElementName(p%Z(i)), ElementMass(p%Z(i)), &
                  POS(p, i, 1), POS(p, i, 2), POS(p, i, 3), &
                  p%g(i), cur_diss, cur_T
          end if
       else
          RAISE_ERROR("Unknown atomic number encountered.", error)
       endif
    enddo

    if (associated(v)) then

       write (un, '(A)')  "<--- Velocities"
       
       do g = 1, p%natloc
          i  = p%global2local(g)
          write(un, '(1X,3ES20.10)')  VEC(v, i, 1), VEC(v, i, 2), VEC(v, i, 3)
       enddo

    endif

    write (un, '(A)')  "<--- Forces"

    do g = 1, p%natloc
       i  = p%global2local(g)
       if (associated(f)) then
          write(un, '(1X,3ES20.10)')  VEC(f, i, 1), VEC(f, i, 2), VEC(f, i, 3)
       else
          write(un, '(1X,3ES20.10)')  0.0_DP, 0.0_DP, 0.0_DP
       endif
    enddo

    !
    ! Everything else ... dynamic!
    !

    do i = 1, p%data%n_real_attr
       write (un, '(1X,A,A)')      "<--- ", trim(p%data%name_real_attr(i))
       write (un, '(1X,ES20.10)')  p%data%data_real_attr(i)
    enddo

    do i = 1, p%data%n_real3_attr
       write (un, '(1X,A,A)')       "<--- ", trim(p%data%name_real3_attr(i))
       write (un, '(1X,3ES20.10)')  p%data%data_real3_attr(:, i)
    enddo

    do i = 1, p%data%n_real3x3_attr
       write (un, '(1X,A,A)')       "<--- ", trim(p%data%name_real3x3_attr(i))
       write (un, '(1X,3ES20.10)')  ( p%data%data_real3x3_attr(:, j, i), j = 1, 3 )
    enddo

    do i = 1, p%data%n_real
       if (iand(p%data%tag_real(i), F_RESTART) /= 0) then
          write (un, '(1X,A,A)')  "<--- ", trim(p%data%name_real(i))
          write (un, '(1X,ES20.10)')  ( p%data%data_real(p%global2local(j), i), j = 1, p%natloc )
       endif
    enddo

    do i = 1, p%data%n_integer
       if (iand(p%data%tag_integer(i), F_RESTART) /= 0) then
          write (un, '(1X,A,A)')  "<--- ", trim(p%data%name_integer(i))
          write (un, '(1X,I10)')  ( p%data%data_integer(p%global2local(j), i), j = 1, p%natloc )
       endif
    enddo

    do i = 1, p%data%n_real3
       if (iand(p%data%tag_real3(i), F_RESTART) /= 0 .and. &
           uppercase(trim(p%data%name_real3(i))) /= "VELOCITIES" .and. &
           uppercase(trim(p%data%name_real3(i))) /= "FORCES") then
          write (un, '(1X,A,A)')  "<--- ", trim(p%data%name_real3(i))
          write (un, '(1X,3ES20.10)') &
               ( ( &
               p%data%data_real3(k, p%global2local(j), i), &
               k = 1, 3 ), &
               j = 1, p%natloc )
       endif
    enddo

    do i = 1, p%data%n_real3x3
       if (iand(p%data%tag_real3x3(i), F_RESTART) /= 0) then
          write (un, '(1X,A,A)')  "<--- ", trim(p%data%name_real3x3(i))
          write (un, '(1X,3ES20.10)') &
               ( ( ( &
               p%data%data_real3x3(k, l, p%global2local(j), i), &
               k = 1, 3 ), &
               l = 1, 3 ), &
               j = 1, p%natloc )
       endif
    enddo

    call fclose(un)

    call timer_stop("write_atoms")

  endsubroutine native_io_write_atoms

#endif

endmodule native_io
