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
#include "filter.inc"

!>
!! Molecular information
!!
!! Contains information on the (for example molecular) composition of
!! the system. It contains the Verlet-type lists "next" and "head".
!! For example for a system composed of two water molecules, the lists
!! would look like this:
!!
!!        O  H  O  H  H  H
!!
!! atom:  1  2  3  4  5  6
!!
!! next:  4  6  2  5  0  0
!!
!! head:  1  3
!!
!! imol:  1  2  2  1  1  2
!!
!! Note: The head-list index is the molecule index given in imol.
!!
!! Note that the next list is the main data structure. The other lists can be
!! updated based on it, using molecules_update_head. One should take care that
!! the head list is updated when it's needed, as it's not done automatically.
!! On the other hand, all routines here should also update the head list on the
!! fly, so it should be enough to call molecules_update_head just
!! once, unless you modify the arrays directly.
!!
!! Note: Even if the routines here keep the head list up to date,
!! calling molecules_update_head will probably result in a different
!! (equally valid) head list!
!!
!! For atoms that are for some reason not part of the molecule structure, the next
!! array should be -1 and imol 0.
!!
!! To loop over atoms in molecule i, do:
!!
!!
!! iat = head(i)
!!
!! do while(iat > 0)
!!
!!   ...
!!
!!   iat = next(iat)
!!
!! end do
!!
!!
!!
!! NOTE: Because of the above, a good coding convention would be to only handle
!!       other arrays than next in this module, so that possible changes in the
!!       internals will not affect things outside, and so that the other arrays
!!       get automatically updated. This is already broken in molecular
!!       decomposed TB, in the molecule separation routine.
!!
!! TODO: Add routines to construct the molecule information based on clustering
!! and other criteria.
!!
!! (Tommi Jaervi, 2009)
!<
module molecules
  use libAtoms_module

  use logging

  use misc
  use timer

  use data
  use particles
  use filter

  implicit none

  private

  public :: init
  interface init
     module procedure molecules_init
  end interface

  public :: del
  interface del
     module procedure molecules_del
  end interface

  public :: register_data
  interface register_data
     module procedure molecules_register_data
  end interface

  integer :: molecules_str_i = 0                              !< (Internal use: For generating independent Id strings for pointers next, imol => particles%data)

  public :: molecules_t
  type molecules_t

     integer                            :: n_atoms = 0        !< Number of atoms in molecules (only including atoms actually referenced in molecules => may be smaller than the length of next etc. arrays!)
     integer                            :: n_molec = 0        !< Number of molecules
     integer, dimension(:), pointer     :: next => NULL()     !< next(i), i in [1,#atoms], is the next atom in the molecule from atom i (should use global indices)
     character(MAX_NAME_STR)            :: next_str           !< Id for pointer: next => particles%data
     integer, dimension(:), allocatable :: head               !< head(i), i in [1,n_molec] is the first atom in the next-array sequence of molecule i

     logical                            :: use_imol = .false. !< imol not allocated by default
     integer, dimension(:), pointer     :: imol => NULL()     !< Molecule atom belongs to
     character(MAX_NAME_STR)            :: imol_str           !< Id for pointer: imol => particles%data

  end type molecules_t

  ! Tommi, should these guys become interfaces?
  public :: molecules_assign_one, molecules_count_atoms, molecules_join_list, molecules_output_arrays
  public :: molecules_update_head, molecules_verify, molecules_copy, molecules_swap
  public :: molecules_separate

contains


  subroutine molecules_copy(this, from, ierror)
    implicit none

    type(molecules_t), intent(inout)   :: this    !< This molecules object
    type(molecules_t), intent(in)      :: from    !< Where to copy data from
    integer, intent(inout), optional   :: ierror  !< Error signals

    ! ---

    integer        :: i

    ! ---

    ! - checks

    if(.not. associated(this%imol) .or. .not. associated(this%next) .or. .not. associated(from%imol) .or. .not. associated(from%next)) then
       write(ilog,*), "pointers and strings, this:", this%next_str, this%next, this%imol_str, this%imol
       write(ilog,*), "pointers and strings, from:", from%next_str, from%next, from%imol_str, from%imol
       RAISE_ERROR("molecules_copy(): One of source or destination arrays in particles object not allocated.", ierror)
    end if

    ! - copy

    ! variables
    !   note: do not copy *_str, because then the storage location of our data would change
    this%n_atoms = from%n_atoms
    this%n_molec = from%n_molec
    this%use_imol = from%use_imol

    ! arrays
    if(allocated(this%head)) deallocate(this%head)
    if(size(from%head) > 0) then
       allocate(this%head(size(from%head)))
    end if

    do i = 1, size(this%next)
       this%next(i) = from%next(i)
       if(this%use_imol) this%imol(i) = from%imol(i)
    end do
    do i = 1, size(this%head)
       this%head(i) = from%head(i)
    end do

  end subroutine molecules_copy


  !>
  !! Register per atom data
  !!
  !! Register per atom data
  !<
  subroutine molecules_register_data(this, p, use_imol)
    implicit none

    type(molecules_t), intent(inout) :: this      !< The molecules object
    type(particles_t), intent(inout) :: p         !< Particles object
    logical, intent(in), optional    :: use_imol  !< Allocate imol-array?

    ! - set use_imol
    if(present(use_imol)) then
       this%use_imol = use_imol
    end if

    ! - arrays
    molecules_str_i = molecules_str_i + 1

    ! next
    this%next_str = "molecules_next_" // molecules_str_i // ""
    call add_integer( &
         p%data, &
         this%next_str, &
         F_CONSTANT + F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS )

    ! imol
    if (present(use_imol)) then
       if (use_imol) then
          this%imol_str = "molecules_imol_" // molecules_str_i // ""
          call add_integer( &
               p%data, &
               this%imol_str, &
               F_CONSTANT + F_VERBOSE_ONLY + F_COMMUNICATE + F_COMM_GHOSTS )
       end if
    endif

  end subroutine molecules_register_data


  !>
  !! Constructor
  !!
  !! Initialize the molecules object and register the next array in particles%data.
  !<
  subroutine molecules_init(this, p, use_imol)
    implicit none

    type(molecules_t), intent(inout)    :: this      !< The molecules object
    type(particles_t), intent(inout)    :: p         !< Particles object
    logical, intent(in), optional    :: use_imol  !< Allocate imol-array

    call molecules_register_data(this, p, use_imol)

  end subroutine molecules_init


  !>
  !! Destructor
  !!
  !! Destroy the molecules object
  !<
  subroutine molecules_del(this)
    implicit none

    type(molecules_t), intent(inout) :: this  !< The molecules object

    ! head
    if(allocated(this%head)) then
       deallocate(this%head)
    end if

    ! others
    ! Tommi: XXX
    this%next => NULL()
    this%use_imol = .false.
    this%imol => NULL()

  end subroutine molecules_del


  !>
  !! Verify pointers
  !!
  !! Verify pointers. Should be called efore trying to access the object data.
  !! (Called internally but also from, for example,
  !! native_io_read_atoms to verify the object.)
  !<
  subroutine molecules_verify(this, p, ierror)
    implicit none

    type(molecules_t), intent(inout) :: this    !< This molecules object
    type(particles_t), target        :: p       !< Particles object
    integer, intent(inout), optional :: ierror  !< Error passing

    ! ---

    ! next
    if(.not. associated(this%next)) then
       call ptr_by_name(p%data, this%next_str, this%next, ierror=ierror)
       PASS_ERROR(ierror)
    end if

    ! imol
    if(this%use_imol .and. .not. associated(this%imol)) then
       call ptr_by_name(p%data, this%imol_str, this%imol, ierror=ierror)
       PASS_ERROR(ierror)
    end if

  end subroutine molecules_verify


  !>
  !! Update molecule head, imol, n_atoms
  !!
  !! Update molecule head list and n_molec from next list only. If the allocated
  !! array is too small, it is automatically reallocated. Also imol is updated
  !! if it's used.
  !!
  !! Note: head array produced by this routine should contain atom
  !!       indices in increasing order.
  !<
  subroutine molecules_update_head(this, p, ierror)
    implicit none

    type(molecules_t), intent(inout)  :: this    !< The molecules object
    type(particles_t), target         :: p       !< The particles object (only number of atoms used)
    integer, intent(inout), optional  :: ierror  !< Error signals

    ! ---

    integer                   :: i         ! Loops
    logical, dimension(p%nat) :: startmol  ! Does this atom start a molecule (i.e., next-array sequence)?
    integer                   :: sizeadd = 1000  ! How much size to add to the head array with one command, when it gets too small

    ! check
    integer :: j, k
    logical   :: done(p%nat)     ! has this atom been assigned already? (just for checks, could be removed too)

    ! ---

    call timer_start('molecules_update_head')

    ! Make sure pointers initialized and reset
    call molecules_verify(this, p, ierror)
    PASS_ERROR(ierror)
    this%head = 0

    ! Check which atoms start a new molecule
    startmol = .true.
    do i = 1, p%nat
       ! the next atom from this is not starting a molecule
       if(this%next(i) > 0) then
          startmol(this%next(i)) = .false.
       endif
       ! atom not included in molecular structure
       if(this%next(i) == -1) then
          startmol(i) = .false.
       end if

       if(this%next(i) == i) then
          RAISE_ERROR("molecules_update_head: Next atom in molecule is atom itself for atom " // i // ".", ierror)
       end if
    enddo

    ! Construct head-list and count atoms in molecule structure
    this%n_molec = 0
    this%n_atoms = 0
    do i = 1, p%nat
       if(this%next(i) /= -1) this%n_atoms = this%n_atoms + 1
       if(startmol(i)) then
          this%n_molec = this%n_molec + 1
          if(size(this%head) < this%n_molec) then
             !write(ilog,*), "molecules_update_head: Resizing head list to", this%n_molec + sizeadd
             call resize(this%head, this%n_molec + sizeadd)
          endif
          this%head(this%n_molec) = i
       end if
    end do

    ! Update imol list
    this%imol = 0  ! 0 for atoms not in any molecule (next=-1)
    done = .false.
    if(this%use_imol) then
       do i = 1, this%n_molec
          ! loop over atoms in this molecule
          j = this%head(i)
          do while(j > 0)
             ! check
             if (done(j)) then
                write(ilog, *), "molecules_update_head: Same atom belongs to more than one molecule."
                write(ilog, *), "Current molec, first atom in molec, atom: ", i, this%head(i), j
                write(ilog, *), "Next, done and startmol arrays, and atom Z:"
                do k = 1, p%nat
                   write(ilog, *), k, this%next(k), done(k), startmol(k), p%Z(k)
                enddo
                RAISE_ERROR("molecules_update_head: Same atom belongs to more than one molecule.", ierror)
             endif

             this%imol(j) = i
             done(j) = .true.
             j = this%next(j)
          end do  ! end of loop over atoms in this molecule
       end do  ! end of loop over molecules
    end if

    call timer_stop('molecules_update_head')

  end subroutine molecules_update_head


  !>
  !! Assign all atoms to one molecule
  !!
  !! Assign all atoms to one molecule and updates the head list accordingly.
  !<
  subroutine molecules_assign_one(this, p, f, ierror)
    implicit none

    type(molecules_t), intent(inout)   :: this    !< This molecules object
    type(particles_t), target          :: p       !< Particles object (only used for number of atoms)
    integer, intent(in), optional      :: f       !< Filter for atoms to assign to molecule
    integer, intent(inout), optional   :: ierror  !< Error passing

    ! ---

    integer :: i      ! loops
    logical :: first  ! first atom to be added to molecule?
    integer :: prev   ! previous atom added to molecule
    integer :: filter ! filter used
    integer :: nf     ! number of local atoms matching filter

    ! ---

    ! make sure pointers initialized
    call molecules_verify(this, p, ierror)
    PASS_ERROR(ierror)

    ! init filter
    if(present(f)) then
       filter = f
    else
       filter = filter_from_string("*", p, ierror=ierror)
       PASS_ERROR(ierror)
    end if
    nf = filter_count(filter, p)

    ! head array
    this%n_molec = 1
    if(.not. allocated(this%head) .or. size(this%head) < this%n_molec) then
       call resize(this%head, 10)
    endif

    ! add atoms to next array, construct head and imol
    this%n_atoms = 0
    first = .true.
    prev = -1
    do i = 1, p%nat - 1
       if(IS_EL(filter, p, i)) then
          this%n_atoms = this%n_atoms + 1
          ! first atom, set head
          if(first) then
             this%head(1) = i
             first = .false.
             prev = i
          else
             this%next(prev) = i
             if(this%use_imol) this%imol(prev) = 1
             prev = i
          end if
       else
          this%next(i) = -1
          if(this%use_imol) this%imol(i) = 0
       end if
    end do
    ! last atom
    if(IS_EL(filter, p, p%nat)) then
       this%n_atoms = this%n_atoms + 1
       if (prev /= -1) then
          this%next(prev) = p%nat
          if(this%use_imol) this%imol(prev) = 1
       endif
       this%next(p%nat) = 0
       if(this%use_imol) this%imol(p%nat) = 1
    else
       if(prev /= -1) then
          this%next(prev) = 0
          if(this%use_imol) this%imol(prev) = 1
       end if
       this%next(p%nat) = -1
       if(this%use_imol) this%imol(p%nat) = 0
    end if

    ! checks
    if(this%n_atoms==0) then
       RAISE_ERROR("molecules_assign_one: No atoms seem to belong to molecules, according to filter: " // filter, ierror)
    end if
    if(this%n_atoms /= nf) then
       RAISE_ERROR("molecules_assign_one: Bug: Number of atoms in constructed structure doesn't match that from filter.", ierror)
    end if

  end subroutine molecules_assign_one


  !>
  !! Separate molecules
  !!
  !! Separate molecule according to input list of the following format:
  !!
  !! list: ns a11 a12 ... 0 a21 a22 ...0
  !!
  !! For example, separating a water dimer to two molecules (atoms 1,2,3 and 4,5,6)
  !! would look like this:
  !!
  !! list: 2 | 2 3 1 0 | 4 6 5 0
  !!
  !! The last new molecule will replace the old one (in molecule id, which is the head-list
  !! index or imol entry). The subsequent ones will be in order at the end of the
  !! head array.
  !!
  !! Note: If the imol array is used (use_imol=.true.), a check is made if an atom in the list
  !! actually belongs to the molecule to be separated. If imol is not used, no checks are made!
  !!
  !! Note: Since the new molecules are put where the separated one
  !! was, and at the end of the head array, it's safe to have a loop
  !! of the type
  !!
  !!    do i = 1, nmol    ! <- not mol%n_molec because it keeps changing!
  !!
  !!      call molecules_separate(mol, i, list(i,:), ierror)
  !!
  !!    end do
  !!
  !! but note that mol%n_molec cannot be used in the loop since it
  !! changes!
  !<
  subroutine molecules_separate(this, imol, list, ierror)
    implicit none

    type(molecules_t), intent(inout)   :: this    !< This molecules object
    integer, intent(in)                :: imol    !< Molecule to separate
    integer, dimension(:), intent(in)  :: list    !< List to separate after
    integer, intent(inout), optional   :: ierror  !< Error signals

    ! ---

    integer           :: head    ! new head array position
    integer           :: nread   ! number of molecules read so far
    logical           :: start   ! starting a new molecule?
    integer           :: pos     ! current position in list
    integer           :: at      ! current atom
    integer           :: nmolec  ! number of atoms in current molecule
    integer           :: prev    ! previous atom in molecule

    ! ---

    call timer_start('molecules_separate')

    ! loop over new molecules
    nread = 0
    start = .true.
    pos = 2
    do while(nread < list(1))

       ! set head array position for new molecule
       if(start) then
          start = .false.
          nmolec = 0
          if(nread==list(1)-1) then
             head = imol
          else if(nread==0) then
             head = this%n_molec + 1
          else
             head = head + 1
          end if
       end if

       ! make sure head is large enough
       if(size(this%head) < head) then
          call resize(this%head, head + 100)  ! XXX: hard-coding
       end if

       ! check
       if(this%use_imol .and. list(pos) /= 0) then
          if(this%imol(list(pos)) /= imol) then
             RAISE_ERROR("molecules_separate: Trying to separate an atom from a molecule where it's not present.", ierror)
          end if
       end if

       ! check if we're done for this molecule
       if(list(pos)==0) then
          start = .true.
          this%next(prev) = 0
          nread = nread + 1
       end if

       ! add atom
       if(list(pos)/=0) then
          at = list(pos)
          nmolec = nmolec + 1
          if(nmolec==1) then
             this%head(head) = at
          else
             this%next(prev) = at
          end if
          prev = at

          ! imol array
          if(this%use_imol) then
             this%imol(at) = head
          end if
       end if

       pos = pos + 1
    end do  ! end of loop over molecules

    ! update number of molecules
    this%n_molec = this%n_molec + nread - 1

    call timer_stop('molecules_separate')

  end subroutine molecules_separate


  !>
  !! Join molecules based on list
  !!
  !! Join molecules based on list. Use this to safely join many
  !! molecule pairs at once. Each molecule can be joined with
  !! more than one other. The joined molecule will appear with the
  !! index which is the lower of the two molecules joined. (This is
  !! important for the functioning of some parts of mdcore.)
  !<
  subroutine molecules_join_list(this, nmol, list1, list2, ierror)

    type(molecules_t), intent(inout)              :: this          !< This molecules object
    integer, intent(in)                           :: nmol          !< Number of molecule pairs to be joined
    integer, intent(inout), dimension(:)          :: list1, list2  !< Join molecules pairs list1(i), list2(i)
    integer, optional, intent(inout)              :: ierror        !< Error signals

    ! ---

    integer           :: i, j           ! loops
    integer           :: i1, i2         ! current molecules to be joined

    ! ---

    ! loop over pairs to be joined
    do i = 1, nmol
       i1 = min(list1(i), list2(i))
       i2 = max(list1(i), list2(i))

       if(i1==i2) then
          RAISE_ERROR("molecules_join_list: Trying to join a molecule with itself.", ierror)
       end if
       if(i1 > this%n_molec .or. i2 > this%n_molec) then
          write(ilog,*), ""
          write(ilog,*), "Molecule pairs listed for joining:"
          do j = 1, nmol
             write(ilog,'(3I10)'), j, list1(j), list2(j)
          end do
          RAISE_ERROR("molecules_join_list: Trying to join molecules " // i1 // " and " // i2 // " but object only contains " // this%n_molec // " molecules.", ierror)
       end if

       call molecules_join(this, i1, i2, ierror)
       PASS_ERROR(ierror)

       ! correct lists
       !   new molecule is in i1 and the last one has been moved to i2
       !   -> i1 stays at the same index -> no changes needed
       !   -> change all references of i2 to i1
       !   -> in lists, change index of last molecule to i2
       !      (==n_molec+1 because molecules_join call above reduces
       !       n_molec by 1)
       do j = i+1, nmol
          if(list1(j)==i2) list1(j) = i1
          if(list2(j)==i2) list2(j) = i1
          if(i2 /= this%n_molec+1) then
             if(list1(j)==this%n_molec+1) list1(j) = i2
             if(list2(j)==this%n_molec+1) list2(j) = i2
          end if
       end do
    end do

  end subroutine molecules_join_list


  !>
  !! Join molecules
  !!
  !! Join two molecules so that new molecule is at min(n1,n2) and the
  !! last molecule will be moved to n2.
  !! (That's important to keep consistent with joining Hamiltonian
  !! blocks!)
  !<
  subroutine molecules_join(this, n1, n2, ierror)
    implicit none

    type(molecules_t), intent(inout)   :: this    !< This molecules object
    integer, intent(in)                :: n1, n2  !< Molecules to be joined (referring to head-list indices)
    integer, intent(inout), optional   :: ierror  !< Error signals

    ! ---

    integer                            :: k           ! loops
    integer                            :: natmol      ! Number of atoms in current molecule
    integer, dimension(:), allocatable :: atmol       ! Atoms in current molecule (indices) (XXX: hard-coded)
    integer                            :: m1, m2      ! Just n1 and n2 but without intent(in) so can be swapped

    ! ---

    call timer_start('molecules_join')

    m1 = n1
    m2 = n2

    ! - collect atoms in new joined molecule

    natmol = 0
    allocate(atmol(1000))

    ! collect atoms from molecule m1
    k  = this%head(m1)
    do while (k > 0)
       natmol = natmol + 1
       if(size(atmol) < natmol) then
          call resize(atmol, size(atmol) + 1000)  ! XXX: hard-coding
       end if
       atmol(natmol) = k
       k = this%next(k)
    enddo

    ! collect atoms from molecule m2
    k  = this%head(m2)
    do while (k > 0)
       natmol = natmol + 1
       if(size(atmol) < natmol) then
          call resize(atmol, size(atmol) + 1000)  ! XXX: hard-coding
       end if
       atmol(natmol) = k
       k = this%next(k)
    enddo

    ! - atoms collected from two molecules, join

    ! check
    if(natmol == 0) then
       RAISE_ERROR("molecules_join: The joined molecules contain zero atoms in total!", ierror)
    end if

    ! order so that m1 < m2: new molecule will be at m1
    if(m1 > m2) then
       call swap(m1, m2)
    endif

    ! head array
    if(m2==this%n_molec) then
       this%n_molec = this%n_molec - 1
    else
       this%head(m2) = this%head(this%n_molec)
       this%n_molec = this%n_molec - 1

       ! imol array of molecule moved to m2
       k = this%head(m2)
       do while(k>0)
          this%imol(k) = m2
          k = this%next(k)
       end do
    endif
    this%head(m1) = atmol(1)

    ! next array
    do k = 1, natmol-1
       this%next(atmol(k)) = atmol(k+1)
       if(this%use_imol) this%imol(atmol(k)) = m1
    enddo
    this%next(atmol(natmol)) = 0
    if(this%use_imol) this%imol(atmol(natmol)) = m1

    call timer_stop('molecules_join')

  end subroutine molecules_join


  !>
  !! Move all rigid object such that they are not wrapped by periodic
  !! boundaries.
  !!
  !! Move all rigid object such that they are not wrapped by periodic
  !! boundaries.
  !!
  !! XXX: This will fail if the atom starting a molecule (the next-array sequence)
  !!      does not come first in the next-array.
  !<
  subroutine molecules_group_rigid_objects(this, p, ierror)
    implicit none

    type(molecules_t), intent(inout)   :: this    !< Molecules object
    type(particles_t), target          :: p       !< Particles object
    integer, intent(inout), optional   :: ierror  !< Error signals

    ! ---

    real(DP)  :: thres_sq, ref(3), d(3), x(3)

    integer   :: i, j

    logical   :: done(p%natloc)

    ! ---

    ! Make sure pointers initialized
    call molecules_verify(this, p)

    thres_sq  = minval( (/ p%Abox(1, 1)/2, p%Abox(2, 2)/2, p%Abox(3, 3)/2 /) )**2

    done      = .false.

    ! ---

    do i = 1, p%natloc
       if (this%next(i) > 0 .and. .not. done(i)) then
          done(i)  = .true.

          ref      = POS3(p, i)

          j = this%next(i)
          do while (j > 0)
             j = p%global2local(j)

             ! check
             if (done(j)) then
                RAISE_ERROR("Recursive rigid object found. Please check your input files.", ierror)
             endif
             done(j)  = .true.

             d        = POS3(p, j) - POS3(p, i)
             if (dot_product(d, d) > thres_sq) then
                x              = matmul(p%Bbox, d)
                d              = matmul(p%Abox, x - nint(x))
                PNC3(p, j)  = PNC3(p, i) + d
#ifndef IMPLICIT_R
                POS3(p, j)  = POS3(p, i) + d
#endif
             endif

             j  = this%next(j)
          enddo
       endif
    enddo

  endsubroutine molecules_group_rigid_objects


  !>
  !! Output info (for debugging)
  !!
  !! Output next and head arrays for debugging.
  !<
  subroutine molecules_output_arrays(this, unit)
    implicit none

    type(molecules_t), intent(in) :: this  !< This molecules object
    integer                       :: unit  !< Unit to output to

    ! ---

    integer  :: i  ! loops

    ! ---

    write(unit, *), ""
    write(unit, *), "Atoms in molecules structure:", this%n_atoms
    write(unit, *), "Number of molecules:         ", this%n_molec
    write(unit, *), "Head array and sizes:"
    !do i = 1, size(this%head)
    do i = 1, this%n_molec
       write(unit, *), i, this%head(i), molecules_count_atoms(this, i)
    enddo
    write(unit, *), ""
    write(unit, *), "Next and imol arrays:"
    do i = 1, size(this%next)
       write(unit, *), i, this%next(i), this%imol(i)
    enddo
    write(unit, *), ""

  end subroutine molecules_output_arrays


  !>
  !! Atom count
  !!
  !! Count atoms in molecule
  !<
  function molecules_count_atoms(this, i, ierror) result(count)
    implicit none

    type(molecules_t), intent(in)    :: this      !< The molecules object
    integer, intent(in)              :: i         !< Molecule to be counted
    integer, intent(inout), optional :: ierror    !< Error signals

    integer                          :: count     !< Number of atoms in molecule i

    ! ---

    integer        :: k

    ! ---

    ! count
    count = 0
    k = this%head(i)
    do while (k > 0)
       count = count + 1
       k = this%next(k)
    enddo
  end function molecules_count_atoms


  !>
  !! Swap molecules
  !!
  !! Swap two molecules. Swapping molecules i and j swaps head(i) <->
  !! head(j) and for each imol(k)==i -> imol(k)==j and similarly for
  !! j.
  !<
  subroutine molecules_swap(this, i, j, ierror)
    implicit none

    type(molecules_t), intent(inout) :: this      !< The molecules object
    integer, intent(in)              :: i, j      !< Molecules to be swapped
    integer, intent(inout), optional :: ierror    !< Error signals

    ! ---

    integer                          :: k         ! loops, temp

    ! ---

    if(i /= j) then
       ! set imol
       if(this%use_imol) then
          ! change i to j in imol
          k = this%head(i)
          do while(k > 0)
             this%imol(k) = j
             k = this%next(k)
          end do
          ! change j to i in imol
          k = this%head(j)
          do while(k > 0)
             this%imol(k) = i
             k = this%next(k)
          end do
       end if

       ! swap head
       k = this%head(i)
       this%head(i) = this%head(j)
       this%head(j) = k
    end if

  end subroutine molecules_swap


end module molecules
