#define IMPLICIT_R

!>
!! Return a single vector from a vector array, i.e. an array
!! with shape (3, :). Here, *i* is the index of the entry. 
!!
!! VEC3 return the full 3-vector, VEC return a single entry from 
!! the 3-vector (entry *d*).
!<

#define VEC3(v, i)      v(1:3, i)
#define VEC(v, i, d)    v(d, i)


!>
!! Return the position from a particles object. Here, *i*
!! is the index of the particle. Same as VEC3 and VEC just
!! specifically for the particle positions.
!<

#define POS3(p, i)      in_cell(p, p%r_non_cyc(1:3, i))
#define POS(p, i, d)    in_cellc(p, p%r_non_cyc(1:3, i), d)

#define PNC3(p, i)      p%r_non_cyc(1:3, i)
#define PNC(p, i, d)    p%r_non_cyc(d, i)

#if defined(LAMMPS) || defined(MDCORE_PYTHON)
#define PCN3(p, i)      p%r_non_cyc(1:3, i)
#define PCN(p, i, d)    p%r_non_cyc(d, i)
#else
#define PCN3(p, i)      p%r_cont(1:3, i)
#define PCN(p, i, d)    p%r_cont(d, i)
#endif


!>
!! Determine whether the pair given by *i* and *j* is actually a pair,
!! i.e. exclude double summation of pairs.
!!
!! If the system is not periodic this should become *i < j*.
!<

#ifdef LAMMPS

#define IS_PAIR(nl, i, nj, j)  (i < j)

#else

#define IS_PAIR(nl, i, ni, j)  ((DC(VEC(nl%dc, ni, 1), VEC(nl%dc, ni, 2), VEC(nl%dc, ni, 3)) == 0 .and. i < j) .or. DC(VEC(nl%dc, ni, 1), VEC(nl%dc, ni, 2), VEC(nl%dc, ni, 3)) > 0)

#endif


!>
!! Macros computing particle distances.
!!
!! Computation of the vector distance between particle *i* and particle *j*.
!! *ni* is the neighbor list index corresponding to neighbor *j*.
!<

#ifdef LAMMPS

! We need to add an offset of 1 to the atom index. LAMMPS starts counting at 0.
#define GET_NEIGHBOR(nl, ni)  (nl%neighbors(ni)+1)

! LAMMPS does not have the dc array. Particles are always kept as ghosts.
#define GET_DRJ(p, nl, i, j, ni)  ( PNC3(p, i) - PNC3(p, j) )

#else

! Just return the neighbor index
#define GET_NEIGHBOR(nl, ni)  nl%neighbors(ni)

#define GET_DRJ(p, nl, i, j, ni)  ( PNC3(p, i) - PNC3(p, j) + matmul(p%Abox, VEC3(nl%dc, ni)) + p%shear_dx*VEC(nl%dc, ni, 3) )

#endif


!>
!! Make a null-terminated string
!<
#define CSTR(x)  x // C_NULL_CHAR


!>
!! Helper macros for distance computation. They all depend of GET_DRJ which should be
!! adjusted to the specific ementation.
!<

#define GET_ABS_DRJ(p, nl, i, j, ni)  sqrt(dot_product(GET_DRJ(p, nl, i, j, ni), GET_DRJ(p, nl, i, j, ni)))

#define GET_DR(p, nl, i, ni)  GET_DRJ(p, nl, i, GET_NEIGHBOR(nl, ni), ni)
#define GET_ABS_DR_SQ(p, nl, i, ni)  dot_product(GET_DR(p, nl, i, ni), GET_DR(p, nl, i, ni))
#define GET_ABS_DR(p, nl, i, ni)  sqrt(GET_ABS_DR_SQ(p, nl, i, ni))

#define DIST_SQ(p, nl, i, ni, dr, abs_dr)  dr = GET_DR(p, nl, i, ni) ; abs_dr = dot_product(dr, dr)
#define DISTJ_SQ(p, nl, i, ni, j, dr, abs_dr) j = GET_NEIGHBOR(nl, ni) ; dr = GET_DRJ(p, nl, i, j, ni) ; abs_dr = dot_product(dr, dr)

#define DIST(p, nl, i, ni, dr, abs_dr)  dr = GET_DR(p, nl, i, ni) ; abs_dr = sqrt(dot_product(dr, dr))
#define DISTJ(p, nl, i, ni, j, dr, abs_dr) j = GET_NEIGHBOR(nl, ni) ; dr = GET_DRJ(p, nl, i, j, ni) ; abs_dr = sqrt(dot_product(dr, dr))


!>
!! Generate a symmetric pair index of the following form
!!
!!     j   1  2  3  4
!!  i
!!
!!  1      1  2  3  4
!!  2      -  5  6  7
!!  3      -  -  8  9
!!  4      -  -  - 10
!<

#define PAIR_INDEX(i, j, maxval) (1+min((i-1)+(j-1)*maxval, (j-1)+(i-1)*maxval)-min((i-1)*i/2, (j-1)*j/2))



!>
!! Generate a non-symmetric pair index of the following form
!!
!!     j    1   2   3   4
!!  i
!!
!!  1       1   2   3   4
!!  2       5   6   7   8
!!  3       9  10  11  12
!!  4      13  14  15  16
!<

#define PAIR_INDEX_NS(i, j, maxval) (j+(i-1)*maxval)


!>
!! Generate a non-symmetric tribole index
!<

#define TRIPLET_INDEX_NS(i, j, k, maxval) (k+maxval*(j-1+maxval*(i-1)))



#include "error.inc"


!>
!! Assign property *x* to *this%x* if *x* exists. Here, x is an
!! optional argument of a constructor. This is used in the following
!! context to pass arguments to an object upon construction.
!!
!! subroutine class_init(this, ..., x, ...)
!!   ...
!!   real(DP), optional   :: x
!!   ...
!!
!!   ASSIGN_PROPERTY(x)
!!
!!   ...
!! endsubroutine class_init
!!
!! Hence, the property *x* is passed on to the object ONLY IF it is present
!! when init is being called. Otherwise, the default value of *this%x* will
!! be retained.
!<

#define ASSIGN_PROPERTY(x)   if (present(x)) then ; this%x = x ; endif
#define ASSIGN_ARRAY_PROPERTY(src, tar, tar_size)   if (present(src)) then ; tar_size = size(src) ; tar = src ; endif
#define ASSIGN_STRING_ARRAY_PROPERTY(src, tar, tar_size, i)   if (present(src)) then ; tar_size = size(src) ; do i = 1, size(src) ; tar(:, i) = s2a(src(i)) ; enddo ; endif


!>
!! Inlined outer product
!<
#define outer_product(x, y)  ( spread(x, dim=2, ncopies=size(y))*spread(y, dim=1, ncopies=size(x)) )


!>
!! c_loc, but on an array. Circumvents a bug (compiler segfault) in gfortran
!<
#define c_locs(x) c_loc(x(1:1))
#define c_loc1(x) c_loc(x(lbound(x,1)))
#define c_loc11(x) c_loc(x(lbound(x,1),lbound(x,2)))
#define c_loc111(x) c_loc(x(lbound(x,1),lbound(x,2),lbound(x,3)))


!
! Some name mangling
!

#define libAtoms_module supplib


#ifdef LAMMPS

#define SUM_VIRIAL(a, i, b)  a(1, i) = a(1, i) - b(1, 1) ; a(2, i) = a(2, i) - b(2, 2) ; a(3, i) = a(3, i) - b(3, 3) ; a(4, i) = a(4, i) - b(2, 1) ; a(5, i) = a(5, i) - b(3, 1) ; a(6, i) = a(6, i) - b(3, 2)

#else

#define SUM_VIRIAL(a, i, b)  a(1:3, 1:3, i) = a(1:3, 1:3, i) + b

#endif

!>
!! Data type for wavefunctions
!<
#define WF_T real

!>
!! Check association of objects, i.e. raise and error is *a* and *b*
!! are not the same object.
!<
#define ASSERT_ASSOCIATION(a, b, ierror)   if (.not. associated(a, b)) then ; RAISE_ERROR("Internal error: Incompatible object associations.", ierror) ; endif

!>
!! Print a warning
!<
#define WARN(x) \
   write (*, '(4A)')     "Warning (", MODULE_STR, "): ", x ; \
   write (ilog, '(4A)')  "Warning (", MODULE_STR, "): ", x


!>
!! Some LAPACK shortcuts
!<
#define GEMM dgemm
#define MM(n, alpha, A, B, beta, C)  dgemm('N', 'N', n, n, n, alpha, A, n, B, n, beta, C, n, n)


#define DEBUG_WRITE(x)
