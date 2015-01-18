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
!! Anderson mixer
!!
!! Produces the Anderson mixing of input vectors in iterative process.
!! See, e.g. V. EYERT in J. Comp. Phys. 124, 271 (1996) (the notation
!! is almost directly from there)
!! x is the input, and y the output vector of iteration. F is the 
!! history of residuals. Mmax is the maximum number of last iterations
!! taken into account. 
!!
!<

#include "macros.inc"

module anderson_mixer
  use supplib

  implicit none

  private

  public :: anderson_mixer_t
  type anderson_mixer_t

     integer                :: n = -1
     integer                :: M = 3

     real(DP), allocatable  :: x_hist(:, :)
     real(DP), allocatable  :: F_hist(:, :)

     real(DP), allocatable  :: xb(:)
     real(DP), allocatable  :: Fb(:)

  endtype anderson_mixer_t


  public :: init
  interface init
     module procedure anderson_mixer_init
  endinterface

  public :: del
  interface del
     module procedure anderson_mixer_del
  endinterface

  public :: set_dimension
  interface set_dimension
     module procedure anderson_mixer_set_dimension
  endinterface

  public :: mix
  interface mix
     module procedure anderson_mixer_mix
  endinterface

contains

  !>
  !! Constructor
  !!
  !! Constructor
  !<
  subroutine anderson_mixer_init(this, M)
    implicit none

    type(anderson_mixer_t), intent(inout)  :: this
    integer, intent(in), optional          :: M

    ! ---

    this%n  = -1

    if (present(M)) then
       this%M  = M
    endif
    
  endsubroutine anderson_mixer_init


  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine anderson_mixer_del(this)
    implicit none

    type(anderson_mixer_t), intent(inout)  :: this

    ! ---

    if (allocated(this%x_hist))  deallocate(this%x_hist)
    if (allocated(this%F_hist))  deallocate(this%F_hist)

    if (allocated(this%xb))      deallocate(this%xb)
    if (allocated(this%Fb))      deallocate(this%Fb)

    this%n = -1
    
  endsubroutine anderson_mixer_del


  !>
  !! Set the dimension of the input/output vectors
  !!
  !! Set the dimension of the input/output vectors
  !<
  subroutine anderson_mixer_set_dimension(this, n)
    implicit none

    type(anderson_mixer_t), intent(inout)  :: this
    integer, intent(in)                    :: n

    ! ---

    call del(this)

    allocate(this%x_hist(0:this%M, n))
    allocate(this%F_hist(0:this%M, n))

    allocate(this%xb(n))
    allocate(this%Fb(n))

    this%x_hist  = 0.0_DP
    this%F_hist  = 0.0_DP

    this%xb      = 0.0_DP
    this%Fb      = 0.0_DP

    this%n  = n
    
  endsubroutine anderson_mixer_set_dimension


  !>
  !! Mixing iteration
  !!
  !! Mixing iteration.
  !<
  subroutine anderson_mixer_mix(this, it, n, xi, yi, beta, limit, done, mx, error)
    implicit none
    
    type(anderson_mixer_t), intent(inout)  :: this   !< Mixer object
    integer, intent(in)                    :: it     !< Iteration
    integer, intent(in)                    :: n      !< Dimension of xi and yi
    real(DP), intent(inout)                :: xi(n)  !< out: Output vector, in: Previous iteration
    real(DP), intent(in)                   :: yi(n)  !< Input vector
    real(DP), intent(in)                   :: beta   !< beta-parameter
    real(DP), optional, intent(in)         :: limit  !< Convergence criterium
    logical,  optional, intent(out)        :: done   !< Is conv. achieved?
    real(DP), optional, intent(in)         :: mx     !< Maximum change for vector elements
    integer,  optional, intent(out)        :: error

    ! ---

    integer   :: i, j, M, ipiv(this%M), info
    real(DP)  :: A(this%M, this%M), b(this%M), hlp

    ! ---

    INIT_ERROR(error)

    if (this%n < n) then
       call set_dimension(this, n)
    endif

    M                    = min(it-1, this%M)
    this%F_hist(0, 1:n)  = yi(1:n) - xi(1:n)   !current residual
    this%x_hist(0, 1:n)  = xi(1:n)             !current input

    b  = 0.0_DP
    A  = 0.0_DP

    !----------------------------
    ! solve A*z=b -> z-->b
    ! (eq.(4.3) in Eyert)
    !----------------------------
    info = 0
    do i = 1, M
       b(i)  = dot_product( this%F_hist(0, 1:n)-this%F_hist(i, 1:n), this%F_hist(0, 1:n) ) 
       do j = 1, M
          A(i, j)  = dot_product( this%F_hist(0, 1:n)-this%F_hist(i, 1:n), this%F_hist(0, 1:n)-this%F_hist(j, 1:n) )
       enddo
    enddo

!#ifdef _MP
!    call sum_in_place(mod_parallel_3d%mpi, b, error=error)
!    PASS_ERROR(error)
!    call sum_in_place(mod_parallel_3d%mpi, A, error=error)
!    PASS_ERROR(error)
!#endif

!#ifdef _MP
!    if (mod_parallel_3d%mpi%my_proc == ROOT) then
!#endif

    if (M > 0) then
       call dgesv(M, 1, A(1:M, 1:M), M, ipiv(1:M), b(1:M), M, info)
    endif

!#ifdef _MP
!    endif
!
!    call bcast(mod_parallel_3d%mpi, info, ROOT, error=error)
!    PASS_ERROR(error)
!#endif

    if (info == 0 .and. M > 0) then
!#ifdef _MP
!       call bcast(mod_parallel_3d%mpi, b, ROOT, error=error)
!       PASS_ERROR(error)
!#endif

       !-----------------------------
       ! We solved the optimum
       ! linear combination b(:)
       !-----------------------------
       this%xb(1:n)  = xi(1:n)
       this%Fb(1:n)  = this%F_hist(0, 1:n)
       do j = 1, M
          this%xb(1:n)  = this%xb(1:n) + b(j) * ( this%x_hist(j, 1:n) - xi(1:n)             )
          this%Fb(1:n)  = this%Fb(1:n) + b(j) * ( this%F_hist(j, 1:n) - this%F_hist(0, 1:n) )
       end do
       this%xb(1:n)  = this%xb(1:n) + beta * this%Fb(1:n)       !next input
    else 
       !----------------------------
       ! The matrix A was singular:
       ! use simple mixing 
       !----------------------------
       this%xb(1:n)  = (1.0_DP-beta)*this%x_hist(0, 1:n) + beta*yi(1:n) !next input
    endif

    !----------------------------------------
    ! The input must not change more than mx
    ! for all elements
    !----------------------------------------
    if (present(mx)) then
       hlp = maxval(abs(this%xb(1:n) - xi(1:n)))
!#ifdef _MP
!       hlp = max(mod_parallel_3d%mpi, hlp, error=error)
!       PASS_ERROR(error)
!#endif
       xi(1:n)  = xi(1:n) + mx/max(mx, hlp) * (this%xb(1:n) - xi(1:n))
    else
       xi(1:n)  = this%xb(1:n)
    endif

    ! shift history
    M                     = min(it, this%M)
    this%F_hist(1:M, 1:n) = this%F_hist(0:M-1, 1:n)
    this%x_hist(1:M, 1:n) = this%x_hist(0:M-1, 1:n)
!    this%F_hist(1:this%M, 1:n) = this%F_hist(0:this%M-1, 1:n)
!    this%x_hist(1:this%M, 1:n) = this%x_hist(0:this%M-1, 1:n)


    !---------------------------------------
    ! convergence: all components of
    ! residual must be less that 'limit'
    !---------------------------------------
    if (present(limit) .and. present(done)) then
       done = all( abs(this%F_hist(0, 1:n)) < limit )
!#ifdef _MP
!       done = all(mod_parallel_3d%mpi, done, error=error)
!       PASS_ERROR(error)
!#endif
    endif

  endsubroutine anderson_mixer_mix

endmodule anderson_mixer
