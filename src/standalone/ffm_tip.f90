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
!   private
!   classtype:ffm_tip_t classname:FFMTip interface:potentials
! @endmeta

#include "macros.inc"

module ffm_tip
  use supplib

  use io
  use logging
  use timer

  use data
  use particles
  use neighbors
  use dynamics

  implicit none

  private

  public :: ffm_tip_t

  type ffm_tip_t
    real(DP) :: abs_v
    real(DP) :: angle       =  0.0_DP

    real(DP) :: k_xy        =  0.0_DP
    real(DP) :: k_z         =  0.0_DP
    real(DP) :: gamma_xy    =  0.0_DP
    real(DP) :: gamma_z     =  0.0_DP

    real(DP) :: f_z         =  0.0_DP
    real(DP) :: r_ghost(3)
    real(DP) :: v_ghost(3)
    real(DP) :: f_norm(3)
    real(DP) :: mass_ghost  =  1.0_DP
    real(DP) :: gamma_ghost =  0.0_DP
    real(DP) :: sl_dir(3)

    integer              :: top_group
    integer              :: n_top
    integer, allocatable :: top_atoms(:)
    real(DP)             :: mass         = -1.0_DP
    real(DP)             :: r_init(3)

    integer  :: un
    real(DP) :: time     =  0.0_DP
    real(DP) :: log_freq = -1.0_DP

    character(100) :: op_mode
    character(100) :: height_mode

  endtype ffm_tip_t


  public :: del
  interface del
    module procedure ffm_tip_del
  endinterface

  public :: bind_to
  interface bind_to
    module procedure ffm_tip_bind_to
  endinterface

  public :: adjust_velocities_and_forces
  interface adjust_velocities_and_forces
    module procedure ffm_tip_adjust_velocities_and_forces
  endinterface

  public :: energy_and_forces_with_dyn
  interface energy_and_forces_with_dyn
    module procedure ffm_tip_energy_and_forces
  endinterface

  public :: register
  interface register
    module procedure ffm_tip_register
  endinterface

contains

  !>
  !! Destructor
  !!
  !! Destructor
  !<
  subroutine ffm_tip_del(this)
    implicit none

    type(ffm_tip_t), intent(inout)  :: this

    ! ---

    if (allocated(this%top_atoms)) then
      deallocate(this%top_atoms)
    end if

    if (this%log_freq > 0.0_DP) then
      call fclose(this%un)
    end if

  endsubroutine ffm_tip_del


  subroutine ffm_tip_bind_to(this, p, nl, ierror)
    implicit none

    type(ffm_tip_t), intent(inout)    :: this
    type(particles_t), intent(inout)  :: p
    type(neighbors_t), intent(in)     :: nl
    integer, intent(inout), optional  :: ierror

    ! ---

    integer  :: i,j
    real(DP) :: r_cm(3)

    ! ---

    if (trim(this%op_mode) == 'scan') then
      if (trim(this%height_mode) /= 'const_height'        .and. &
          trim(this%height_mode) /= 'const_force'         .and. &
          trim(this%height_mode) /= 'spring_const_height' .and. &
          trim(this%height_mode) /= 'spring_const_force') then
        RAISE_ERROR("Unknown height-mode.", ierror)
      end if
    end if

    this%angle = this%angle/180.0_DP*PI

    this%v_ghost = (/this%abs_v * cos(this%angle), this%abs_v * sin(this%angle), 0.0_DP/)

    this%f_norm = (/0.0_DP, 0.0_DP, this%f_z/)

    if (any(this%v_ghost(:) > 0.0_DP)) then
      this%sl_dir = this%v_ghost / sqrt(dot_product(this%v_ghost, this%v_ghost))
    else
      this%sl_dir = 0.0_DP
    end if

    do i=1,p%nat
      if (p%g(i) == this%top_group) then
        this%n_top = this%n_top + 1
      end if
    end do

    if (this%n_top == 0) then
       RAISE_ERROR("No *top* atoms found.", ierror)
    endif

    allocate(this%top_atoms(this%n_top))

    j = 1
    do i=1,p%nat
      if (p%g(i) == this%top_group) then
        this%top_atoms(j) = i
        r_cm   = r_cm + POS3(p,i)
        j = j+1
      end if
    end do

    r_cm = r_cm / float(this%n_top)

    this%r_init  = r_cm
    this%r_ghost = r_cm



    call prlog("- FFM_tip_init -")
    call prlog(" op_mode      = " // this%op_mode)

    if (trim(this%op_mode) == 'scan') then
      call prlog(" height_mode  = " // this%height_mode)
    end if

    call prlog(" top          = " // this%top_group)

    if (this%mass < 0.0_DP) then
      this%mass = 0.0_DP
      do i=1,this%n_top
        this%mass = this%mass + p%m(this%top_atoms(i))
      end do
    else
      call prlog(" mass         = " // this%mass)
      do i=1,this%n_top
        p%m(this%top_atoms(i)) = this%mass/float(this%n_top)
      end do
    end if

    if (trim(this%op_mode) == 'indent_relax') then
      do i=1,this%n_top
        p%m(this%top_atoms(i)) = p%m(this%top_atoms(i)) * float(this%n_top)
      end do
    end if


    if (trim(this%op_mode) == 'indent' .or. trim(this%op_mode) == 'indent_relax') then
      this%height_mode = 'const_force'
      call prlog(" Fz           = " // this%f_z)
      call prlog(" gamma_z      = " // this%gamma_z)
    end if

    if (trim(this%op_mode) == 'approach') then
      this%height_mode = 'spring_const_force'
      call prlog(" mass_ghost   = " // this%mass_ghost)
      call prlog(" gamma_ghost  = " // this%gamma_ghost)
      call prlog(" Fz           = " // this%f_z)
      call prlog(" k_z          = " // this%k_z)
      call prlog(" gamma_z      = " // this%gamma_z)
    end if

    if (trim(this%op_mode) == 'scan') then
      if (trim(this%height_mode) == 'const_height') then
        call prlog(" v            = " // this%abs_v)
        call prlog(" angle        = " // this%angle)
      end if
      if (trim(this%height_mode) == 'const_force') then
        call prlog(" Fz           = " // this%f_z)
        call prlog(" gamma_z      = " // this%gamma_z)
        call prlog(" v            = " // this%abs_v)
        call prlog(" angle        = " // this%angle)
      end if
      if (trim(this%height_mode) == 'spring_const_height') then
        call prlog(" k_z          = " // this%k_z)
        call prlog(" gamma_z      = " // this%gamma_z)
        call prlog(" v            = " // this%abs_v)
        call prlog(" angle        = " // this%angle)
      end if
      if (trim(this%height_mode) == 'spring_const_force') then
        call prlog(" mass_ghost   = " // this%mass_ghost)
        call prlog(" gamma_ghost  = " // this%gamma_ghost)
        call prlog(" Fz           = " // this%f_z)
        call prlog(" k_z          = " // this%k_z)
        call prlog(" gamma_z      = " // this%gamma_z)
        call prlog(" v            = " // this%abs_v)
        call prlog(" angle        = " // this%angle)
      end if
    end if


    if (this%log_freq > 0.0_DP) then
      this%un = fopen("ffm_tip.out", F_WRITE)
      write(this%un, '(A6,14X,12(4X,A20))')  "#01:ti", "02:rx", "03:ry", "04:rz", "05:vx", "06:vy", "07:vz" , "08:fx", "09:fy", "10:fz", "11:fsx", "12:fsy", "13:fz"
    end if

  endsubroutine ffm_tip_bind_to


  subroutine ffm_tip_adjust_velocities_and_forces(this, p, v, f, ti, dt)
    implicit none

    type(ffm_tip_t), intent(inout)  :: this
    type(particles_t), intent(in)   :: p
    real(DP), intent(inout)         :: v(3, p%maxnatloc)
    real(DP), intent(inout)         :: f(3, p%maxnatloc)
    real(DP), intent(in)            :: ti
    real(DP), intent(in)            :: dt

    ! ---

    integer  :: i,j
    real(DP) :: r_cm(3)
    real(DP) :: v_cm(3)
    real(DP) :: f_cm(3)
    real(DP) :: f_cm_no_mod(3)
    real(DP) :: f_cm_log(3)
    real(DP) :: dr(3)
    real(DP) :: dr_abs
    real(DP) :: dr_n(3)
    real(DP) :: f_harm
    real(DP) :: f_damp(3)
    real(DP) :: f_ghost_z

    ! ---


    r_cm        = 0.0_DP
    v_cm        = 0.0_DP
    f_cm        = 0.0_DP
    f_cm_no_mod = 0.0_DP

    do i=1,this%n_top
      j = this%top_atoms(i)
      r_cm(:) = r_cm(:) + p%r_cont(:,j)
      v_cm(:) = v_cm(:) + v(:,j)
      f_cm(:) = f_cm(:) + f(:,j)
    end do
    f_cm_no_mod = f_cm

    r_cm = r_cm / float(this%n_top)
    v_cm = v_cm / float(this%n_top)

    this%r_ghost    = this%r_ghost + this%v_ghost * dt


    ! Constant force applied on ghost atom in z-direction:
    ! use Euler-algorithm to propagate ghost atom in z-direction
    if (trim(this%height_mode) == 'spring_const_force') then
      f_ghost_z = this%k_z * (r_cm(3) - this%r_ghost(3)) + this%f_norm(3) - this%gamma_ghost * this%v_ghost(3)
      f_cm(3) = f_cm(3) - this%k_z * (r_cm(3) - this%r_ghost(3))

      this%r_ghost(3) = this%r_ghost(3) + this%v_ghost(3) * dt + 0.5_DP * f_ghost_z/this%mass_ghost * dt**2
      this%v_ghost(3) = this%v_ghost(3) + f_ghost_z/this%mass_ghost * dt
    end if


    ! Ghost atom on constant height: spring couples ghost
    ! atom and center of mass of rigid layer in z-direction
    if (trim(this%height_mode) == 'spring_const_height') then
      f_cm(3) = f_cm(3) - this%k_z * (r_cm(3) - this%r_ghost(3))
    end if


    ! distance between ghost point and center of
    ! mass of rigid layer parallel to surface
    dr     = this%r_ghost - r_cm
    dr(3)  = 0.0_DP
    dr_abs = sqrt(dr(1)**2 + dr(2)**2)

    if (dr_abs .gt. 1d-10) then
      dr_n = dr/dr_abs
    else
      dr_n = 0.0_DP
    end if


    f_harm = this%k_xy * dr_abs

    f_damp = (/-this%gamma_xy*v_cm(1), -this%gamma_xy*v_cm(2), -this%gamma_z*v_cm(3)/)


    ! adjust forces on center of mass of rigid layer
    if (trim(this%op_mode) == 'scan') then
      ! movements perpendicular to the sliding direction are forbidden
      f_cm = f_cm + f_harm*dr_n + f_damp
      v_cm = dot_product(v_cm, this%sl_dir) * this%sl_dir + (/0.0_DP, 0.0_DP, v_cm(3)/)
      f_cm = dot_product(f_cm, this%sl_dir) * this%sl_dir + (/0.0_DP, 0.0_DP, f_cm(3)/)
    else
      ! indent or approach - ony movements in z-direction allowed
      f_cm = f_cm + f_damp
      v_cm = (/0.0_DP, 0.0_DP, v_cm(3)/)
      f_cm = (/0.0_DP, 0.0_DP, f_cm(3)/)
    end if


    ! adjust velocities and forces in z-direction
    if (trim(this%height_mode) == 'const_force') then
      ! add constant force in z-direction
      v_cm = (/this%abs_v * cos(this%angle), this%abs_v * sin(this%angle), v_cm(3)/)
      f_cm = f_cm + this%f_norm
    else if (trim(this%height_mode) == 'const_height') then
      ! no movements in z-direction allowed
      v_cm = (/this%abs_v * cos(this%angle), this%abs_v * sin(this%angle), 0.0_DP/)
      f_cm = (/0.0_DP, 0.0_DP, 0.0_DP/)
    end if


    do i=1,this%n_top
      j = this%top_atoms(i)
      v(:,j) = v_cm(:)

      if (trim(this%op_mode) == 'indent_relax') then
        f(:,j) = f_cm(:)
      else
        f(:,j) = f_cm(:)/float(this%n_top)
      end if
    end do


    this%time = this%time + dt
    
    if (trim(this%height_mode) == 'const_height') then
      f_cm_log = f_cm_no_mod
    else
      f_cm_log = f_cm
    end if

    if (this%log_freq > 0.0_dp .and. this%time > this%log_freq) then
      write(this%un, '(E20.10,12(4X,E20.10))')  ti, r_cm(1), r_cm(2), r_cm(3), v_cm(1), v_cm(2), v_cm(3) , f_cm_log(1), f_cm_log(2), f_cm_log(3), f_harm*dr_n(1), f_harm*dr_n(2), this%f_norm(3)
      this%time = 0.0_DP
    end if

  endsubroutine ffm_tip_adjust_velocities_and_forces


  subroutine ffm_tip_energy_and_forces(this, dyn, nl, ierror)
    implicit none

    type(ffm_tip_t), intent(inout)     :: this
    type(dynamics_t), target           :: dyn 
    type(neighbors_t), intent(in)      :: nl
    integer, intent(inout), optional   :: ierror

    ! ---

    call adjust_velocities_and_forces(this, dyn%p, dyn%v, dyn%f, dyn%ti, dyn%dt)

  endsubroutine ffm_tip_energy_and_forces

  subroutine ffm_tip_register(this, cfg, m)
    use, intrinsic :: iso_c_binding

    implicit none

    type(ffm_tip_t), target, intent(inout)  :: this
    type(c_ptr), intent(in)            :: cfg
    type(c_ptr), intent(out)           :: m

    ! ---

    m = ptrdict_register_section(cfg, CSTR("FFMTip"), &
         CSTR("Couple a set of atoms to a cantilever."))

    call ptrdict_register_string_property(m, c_loc(this%op_mode), 100, &
         CSTR("op_mode"), &
         CSTR("Operation mode of the tip."))
    call ptrdict_register_string_property(m, c_loc(this%height_mode), 100, &
         CSTR("height_mode"), &
         CSTR("Height mode of the tip."))

    call ptrdict_register_integer_property(m, c_loc(this%top_group), &
         CSTR("top"), &
         CSTR("Top atoms (group)."))
    call ptrdict_register_real_property(m, c_loc(this%mass), CSTR("mass"), &
         CSTR("Total mass of top atoms."))
    call ptrdict_register_real_property(m, c_loc(this%f_z), CSTR("Fz"), &
         CSTR("Normal force."))
    call ptrdict_register_real_property(m, c_loc(this%abs_v), CSTR("v"), &
         CSTR("Cantilever velocity."))
    call ptrdict_register_real_property(m, c_loc(this%k_xy), CSTR("k_xy"), &
         CSTR("Spring constant (x/y-direction)."))
    call ptrdict_register_real_property(m, c_loc(this%k_z), CSTR("k_z"), &
         CSTR("Spring constant (z-direction)."))
    call ptrdict_register_real_property(m, c_loc(this%gamma_xy), &
         CSTR("gamma_xy"), &
         CSTR("Damping constant (x-/y-direction)."))
    call ptrdict_register_real_property(m, c_loc(this%gamma_z), &
         CSTR("gamma_z"), &
         CSTR("Damping constant (z-direction)."))
    call ptrdict_register_real_property(m, c_loc(this%angle), CSTR("angle"), &
         CSTR("Angle of v relative to the x-axis (in degrees)."))
    call ptrdict_register_real_property(m, c_loc(this%mass_ghost), &
         CSTR("mass_ghost"), &
         CSTR("Mass of ghost atom."))
    call ptrdict_register_real_property(m, c_loc(this%gamma_ghost), &
         CSTR("gamma_ghost"), &
         CSTR("Damping constant for ghost atom movements."))

    call ptrdict_register_real_property(m, c_loc(this%log_freq), &
         CSTR("log_freq"), &
         CSTR("Logfile output frequency."))

  endsubroutine ffm_tip_register

endmodule ffm_tip

