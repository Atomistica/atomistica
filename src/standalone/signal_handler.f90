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
module signal_handler
  implicit none

  interface
     subroutine sigreg(signum, func) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       integer(C_INT), value :: signum
       type(C_PTR),    value :: func
     endsubroutine sigreg

     subroutine sigclear(signum) bind(C)
       use, intrinsic :: iso_c_binding

       implicit none

       integer(C_INT), value :: signum
     endsubroutine sigclear
  endinterface

  !
  ! Run control
  !

  logical, save :: done

  contains

  function handle_signal(signum) bind(C)
    !dir$ attributes default :: hangle_signal
    use, intrinsic :: iso_c_binding

    implicit none

    integer(C_INT), value  :: signum
    integer(C_INT)         :: handle_signal

    ! ---

    write (*, '(A,I2,A)')     "RECEIVED SIGNAL ", signum, "; USER REQUESTED ABORT"

    done  = .true.

    handle_signal  = 1

  endfunction handle_signal

endmodule signal_handler
