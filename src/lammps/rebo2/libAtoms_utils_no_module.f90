!! ======================================================================
!! MDCORE - Interatomic potential library
!! https://github.com/pastewka/mdcore
!! Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
!! See the AUTHORS file in the top-level MDCORE directory.
!!
!! Copyright (2005-2013) Fraunhofer IWM
!! This software is distributed under the GNU General Public License.
!! See the LICENSE file in the top-level MDCORE directory.
!! ======================================================================
#include "error.inc"

function quippy_running()
  use system_module, only: get_quippy_running
  logical quippy_running
  quippy_running = get_quippy_running()
end function quippy_running

! Error handling routines callable from C

subroutine c_push_error_with_info(doc, fn, line, kind)
  use error_module
  character(*), intent(in)       :: doc
  character(*), intent(in)       :: fn
  integer, intent(in)            :: line
  integer, intent(in), optional  :: kind

  call push_error_with_info(doc, fn, line, kind)

end subroutine c_push_error_with_info

subroutine c_error_abort(error)
  use error_module
  integer, intent(inout), optional :: error

  call error_abort(error)
end subroutine c_error_abort

