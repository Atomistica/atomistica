#include "error.inc"

subroutine c_system_initialise(verbosity)
  use system_module
  integer verbosity

  call system_initialise(verbosity)

end subroutine c_system_initialise

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

subroutine c_push_error(fn, line, kind)
  use error_module

  character(*), intent(in)       :: fn
  integer, intent(in)            :: line
  integer, intent(in), optional  :: kind

  call push_error(fn, line, kind)

end subroutine c_push_error

subroutine c_error_abort(error)
  use error_module
  integer, intent(inout), optional :: error

  call error_abort(error)
end subroutine c_error_abort

subroutine c_error_clear_stack()
  use error_module
  call error_clear_stack
end subroutine c_error_clear_stack


