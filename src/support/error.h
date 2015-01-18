/* ======================================================================
   Atomistica - Interatomic potential library and molecular dynamics code
   https://github.com/Atomistica/atomistica

   Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
   See the AUTHORS file in the top-level Atomistica directory.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 2 of the License, or
   (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   ====================================================================== */
#ifndef __ERROR_H
#define __ERROR_H

#define ERROR_NONE                0
#define ERROR_UNSPECIFIED        -1
#define ERROR_IO                 -2
#define ERROR_IO_EOF             -3

#define INIT_ERROR(error) if (error != NULL) *error = ERROR_NONE
#define RAISE_ERROR(error, info, ...) { char error_h_info[1000]; sprintf(error_h_info, info, ## __VA_ARGS__ ); int error_h_line = __LINE__; int error_h_kind = ERROR_UNSPECIFIED; c_push_error_with_info(error_h_info, __FILE__, error_h_line, error_h_kind); if (error != NULL) { *error = error_h_kind; return; } else c_error_abort(error_h_kind); }
#define RAISE_ERROR_WITH_RET(error, x, info, ...) { char error_h_info[1000]; sprintf(error_h_info, info, ## __VA_ARGS__ ); int error_h_line = __LINE__; int error_h_kind = ERROR_UNSPECIFIED; c_push_error_with_info(error_h_info, __FILE__, error_h_line, error_h_kind); if (error != NULL) { *error = error_h_kind; return (x); } else c_error_abort(error_h_kind); }
#define RAISE_ERROR_WITH_KIND(error, kind, info, ...) { char error_h_info[1000]; sprintf(error_h_info, info, ## __VA_ARGS__ ); int error_h_line = __LINE__; int error_h_kind = kind; c_push_error_with_info(error_h_info, __FILE__, error_h_line, error_h_kind); if (error != NULL) { *error = error_h_kind; return; } else c_error_abort(error_h_kind); }
#define PASS_ERROR(error) if (error != NULL && *error != ERROR_NONE) { int error_h_line = __LINE__; c_push_error(__FILE__, error_h_line, *error); return; }
#define PASS_ERROR_WITH_RET(error, x) if (error != NULL && *error != ERROR_NONE) { int error_h_line = __LINE__; c_push_error(__FILE__, error_h_line, *error); return (x); }
#define PASS_PYTHON_ERROR(error, res) { if (!res) { py_to_error(__FILE__, __LINE__, error); return; } }
#define CLEAR_ERROR error_clear_stack();

#ifdef __cplusplus
extern "C" {
#endif

void c_push_error_with_info(const char *, const char *, int, int);
void c_push_error(const char*, int, int);
void error_clear_stack(void);
void c_error_abort(int);

#ifdef __cplusplus
};
#endif

#endif
