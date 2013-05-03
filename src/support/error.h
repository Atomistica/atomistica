/* ======================================================================
   Atomistica - Interatomic potential library
   https://github.com/pastewka/atomistica
   Lars Pastewka, lars.pastewka@iwm.fraunhofer.de, and others
   See the AUTHORS file in the top-level MDCORE directory.

   Copyright (2005-2013) Fraunhofer IWM
   This software is distributed under the GNU General Public License.
   See the LICENSE file in the top-level MDCORE directory.
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
