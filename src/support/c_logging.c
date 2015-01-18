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
#include <stdarg.h>

#include "logging.h"


/*!
 * Record a log message to screen and file
 */
void prscrlog(const char *msg, ...)
{
  char buf[1024];
  va_list args;
  va_start(args, msg);
  vsprintf(buf, msg, args);
  va_end(args);
  c_prscrlog(buf);
}


/*!
 * Record a log message to file only
 */
void prlog(const char *msg, ...)
{
  char buf[1024];
  va_list args;
  va_start(args, msg);
  vsprintf(buf, msg, args);
  va_end(args);
  c_prlog(buf);
}
