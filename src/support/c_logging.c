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
