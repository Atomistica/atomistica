#ifndef __LOGGING_H
#define __LOGGING_H

#ifdef __cplusplus
extern "C" {
#endif

void prscrlog(const char *msg, ...);
void prlog(const char *msg, ...);

void c_prscrlog(const char *msg);
void c_prlog(const char *msg);

#ifdef __cplusplus
}
#endif

#endif
