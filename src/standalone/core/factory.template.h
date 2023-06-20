/*
 * %(disclaimer)s
 */

#ifndef __%(name)s_DISPATCH_H_
#define __%(name)s_DISPATCH_H_

#include "ptrdict.h"

#define N_CLASSES %(n_classes)i

/*
 * Class definition
 */

typedef struct __%(name)s_class_t {

  char name[MAX_NAME+1];
  void (*new_instance)(void *, section_t *, section_t **);

} %(name)s_class_t;

extern %(name)s_class_t %(name)s_classes[N_CLASSES];

#endif


