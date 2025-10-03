/* Compatibility header for NumPy 1.x and 2.x API changes */

#ifndef NUMPY_COMPAT_H
#define NUMPY_COMPAT_H

#include <numpy/arrayobject.h>

/* NumPy 2.0 removed several deprecated constants */
#if NPY_ABI_VERSION < 0x02000000

/* NumPy 1.x - constants exist, no need to redefine */

#else

/* NumPy 2.0+ - define removed constants using new API */

/* NPY_FARRAY was (NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED) */
#ifndef NPY_FARRAY
#define NPY_FARRAY (NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED)
#endif

/* NPY_BEHAVED was (NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE) */
#ifndef NPY_BEHAVED
#define NPY_BEHAVED (NPY_ARRAY_ALIGNED | NPY_ARRAY_WRITEABLE)
#endif

/* NPY_C_CONTIGUOUS is now NPY_ARRAY_C_CONTIGUOUS */
#ifndef NPY_C_CONTIGUOUS
#define NPY_C_CONTIGUOUS NPY_ARRAY_C_CONTIGUOUS
#endif

/* NPY_DEFAULT was removed - use 0 or NPY_ARRAY_ENSUREARRAY */
#ifndef NPY_DEFAULT
#define NPY_DEFAULT 0
#endif

#endif /* NPY_ABI_VERSION */

#endif /* NUMPY_COMPAT_H */
