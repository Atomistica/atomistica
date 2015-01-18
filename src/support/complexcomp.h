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

#ifndef __COMPLEXCOMP_H
#define __COMPLEXCOMP_H

/*
 * Compatibility with different complex data structures
 */

#include <complex>

#ifdef C99_COMPLEX

/*
 * C99 complex type
 */

#include <complex.h>

typedef double complex double_complex;

#define cnorm(x)  (creal(x)*creal(x)+cimag(x)*cimag(x))

#define COMPLEX_NUMBER(a, b)  ( (a) + I*(b) )

#else

#ifdef FFTW_COMPLEX

/*
 * FFTW complex type
 */

#include <fftw3.h>

typedef fftw_complex double_complex;

#define COMPLEX_NUMBER(a, b)  { (a), (b) }

#define creal(x)  (x[0])
#define cimag(x)  (x[1])
#define cnorm(x)  (creal(x)*creal(x)+cimag(x)*cimag(x))
#define cabs(x)   sqrt(cnorm(x))
#define conj(x)   COMPLEX_NUMBER(creal(x), -cimag(x))

#else

#ifdef HAVE_MKL

#include <mkl.h>

typedef MKL_Complex16 double_complex;

#define COMPLEX_NUMBER(a, b)  ((double_complex) { (a), (b) })

inline double creal(double_complex x) { return x.real; };
inline double cimag(double_complex x) { return x.imag; };
inline double cnorm(double_complex x) { return creal(x)*creal(x)+cimag(x)*cimag(x); };
inline double cabs(double_complex x) { return sqrt(cnorm(x)); };
inline double abs(double_complex x) { return cabs(x); };
inline double_complex conj(double_complex x) {
  return COMPLEX_NUMBER(creal(x), -cimag(x));
};

inline double_complex operator+(double_complex x, double_complex y) {
  return COMPLEX_NUMBER(creal(x)+creal(y), cimag(x)+cimag(y));
}
inline double_complex operator-(double_complex x, double_complex y) {
  return COMPLEX_NUMBER(creal(x)-creal(y), cimag(x)-cimag(y));
}
inline double_complex operator*(double_complex x, double_complex y) {
  return COMPLEX_NUMBER(creal(x)*creal(y)+cimag(x)*cimag(y),
			-creal(x)*cimag(y)-cimag(x)*creal(y));
}

#else

#ifdef HAVE_CUDA

#include "cublas.h"

typedef cuDoubleComplex double_complex;

#define COMPLEX_NUMBER(a, b)  ((double_complex) { (a), (b) })

inline double creal(double_complex x) { return x.x; };
inline double cimag(double_complex x) { return x.y; };
inline double cnorm(double_complex x) { return creal(x)*creal(x)+cimag(x)*cimag(x); };
inline double cabs(double_complex x) { return sqrt(cnorm(x)); };
inline double abs(double_complex x) { return cabs(x); };
inline double_complex conj(double_complex x) {
  return COMPLEX_NUMBER(creal(x), -cimag(x));
};

inline double_complex operator+(double_complex x, double_complex y) {
  return COMPLEX_NUMBER(creal(x)+creal(y), cimag(x)+cimag(y));
}
inline double_complex operator-(double_complex x, double_complex y) {
  return COMPLEX_NUMBER(creal(x)-creal(y), cimag(x)-cimag(y));
}
inline double_complex operator*(double_complex x, double_complex y) {
  return COMPLEX_NUMBER(creal(x)*creal(y)+cimag(x)*cimag(y),
			-creal(x)*cimag(y)-cimag(x)*creal(y));
}

#else

/*
 * C++ complex type is default
 */

#include <complex>

typedef std::complex<double> double_complex;

#define cabs(x)   std::abs(x)
#define creal(x)  std::real(x)
#define cimag(x)  std::imag(x)
#define conj(x)   std::conj(x)
#define cnorm(x)  std::norm(x)
#define cexp(x)   std::exp(x)

#define I double_complex(0.0, 1.0);

#define COMPLEX_NUMBER(a, b)  double_complex(a, b)

#endif

#endif

#endif

#endif

#endif
