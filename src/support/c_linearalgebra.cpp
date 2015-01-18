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

#include <stdio.h>
#include <string.h>

#include "error.h"
#include "logging.h"

#include "linearalgebra.h"

#ifdef HAVE_MKL
#include "mkl_lapack.h"
#endif

/* ----------------------------------------------------------------------
 * print a complex vector to screen
 * --------------------------------------------------------------------*/

void printvec(int dim, double_complex *m)
{
  int dim_sq = dim*dim;
  char str[1024];

  bool real_only = false;

  for (int i = 0; i < dim_sq; i++) {
    if (fabs(cimag(m[i])) > 1e-12)
      real_only = false;
  }

  strcpy(str, " { ");

  for (int l = 0; l < dim; l++) {
    if (real_only) {
      sprintf(str, "%s%10.3e", str, creal(m[l]));
    }
    else {
      sprintf(str, "%s%10.3e+I*%10.3e", str, creal(m[l]), cimag(m[l]));
    }
    if (l != dim-1)  sprintf(str, "%s, ", str);
  }

  printf("%s },\n", str);
}


/* ----------------------------------------------------------------------
 * print a complex matrix to screen
 * --------------------------------------------------------------------*/

void printmat(int dim, double_complex *m)
{
  int dim_sq = dim*dim;
  char str[1024];

  bool real_only = false;

  for (int i = 0; i < dim_sq; i++) {
    if (fabs(cimag(m[i])) > 1e-12)
      real_only = false;
  }

  for (int k = 0; k < dim; k++) {
    if (k == 0) 
      strcpy(str, "{{ ");
    else
      strcpy(str, " { ");

    for (int l = 0; l < dim; l++) {
      if (real_only) {
	sprintf(str, "%s%10.3e", str, creal(m[_IDX2(dim, k, l)]));
      }
      else {
	sprintf(str, "%s%10.3e+I*%10.3e", str, creal(m[_IDX2(dim, k, l)]),
		cimag(m[_IDX2(dim, k, l)]));
      }
      if (l != dim-1)  sprintf(str, "%s, ", str);
    }

    if (k == dim-1)
      printf("%s }}\n", str);
    else
      printf("%s },\n", str);
  }
}


/*!
 * Iterative matrix inversion
 *
 * Invert a matrix using an iterative process. If prev==true, use
 * matrix given in invmat as starting point.
 */
extern "C"
void iterative_matrix_inverse(double *matptr, double *invmatptr, int n,
			      _Bool prev, double epsilon, double *work1,
			      double *work2, int *error,
			      cublasHandle_t cublas_handle, int *nit_out)
{
  INIT_ERROR(error);

  mat<double> matr(n, matptr, cublas_handle);
  mat<double> invmat(n, invmatptr, cublas_handle);
  /* Will allocate and release upon destruction if work1, work2 == NULL */
  mat<double> help1(n, work1, cublas_handle);
  mat<double> help2(n, work2, cublas_handle);

  /*
   * - Initialize inverse matrix if previous not used
   *   The starting invmat has to be small enough so that the iteration
   *   won't start running to infinity
   */

#if 0
  mat<double> dummy(n);
  dummy = matr;

  printf("dummy.data() = %p\n", dummy.data());
  printf("matr.data() = %p\n", matr.data());
  printf("dummy.on_host() = %i\n", dummy.on_host(error));
  PASS_ERROR(error);
  printf("matr.on_host() = %i\n", matr.on_host(error));
  PASS_ERROR(error);
  printf("sum = %f %f\n", dummy.sum(), matr.sum());
  printf("max = %f %f\n", dummy.max(), matr.max());
  printf("min = %f %f\n", dummy.min(), matr.min());
  printf("amax = %f %f\n", dummy.amax(), matr.amax());
  printf("amin = %f %f\n", dummy.amin(), matr.amin());
#endif

  if (!prev) {
    double smin, smax;
    ev_bounds(n, matptr, &smin, &smax, error);
    PASS_ERROR(error);
    mat_mul_sca(1.0/(n*MAX(fabs(smin), fabs(smax))), matr, invmat, error);
    PASS_ERROR(error);
  }

  /*
   * Find inverse via S^-1 = 2 S^1 - S^-1 S S^-1
   */

  double sigma = epsilon + 1.0;
  int i = 0;
  while (sigma > epsilon) {

    /*
     * help1 = matr.invmat
     */

    gemm(OP_N, OP_N, 1.0, matr, invmat, 0.0, help1, error);
    PASS_ERROR(error);

    help2 = invmat;

    /*
     * invmat = -help2.help1 + 2*invmat
     */

    gemm(OP_N, OP_N, -1.0, help2, help1, 2.0, invmat, error);
    PASS_ERROR(error);

    mat_mul_sca(1.0, help2, -1.0, invmat, help1, error);
    PASS_ERROR(error);

    sigma = help1.amax(error);
    PASS_ERROR(error);
    i = i+1;

    if (i % 100 == 0) {
      prscrlog("iterative_matrix_inverse: No convergence after %i iterations.",
	       i);
    }

  }

  if (nit_out) {
    *nit_out = i;
  }
}


extern "C"
void dev_bounds(int n, double *H, double *l, double *u)
{
  ev_bounds(n, H, l, u);
}

