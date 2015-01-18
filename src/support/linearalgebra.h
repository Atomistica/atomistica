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
#ifndef __LINEARALGEBRA_H
#define __LINEARALGEBRA_H

#include <math.h>
#include <stdbool.h>

#include "error.h"
#include "complexcomp.h"

#ifdef HAVE_MKL
#include <mkl.h>
#else

/* FIXME: This the right place? */
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

/* FIXME: When do we actually need to do this? */
#define dgemm dgemm_
#define zgemm zgemm_

extern "C" void dgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                      const double *alpha, double *A, const int *ldA, double *B, const int *lbB, double *beta,
                      double *C, const int *ldC);
extern "C" void zgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
                      const double_complex *alpha, double_complex *A, const int *ldA, double_complex *B,
                      const int *lbB, double_complex *beta, double_complex *C, const int *ldC);
#endif

#include "cu_vec.h"
#include "mat.h"
#include "vec.h"

enum op_t { OP_N = 0, OP_T = 1, OP_C = 2 };
const char blas_op[][2]  = { "N", "T", "C" };
#ifdef HAVE_CUDA
const cublasOperation_t cublas_op[] = { CUBLAS_OP_N, CUBLAS_OP_T, CUBLAS_OP_C };
#endif

/*!
 * invert matrix
 */
extern "C"
void iterative_matrix_inverse(double *mat, double *invmat, int n, _Bool prev,
			      double epsilon, double *work1=NULL,
			      double *work2=NULL, int *error=NULL,
			      cublasHandle_t cublas_handle=NULL,
			      int *nit_out=NULL);

/*!
 * get bounds on the eigenvalues
 */
extern "C"
void dev_bounds(int n, double *H, double *l, double *u);
extern "C"
void zev_bounds(int n, double_complex *H, double *l, double *u);

/*!
 * inlined method, invert a 3x3 matrix
 */
template<typename T>
inline void invert3x3(T *A, int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(A)) {
#endif
    T Ainv[9];

#define I3(i, j)  ((j-1)*3+(i-1))

    Ainv[I3(1,1)] = A[I3(2,2)]*A[I3(3,3)] - A[I3(3,2)]*A[I3(2,3)];
    Ainv[I3(1,2)] = A[I3(3,2)]*A[I3(1,3)] - A[I3(1,2)]*A[I3(3,3)];
    Ainv[I3(1,3)] = A[I3(1,2)]*A[I3(2,3)] - A[I3(1,3)]*A[I3(2,2)];
    T detA = Ainv[I3(1,1)]*A[I3(1,1)] + Ainv[I3(1,2)]*A[I3(2,1)] +
      Ainv[I3(1,3)]*A[I3(3,1)];

    Ainv[I3(1,1)] = Ainv[I3(1,1)]/detA;
    Ainv[I3(1,2)] = Ainv[I3(1,2)]/detA;
    Ainv[I3(1,3)] = Ainv[I3(1,3)]/detA;

    Ainv[I3(2,1)] = (A[I3(2,3)]*A[I3(3,1)] - A[I3(2,1)]*A[I3(3,3)])/detA;
    Ainv[I3(2,2)] = (A[I3(1,1)]*A[I3(3,3)] - A[I3(3,1)]*A[I3(1,3)])/detA;
    Ainv[I3(2,3)] = (A[I3(2,1)]*A[I3(1,3)] - A[I3(1,1)]*A[I3(2,3)])/detA;

    Ainv[I3(3,1)] = (A[I3(2,1)]*A[I3(3,2)] - A[I3(2,2)]*A[I3(3,1)])/detA;
    Ainv[I3(3,2)] = (A[I3(3,1)]*A[I3(1,2)] - A[I3(1,1)]*A[I3(3,2)])/detA;
    Ainv[I3(3,3)] = (A[I3(1,1)]*A[I3(2,2)] - A[I3(1,2)]*A[I3(2,1)])/detA;

    for (int i = 0; i < 9; i++)
      A[i] = Ainv[i];

#undef I3
#ifdef HAVE_CUDA
  }
  else {
    RAISE_ERROR(error, "On-device *invert3x3* not yet supported.");
  }
#endif
}


/*!
 * inlined method, invert a 3x3 matrix
 */
template<typename T>
inline void invert3x3(mat<T> &A, int *error=NULL)
{
  invert3x3(A.data(), error);
}


/*!
 * inlined method, transpose of a matrix
 */
template<typename T>
inline void transpose(int dim, T *mat1, T *mat2, int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(mat1) && cu_on_host(mat2)) {
#endif
    int m = 0;
    for (int i = 0; i < dim; i++){
      int n = i;
      for (int j = 0; j < dim; j++, m++, n+=dim)
	mat1[m] = mat2[n];
    }
#ifdef HAVE_CUDA
  }
  else {
    RAISE_ERROR(error, "On-device *transpose* not yet supported.");
  }
#endif
}


/*!
 * inlined method, transpose of a matrix
 */
template<typename T>
inline void transpose(mat<T> &A, mat<T> &B, int *error=NULL)
{
  INIT_ERROR(error);
  if (A.dim() != B.dim()) {
    RAISE_ERROR(error, "Matrices not aligned.");
  }
  transpose(A.dim(), A.data(), B.data(), error);
  PASS_ERROR(error);
}


/*!
 * inlined method, to do vector dot product
 */
template<typename T>
inline T dot(int dim, T *v1, T *v2, int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(v1) && cu_on_host(v2)) {
#endif
    T d = 0.0;
    for (int i = 0; i < dim; i++) {
      d += v1[i]*v2[i];
    }
    return d;
#ifdef HAVE_CUDA
  }
  else {
    RAISE_ERROR(error, "On-device *dot* not yet supported.");
  }
#endif
}


/*!
 * inlined method, to do matrix-vector multiplication for
 * matrix and vector; square matrix is required and output vector is
 * overwritten.
 */
template<typename T>
inline void gemv(int dim, T alpha, T *Mat, T *Vin, T beta, T *Vout,
		 int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(Mat) && cu_on_host(Vin) && cu_on_host(Vout)) {
#endif
    int m=0;
    for (int i=0; i<dim; i++){
      Vout[i] = 0.0;
      for (int j=0; j<dim; j++) Vout[i] = alpha*Mat[m++]*Vin[j] + beta*Vout;
    }
#ifdef HAVE_CUDA
  }
  else {
    RAISE_ERROR(error, "On-device *gemv* not yet supported.");
  }
#endif
}


/*!
 * inlined method, to do matrix-matrix multiplication for
 * matrix; square matrix is required. For double data type.
 */
inline void gemm(op_t op_A, op_t op_B, int dim, double alpha, double *A,
		 double *B, double beta, double *C, int *error = NULL,
		 cublasHandle_t cublas_handle = NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  bool hA = cu_on_host(A), hB = cu_on_host(B), hC = cu_on_host(C);
  if (hA && hB && hC) {
#endif
    dgemm(blas_op[op_A], blas_op[op_B], &dim, &dim, &dim, &alpha, A, &dim, B,
	  &dim, &beta, C, &dim);
#ifdef HAVE_CUDA
  }
  else if (!hA && !hB && !hC) {
    if (!cublas_handle) {
      RAISE_ERROR(error, "Please provide a CUBLAS handle for on-device gemm.");
    }
    PASS_CUBLAS_ERROR( error,
		       cublasDgemm(cublas_handle, cublas_op[op_A],
				   cublas_op[op_B], dim, dim, dim, &alpha, A,
				   dim, B, dim, &beta, C, dim) );
  }
  else {
    RAISE_ERROR(error, "Mixed device/host *gemm* not yet supported.");
  }
#endif
}


/*!
 * inlined method, to do matrix-matrix multiplication for
 * matrix; square matrix is required. For double complex data type.
 */
inline void gemm(op_t op_A, op_t op_B, int dim, double_complex alpha,
		 double_complex *A, double_complex *B, double_complex beta,
		 double_complex *C, int *error = NULL,
		 cublasHandle_t cublas_handle = NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  bool hA = cu_on_host(A), hB = cu_on_host(B), hC = cu_on_host(C);
  if (hA && hB && hC) {
#endif
    zgemm(blas_op[op_A], blas_op[op_B], &dim, &dim, &dim, &alpha, A, &dim, B,
	  &dim, &beta, C, &dim);
#ifdef HAVE_CUDA
  }
  else if (!hA && !hB && !hC) {
    if (!cublas_handle) {
      RAISE_ERROR(error, "Please provide a CUBLAS handle for on-device gemm.");
    }
    PASS_CUBLAS_ERROR( error,
		       cublasZgemm(cublas_handle, cublas_op[op_A],
				   cublas_op[op_B], dim, dim, dim,
				   (cuDoubleComplex *) &alpha,
				   (cuDoubleComplex *) A, dim, 
				   (cuDoubleComplex *) B, dim,
				   (cuDoubleComplex *) &beta,
				   (cuDoubleComplex *) C, dim) );
  }
  else {
    RAISE_ERROR(error, "On-device *gemm* not yet supported.");
  }
#endif
}


/*!
 * inlined method, to do matrix-matrix multiplication for
 * matrix; square matrix is required.
 */
template<typename T>
inline void gemm(op_t op_A, op_t op_B, T alpha, mat<T> &A, mat<T> &B, T beta,
		 mat<T> &C, int *error=NULL)

{
  INIT_ERROR(error);
  if (A.dim() != B.dim() || A.dim() != C.dim()) {
    RAISE_ERROR(error, "Matrices not aligned.");
  }
  if (A.cublas_handle() != B.cublas_handle() || 
      A.cublas_handle() != C.cublas_handle()) {
    RAISE_ERROR(error, "CUBLAS handles do not match.");
  }
  gemm(op_A, op_B, A.dim(), alpha, A.data(), B.data(), beta, C.data(), error,
       A.cublas_handle());
  PASS_ERROR(error);
}


/*!
 * element wise multiplication
 */
void cu_elwise_mul(op_t op_A, op_t op_B, int dim, double *A, double *B,
		   double *C, int *error=NULL);


/*!
 * element wise multiplication
 */
template<typename T>
inline void elwise_mul(op_t op_A, op_t op_B, int dim, T *A, T *B, T *C,
		       int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(A)) {
#endif
    int size = dim*dim;
    if (op_A == op_B) {
      for (int i = 0; i < size; i++, A++, B++, C++) {
	*C = (*A) * (*B);
      }
    }
    else {
      for (int i = 0; i < dim; i++) {
	T *Aptr = &A[i];
	for (int j = 0; j < dim; j++, Aptr+=dim, B++, C++) {
	  *C = (*Aptr) * (*B);
	}
      }
    }
#ifdef HAVE_CUDA
  }
  else {
    cu_elwise_mul(op_A, op_B, dim, A, B, C, error);
    PASS_ERROR(error);
  }
#endif
}


/*!
 * element wise multiplication
 */
template<typename T>
inline void elwise_mul(op_t op_A, op_t op_B, mat<T> A, mat<T> B, mat<T> C,
		       int *error=NULL)
{
  INIT_ERROR(error);
  if (A.dim() != B.dim() || A.dim() != C.dim()) {
    RAISE_ERROR(error, "Matrices not aligned.");
  }
  elwise_mul(op_A, op_B, A.dim(), A.data(), B.data(), C.data(), error);
  PASS_ERROR(error);
}


/*!
 * inlined method, to do matrix-scalar multiplication; square matrix is
 * required and output vector is overwritten. CUDA device code.
 */
void cu_mat_mul_sca(int size, double alpha, double *A, double *B,
		    int *error=NULL);


/*!
 * inlined method, to do matrix-scalar multiplication; square matrix is
 * required and output vector is overwritten. Host and dispatch code.
 */
template<typename T>
inline void mat_mul_sca(int size, T alpha, T *A, T *B, int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(A)) {
#endif
    for (int i = 0; i < size; i++, A++, B++) {
      *B = alpha*(*A);
    }
#ifdef HAVE_CUDA
  }
  else {
    cu_mat_mul_sca(size, alpha, A, B, error);
    PASS_ERROR(error);
  }
#endif
}


/*!
 * inlined method, to do matrix-scalar multiplication; square matrix is
 * required and output vector is overwritten.
 */
template<typename T>
inline void mat_mul_sca(T alpha, mat<T> &A, mat<T> &B, int *error=NULL)
{
  INIT_ERROR(error);
  if (A.dim() != B.dim()) {
    RAISE_ERROR(error, "Matrices not aligned.");
  }
  mat_mul_sca(A.size(), alpha, A.data(), B.data(), error);
  PASS_ERROR(error);
}


/*!
 * Element wise multiplication and addition. CUDA device code.
 */
void cu_mat_mul_sca(int size, double alpha, double *A, double beta, double *B,
		    double *C, int *error=NULL);


/*!
 * Element wise multiplication and addition. Dispatch and host code.
 */
template<typename T>
inline void mat_mul_sca(int size, T alpha, T *A, T beta, T *B, T *C,
			int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  bool hA = cu_on_host(A), hB = cu_on_host(B), hC = cu_on_host(C);
  if (hA && hB && hC) {
#endif
    for (int i = 0; i < size; i++, A++, B++, C++) {
      *C = alpha*(*A) + beta*(*B);
    }
#ifdef HAVE_CUDA
  }
  else if (!hA && !hB && !hC) {
    cu_mat_mul_sca(size, alpha, A, beta, B, C, error);
    PASS_ERROR(error);
  }
  else {
    RAISE_ERROR(error, "Mixed device/host *mat_mul_sca* not yet supported.");
  }
#endif
}


/*!
 * element wise multiplication
 */
template<typename T>
inline void mat_mul_sca(T alpha, mat<T> &A, T beta, mat<T> &B, mat<T> &C,
			int *error=NULL)
{
  INIT_ERROR(error);
  if (A.dim() != B.dim() || A.dim() != C.dim()) {
    RAISE_ERROR(error, "Matrices not aligned.");
  }
  mat_mul_sca(A.size(), alpha, A.data(), beta, B.data(), C.data(), error);
  PASS_ERROR(error);
}


/*!
 * Element wise multiplication and addition. CUDA device code.
 */
void cu_mat_mul_sca(int size, double alpha, double *A, double beta, double *B,
		    double gamma, double *C, double *D, int *error=NULL);


/*!
 * element wise multiplication and addition
 */
template<typename T>
inline void mat_mul_sca(int size, T alpha, T *A, T beta, T *B, T gamma,
			T *C, T *D, int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  bool hA = cu_on_host(A), hB = cu_on_host(B), hC = cu_on_host(C);
  bool hD = cu_on_host(D);
  if (hA && hB && hC && hD) {
#endif
    for (int i = 0; i < size; i++, A++, B++, C++, D++) {
      *D = alpha*(*A) + beta*(*B) + gamma*(*C);
    }
#ifdef HAVE_CUDA
  }
  else if (!hA && !hB && !hC && !hD) {
    cu_mat_mul_sca(size, alpha, A, beta, B, gamma, C, D, error);
    PASS_ERROR(error);
  }
  else {
    RAISE_ERROR(error, "Mixed device/host *mat_mul_sca* not yet supported.");
  }
#endif
}


/*!
 * element wise multiplication
 */
template<typename T>
inline void mat_mul_sca(T alpha, mat<T> &A, T beta, mat<T> &B, T gamma,
			mat<T> &C, mat<T> &D, int *error=NULL)
{
  INIT_ERROR(error);
  if (A.dim() != B.dim() || A.dim() != C.dim() || A.dim() != D.dim()) {
    RAISE_ERROR(error, "Matrices not aligned.");
  }
  mat_mul_sca(A.size(), alpha, A.data(), beta, B.data(), gamma, C.data(),
	      D.data(), error);
  PASS_ERROR(error);
}


/*!
 * Eigenvalue bounds
 *
 * Determine (conservative) lower and upper bounds for the
 * eigenvalue spectrum of a matrix. Host code.
 */
template<typename T>
inline void host_ev_bounds(int n, T *H, double *l, double *u, int *error=NULL)
{
  double lb, ub;

  for (int i = 0; i < n; i++) {
    double lh = std::abs(H[_IDX2(n, i, i)]);
    double uh = 0.0;

    for (int j = 0; j < n-1; j++) {
      if (i != j) {
	lh = lh - std::abs(H[_IDX2(n, i, j)]);
	uh = uh + std::abs(H[_IDX2(n, i, j)]);
      }
    }
          
    if (i == 0) {
      lb = lh;
      ub = uh;
    }
    else {
      lb = MIN(lb, lh);
      ub = MAX(ub, uh);
    }
  }

  *l = lb;
  *u = ub;
}


/*!
 * Eigenvalue bounds
 *
 * Determine (conservative) lower and upper bounds for the
 * eigenvalue spectrum of a matrix. CUDA device code.
 */
void cu_ev_bounds(int n, double *H, double *l, double *u, int *error=NULL);


/*!
 * Eigenvalue bounds
 *
 * Determine (conservative) lower and upper bounds for the
 * eigenvalue spectrum of a matrix. Dispatch code.
 */
inline void ev_bounds(int n, double *H, double *l, double *u, int *error=NULL)
{
  INIT_ERROR(error);
#ifdef HAVE_CUDA
  if (cu_on_host(H)) {
#endif
    host_ev_bounds(n, H, l, u, error);
    PASS_ERROR(error);
#ifdef HAVE_CUDA
  }
  else {
    cu_ev_bounds(n, H, l, u, error);
    PASS_ERROR(error);
  }
#endif
}

#endif
