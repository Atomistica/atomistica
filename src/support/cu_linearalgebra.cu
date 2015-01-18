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

#include "cu_vec.h"
#include "linearalgebra.h"


/*!
 * Element wise multiplication. CUDA kernel.
 */
__global__ void _elwise_mul(int size, double *A, double *B, double *C)
{
  int i = (blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;

  if (i < size) {
    C[i] = A[i]*B[i];
  }
}


/*!
 * Element wise multiplication with transpose. CUDA kernel.
 */
__global__ void _elwise_mul_T(int dim, int size, double *A, double *B,
			      double *C)
{
  int k = (blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;

  if (k < size) {
    int i = k % dim;
    int j = k / dim;

    C[k] = A[i*dim+j]*B[k];
  }
}


/*!
 * Element wise multiplication
 */
void cu_elwise_mul(op_t op_A, op_t op_B, int dim, double *A, double *B,
		   double *C, int *error)
{
  INIT_ERROR;
  int size = dim*dim;
  int num_blocks = size / THREADS_PER_BLOCK;
  if (num_blocks*THREADS_PER_BLOCK < size)  num_blocks++;

  if (op_A == op_B) {
    _elwise_mul<<<num_blocks, THREADS_PER_BLOCK>>>(size, A, B, C);
    PASS_CUDA_ERROR;
  }
  else {
    _elwise_mul_T<<<num_blocks, THREADS_PER_BLOCK>>>(dim, size, A, B, C);
    PASS_CUDA_ERROR;
  }
}


/*!
 * Multiply matrix with a scalar
 */
__global__ void _mat_mul_sca(int size, double alpha, double *A, double *B)
{
  int i = (blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;

  if (i < size) {
    B[i] = alpha*A[i];
  }
}


/*!
 * Multiply matrix with a scalar
 */
void cu_mat_mul_sca(int size, double alpha, double *A, double *B, int *error)
{
  INIT_ERROR;
  int num_blocks = size / THREADS_PER_BLOCK;
  if (num_blocks*THREADS_PER_BLOCK < size)  num_blocks++;

  _mat_mul_sca<<<num_blocks, THREADS_PER_BLOCK>>>(size, alpha, A, B);
  PASS_CUDA_ERROR;
}


/*!
 * Multiply matrix with a scalar and add elements
 */
__global__ void _mat_mul_sca(int size, double alpha, double *A, double beta,
			     double *B, double *C)
{
  int i = (blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;

  if (i < size) {
    C[i] = alpha*A[i] + beta*B[i];
  }
}


/*!
 * Multiply matrix with a scalar and add elements
 */
void cu_mat_mul_sca(int size, double alpha, double *A, double beta, double *B,
		    double *C, int *error)
{
  INIT_ERROR;
  int num_blocks = size / THREADS_PER_BLOCK;
  if (num_blocks*THREADS_PER_BLOCK < size)  num_blocks++;

  _mat_mul_sca<<<num_blocks, THREADS_PER_BLOCK>>>(size, alpha, A, beta, B, C);
  PASS_CUDA_ERROR;
}


/*!
 * Multiply matrix with a scalar and add elements
 */
__global__ void _mat_mul_sca(int size, double alpha, double *A, double beta,
			     double *B, double gamma, double *C, double *D)
{
  int i = (blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;

  if (i < size) {
    C[i] = alpha*A[i] + beta*B[i] + gamma*D[i];
  }
}


/*!
 * Multiply matrix with a scalar and add elements
 */
void cu_mat_mul_sca(int size, double alpha, double *A, double beta, double *B,
		    double gamma, double *C, double *D, int *error)
{
  INIT_ERROR;
  int num_blocks = size / THREADS_PER_BLOCK;
  if (num_blocks*THREADS_PER_BLOCK < size)  num_blocks++;

  _mat_mul_sca<<<num_blocks, THREADS_PER_BLOCK>>>(size, alpha, A, beta, B,
						  gamma, C, D);
  PASS_CUDA_ERROR;
}


/*
 * FIXME: This is a hack. Copies matrix to host and then does the eigenvalue
 * bounds. Need to move this to the GPU. However, this is not the most
 * time consuming part of the calculation.
 */
void cu_ev_bounds(int n, double *H, double *l, double *u, int *error)
{
  INIT_ERROR;

  mat<double> Htmp(n);

  /* Copy from device to host */
  Htmp = H;

  host_ev_bounds(n, Htmp.data(), l, u, error);
  PASS_ERROR;
}
