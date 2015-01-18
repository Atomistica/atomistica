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

#include "cu_vec.h"

const char cublas_error_strings[][80] = {
  "Unknown CUBLAS error.",
  "CUBLAS_STATUS_NOT_INITIALIZED",
  "CUBLAS_STATUS_ALLOC_FAILED",
  "CUBLAS_STATUS_INVALID_VALUE",
  "CUBLAS_STATUS_ARCH_MISMATCH",
  "CUBLAS_STATUS_MAPPING_ERROR",
  "CUBLAS_STATUS_EXECUTION_FAILED",
  "CUBLAS_STATUS_INTERNAL_ERROR"
};

const char *get_cublas_error_string(cublasStatus_t err)
{
  switch (err) {
  case CUBLAS_STATUS_NOT_INITIALIZED:
    return cublas_error_strings[1];
  case CUBLAS_STATUS_ALLOC_FAILED:
    return cublas_error_strings[2];
  case CUBLAS_STATUS_INVALID_VALUE:
    return cublas_error_strings[3];
  case CUBLAS_STATUS_ARCH_MISMATCH:
    return cublas_error_strings[4];
  case CUBLAS_STATUS_MAPPING_ERROR:
    return cublas_error_strings[5];
  case CUBLAS_STATUS_EXECUTION_FAILED:
    return cublas_error_strings[6];
  case CUBLAS_STATUS_INTERNAL_ERROR:
    return cublas_error_strings[7];
  }

  return cublas_error_strings[0];
}

/*!
 * Reduction kernel. Template arguments the number of threads per block,
 * input type (*itype*), output type (*otype*) and the reduction and map
 * operations.
 *
 * This is taken from the NVIDIA CUDA reduction slides.
 *
 * Reduction and map operations should be defined as functors (see below):
 *   struct sum { double operator()(double a, double b) { return a+b }; }
 */
template<int threads_per_block, typename itype, typename otype,
	 typename reduce_type, typename map_type>
__global__ void _reduction(int n, itype *in, otype *out, reduce_type reduce_op,
			   map_type map_op)
{								       
  extern __shared__ otype sdata[];

  int tid = threadIdx.x;
  int i = blockIdx.x*(threads_per_block*2) + tid;
  int grid_size = threads_per_block*2*gridDim.x;

  otype acc = reduce_op.neutral(map_op(n, in, i));

  while (i < n) {
    acc = reduce_op(acc, map_op(n, in, i));
    if (i+threads_per_block < n) {
      acc = reduce_op(acc, map_op(n, in, i+threads_per_block));
    }
    i += grid_size;
  }

  sdata[tid] = acc;

  __syncthreads();

  if (threads_per_block >= 512) {
    if (tid < 256) { sdata[tid] = reduce_op(sdata[tid], sdata[tid + 256]); }
    __syncthreads();
  }

  if (threads_per_block >= 256) {
    if (tid < 128) { sdata[tid] = reduce_op(sdata[tid], sdata[tid + 128]); }
    __syncthreads();
  }

  if (threads_per_block >= 128) {
    if (tid < 64) { sdata[tid] = reduce_op(sdata[tid], sdata[tid + 64]); }
    __syncthreads();							
  }

  if (tid < 32) {
    volatile otype *smem = sdata;
    if (threads_per_block >= 64) smem[tid] = reduce_op(smem[tid], smem[tid+32]);
    if (threads_per_block >= 32) smem[tid] = reduce_op(smem[tid], smem[tid+16]);
    if (threads_per_block >= 16) smem[tid] = reduce_op(smem[tid], smem[tid+ 8]);
    if (threads_per_block >=  8) smem[tid] = reduce_op(smem[tid], smem[tid+ 4]);
    if (threads_per_block >=  4) smem[tid] = reduce_op(smem[tid], smem[tid+ 2]);
    if (threads_per_block >=  2) smem[tid] = reduce_op(smem[tid], smem[tid+ 1]);
  }

  if (tid == 0) out[blockIdx.x] = sdata[0];
}


#define MAX_BLOCKS_PER_GRID 1024


/*!
 * Launcher for the reduction kernel
 */
template<typename itype, typename otype, typename reduce_type,
	 typename map_type>
otype cu_reduction(int n, itype *in, reduce_type reduce_op, map_type map_op,
		   int *error, int threads_per_block)
{
  INIT_ERROR;

  /*
   * If threads_per_block < n, then we need to reduce the number of 
   * threads_per_block
   */

  while (threads_per_block > n)  threads_per_block /= 2;

  /*
   * Compute grid and block size
   */

  int blocks_per_grid = n / threads_per_block;
  if (blocks_per_grid*threads_per_block < n)  blocks_per_grid++;
  blocks_per_grid = MIN(blocks_per_grid, MAX_BLOCKS_PER_GRID);
  int smem = threads_per_block*sizeof(otype);

  /*
   * Output buffer on device
   */

  otype *out;
  
  cudaMalloc(&out, blocks_per_grid*sizeof(otype));
  PASS_CUDA_ERROR_WITH_RET(0.0);

  /*
   * Launch kernel
   */

  switch (threads_per_block) {
  case 512:
    _reduction<512, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case 256:
    _reduction<256, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case 128:
    _reduction<128, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case 64:
    _reduction< 64, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case 32:
    _reduction< 32, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case 16:
    _reduction< 16, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case  8:
    _reduction<  8, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case  4:
    _reduction<  4, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case  2:
    _reduction<  2, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  case  1:
    _reduction<  1, itype, otype, reduce_type, map_type>
      <<<blocks_per_grid,threads_per_block,smem>>>(n, in, out, reduce_op,
						   map_op);
    break;
  }
  PASS_CUDA_ERROR_WITH_RET(0.0);

  /*
   * Final reduction on CPU
   */

  otype hout[blocks_per_grid];
  cudaMemcpy(hout, out, blocks_per_grid*sizeof(otype), cudaMemcpyDeviceToHost);
  PASS_CUDA_ERROR_WITH_RET(0.0);
  cudaFree(out);
  PASS_CUDA_ERROR_WITH_RET(0.0);

  double acc = hout[0];
  for (int i = 1; i < blocks_per_grid; i++) {
    acc = reduce_op(acc, hout[i]);
  }

  return acc;
}


/*
 * Specific reduction operations
 */

template<typename T>
struct _identity {
  __device__ __host__ T operator()(int n, T *in, int i) const {
    return in[i];
  }
};

template<typename T>
struct _abs {
  __device__ __host__ T operator()(int n, T *in, int i) const {
    return fabs(in[i]);
  }
};

template<typename T>
struct _diagonal {
  __device__ __host__ T operator()(int n, T *in, int i) const {
    return in[i*(n+1)];
  }
};

template<typename T>
struct _sum {
  __device__ __host__ T operator()(T a, T b) const {
    return a+b;
  }
  __device__ __host__ T neutral(T a) const {
    return 0.0;
  }
};

template<typename T>
struct _max {
  __device__ __host__ T operator()(T a, T b) const {
    return MAX(a, b);
  }
  __device__ __host__ T neutral(T a) const {
    return a;
  }
};

template<typename T>
struct _min {
  __device__ __host__ T operator()(T a, T b) const {
    return MIN(a, b);
  }
  __device__ __host__ T neutral(T a) const {
    return a;
  }
};


/*!
 * Summation
 */
template<typename T>
T _cu_sum(int n, T *A, int *error, int threads_per_block)
{
  return cu_reduction< T, T, _sum<T>, _identity<T> >
    (n, A, _sum<T>(), _identity<T>(), error, threads_per_block);
}

double cu_sum(int n, double *A, int *error, int threads_per_block)
{
  return _cu_sum(n, A, error, threads_per_block);
}


/*!
 * Maximum value
 */
template<typename T>
T _cu_max(int n, T *A, int *error, int threads_per_block)
{
  return cu_reduction< T, T, _max<T>, _identity<T> >
    (n, A, _max<T>(), _identity<T>(), error, threads_per_block);
}

double cu_max(int n, double *A, int *error, int threads_per_block)
{
  return _cu_max(n, A, error, threads_per_block);
}


/*!
 * Minimum value
 */
template<typename T>
T _cu_min(int n, T *A, int *error, int threads_per_block)
{
  return cu_reduction< T, T, _min<T>, _identity<T> >
    (n, A, _min<T>(), _identity<T>(), error, threads_per_block);
}

double cu_min(int n, double *A, int *error, int threads_per_block)
{
  return _cu_min(n, A, error, threads_per_block);
}


/*!
 * Absolute ,aximum value
 */
template<typename T>
T _cu_amax(int n, T *A, int *error, int threads_per_block)
{
  return cu_reduction< T, T, _max<T>, _abs<T> >
    (n, A, _max<T>(), _abs<T>(), error, threads_per_block);
}

double cu_amax(int n, double *A, int *error, int threads_per_block)
{
  return _cu_amax(n, A, error, threads_per_block);
}


/*!
 * Absolute minimum value
 */
template<typename T>
T _cu_amin(int n, T *A, int *error, int threads_per_block)
{
  return cu_reduction< T, T, _min<T>, _abs<T> >
    (n, A, _min<T>(), _abs<T>(), error, threads_per_block);
}

double cu_amin(int n, double *A, int *error, int threads_per_block)
{
  return _cu_amin(n, A, error, threads_per_block);
}


/*!
 * Trace
 */
template<typename T>
T _cu_trace(int n, T *A, int *error, int threads_per_block)
{
  return cu_reduction< T, T, _sum<T>, _diagonal<T> >
    (n, A, _sum<T>(), _diagonal<T>(), error, threads_per_block);
}

double cu_trace(int n, double *A, int *error, int threads_per_block)
{
  return _cu_trace(n, A, error, threads_per_block);
}
