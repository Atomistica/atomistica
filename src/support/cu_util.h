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
#ifndef __CU_UTIL_H
#define __CU_UTIL_H

#include <stdio.h>
#include <string.h>

/*
 * Minimum and maximum
 */
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

#ifdef HAVE_CUDA

extern "C" {
#include "libatoms.h"
}

#include "cuda_runtime.h"
#include "cublas_v2.h"
#define PASS_CUDA_ERROR(error) { cudaError_t e = cudaGetLastError(); if (e != cudaSuccess) { RAISE_ERROR(error, cudaGetErrorString(e)); } }
#define PASS_CUDA_ERROR_WITH_RET(error, r) { cudaError_t e = cudaGetLastError(); if (e != cudaSuccess) { RAISE_ERROR_WITH_RET(error, r, cudaGetErrorString(e)); } }
#define PASS_CUBLAS_ERROR(error, x) { cublasStatus_t e = x; if (e != CUBLAS_STATUS_SUCCESS) { RAISE_ERROR(error, get_cublas_error_string(e)); } }
#define PASS_CUBLAS_ERROR_WITH_RET(error, r, x) { if (x != CUBLAS_STATUS_SUCCESS) { RAISE_ERROR_WITH_RET(error, r, "Call to CUBLAS routine failed."); } }

/*!
 * Check the CUDA memory type for this pointer
 */
inline cudaMemoryType cu_memory_type(void *data, int *error=NULL) {
  INIT_ERROR(error);
  cudaPointerAttributes attr;
  cudaError_t err = cudaPointerGetAttributes(&attr, data);
  /* cudaPointerGetAttributes seems to generate an error for malloc'd ptrs */
  if (err == cudaErrorInvalidValue) {
    /* Clear error state */
    cudaGetLastError();
    return cudaMemoryTypeHost;
  }
  else {
    PASS_CUDA_ERROR_WITH_RET(error, cudaMemoryTypeHost);
  }
  return attr.memoryType;
}

/*!
 * Check whether this matrix resides on the device
 */
inline bool cu_on_device(void *data, int *error=NULL) {
  return cu_memory_type(data, error) == cudaMemoryTypeDevice;
}

/*!
 * Check whether this matrix resides on the device
 */
inline bool cu_on_host(void *data, int *error=NULL) {
  return cu_memory_type(data, error) == cudaMemoryTypeHost;
}

const char *get_cublas_error_string(cublasStatus_t err);

#define THREADS_PER_BLOCK 512

template<typename T>
inline void destroy(T *&ptr) {
  if (ptr) {
    if (cu_on_host(ptr)) {
      free(ptr);
    }
    else {
      cudaFree(ptr);
    }
    ptr = NULL;
  }
}

#else

typedef void *cublasHandle_t;

template<typename T>
inline void destroy(T *&ptr) {
  if (ptr) {
    free(ptr);
    ptr = NULL;
  }
}

#endif

#endif
