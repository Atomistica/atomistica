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
#ifndef __MAT_H
#define __MAT_H

#include "error.h"
#include "complexcomp.h"
#include "cu_mat.h"
#include "cu_vec.h"

#define _IDX2(dim, i, j)  ( (i)*(dim)+(j) )

/*
 * Simple matrix class, including support for addition and multiplication
 */

template<typename T>
class mat {
 public:
  mat(int dim, T *data = NULL, cublasHandle_t cublas_handle = NULL) {
    dim_ = dim;
    size_ = dim*dim;
    data_ = data;
    own_data_ = false;
    if (!data_) {
#ifdef HAVE_CUDA
      if (cublas_handle) {
	cudaMalloc(&data_, size_*sizeof(T));
	cudaMemset(data_, 0, size_*sizeof(T));
      }
      else {
#endif
	data_ = new T[size_];
	memset(data_, 0, size_*sizeof(T));
#ifdef HAVE_CUDA
      };
#endif
      own_data_ = true;
    }
    cublas_handle_ = cublas_handle;
  }

  ~mat() {
    if (own_data_) {
#ifdef HAVE_CUDA
      if (on_host()) {
#endif
	delete [] data_;
#ifdef HAVE_CUDA
      }
      else {
	cudaFree(data_);
      }
#endif
    }
  }

  /*
   * Data pointer, matrix dimension and size
   */

  T *data() {
    return data_;
  }

  const T *data() const {
    return data_;
  }

  int dim() const {
    return dim_;
  }  

  int size() const {
    return size_;
  }

  /*
   * Assignment and element access operators
   */

  void fill_with(T value) {
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      for (int i = 0; i < size_; i++){
	data_[i] = value;
      }
#ifdef HAVE_CUDA
    }
    else {
      cu_fill_with(data_, value, size_);
    }
#endif
  }

  mat<T> &operator=(const mat<T> &other) {
#ifdef HAVE_CUDA
    cudaMemcpy(data_, other.data_, size_*sizeof(T), cudaMemcpyDefault);
#else
    memcpy(data_, other.data_, size_*sizeof(T));
#endif
    return *this;
  }

  mat<T> &operator=(T *other) {
#ifdef HAVE_CUDA
    cudaMemcpy(data_, other, size_*sizeof(T), cudaMemcpyDefault);
#else
    memcpy(data_, other, size_*sizeof(T));
#endif
    return *this;
  }

  mat<T> &operator=(T data) {
    fill_with(data);
    return *this;
  }

  mat<T> &operator+=(T v) {
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      T *ptr = data_;
      for (int i = 0; i < dim_; i++, ptr+=dim_+1) {
	*ptr += v;
      }
#ifdef HAVE_CUDA
    }
    else {
      cu_add_to_diagonal(dim_, data_, v);
    }
#endif
    return *this;
  }

  bool almost_equal(const mat<T> &other, T tol=1e-12, int *error=NULL) const {
    INIT_ERROR(error);
    bool equal = true;
    if (size_ != other.size_)  return false;
#ifdef HAVE_CUDA
    bool ht = on_host(), ho = other.on_host();
    if (ht && ho) {
#endif
      T *ptrt = data_, *ptro = other.data_;
      for (int i = 0; i < size_ && equal; i++, ptrt++, ptro++) {
	if (abs(*ptrt - *ptro) > tol)  equal = false;
      }
#ifdef HAVE_CUDA
    }
    else {
      RAISE_ERROR_WITH_RET(error, false,
			   "On-device *almost_equal* not yet supported.");
    }
#endif
    return equal;
  }

  T *operator[](int x) {
    return &data_[_IDX2(dim_, x, 0)];
  }

  T operator()(int x, int y, int *error=NULL) {
    INIT_ERROR(error);
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      return data_[_IDX2(dim_, x, y)];
#ifdef HAVE_CUDA
    }
    else {
      T v;
      cudaMemcpy(&v, &data_[_IDX2(dim_, x, y)], sizeof(T),
		 cudaMemcpyDeviceToHost);
      PASS_CUDA_ERROR_WITH_RET(error, 0.0);
      return v;
    }
#endif
  }

  void set(int x, int y, T v, int *error=NULL) {
    INIT_ERROR(error);
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      data_[_IDX2(dim_, x, y)] = v;
#ifdef HAVE_CUDA
    }
    else {
      cudaMemcpy(&data_[_IDX2(dim_, x, y)], &v, sizeof(T),
		 cudaMemcpyHostToDevice);
      PASS_CUDA_ERROR(error);
    }
#endif
  }

  /*
   * Wrappers for CUBLAS functions
   */

  void axpy(double alpha, const double *A, int *error=NULL) {
    INIT_ERROR(error);
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      for (int i = 0; i < size_; i++) {
	data_[i] += alpha*A[i];
      }
#ifdef HAVE_CUDA
    }
    else {
      PASS_CUBLAS_ERROR( error,
			 cublasDaxpy(cublas_handle_, size_, &alpha, A, 1,
				     data_, 1) );
    }
#endif
  }

  void axpy(double_complex alpha, const double_complex *A, int *error=NULL) {
    INIT_ERROR(error);
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      for (int i = 0; i < size_; i++) {
	data_[i] += alpha*A[i];
      }
#ifdef HAVE_CUDA
    }
    else {
      PASS_CUBLAS_ERROR( error,
			 cublasZaxpy(cublas_handle_, size_, &alpha, A, 1,
				     data_, 1) );
    }
#endif
  }

  void axpy(T alpha, const mat<T> &A, int *error=NULL) {
    INIT_ERROR(error);
    if (A.dim() != dim_) {
      RAISE_ERROR(error, "Matrices not aligned");
    }
    axpy(alpha, A.data(), error);
    PASS_ERROR(error);
  }

  /*
   * In-place operators
   */

  void operator+=(const mat<T> &A) {
    axpy(1.0, A);
  }

  void operator+=(T *A) {
    axpy(1.0, A);
  }

  void operator-=(const mat<T> &A) {
    axpy(-1.0, A);
  }

  void operator-=(T *A) {
    axpy(-1.0, A);
  }

  /*!
   * Return maximum of absolute value over all matrix elements
   */
  T amax(int *error=NULL) {
    INIT_ERROR(error);
    T v;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      T *p = data_;
      v = fabs(*p);
      for (int i = 0; i < size_; i++, p++) {
	v = MAX(v, fabs(*p));
      }
#ifdef HAVE_CUDA
    }
    else {
      v = cu_amax(size_, data_, error);
      PASS_ERROR_WITH_RET(error, 0.0);
    }
#endif
    return v;
  }

  /*!
   * Return minium of absolute value over all matrix elements
   */
  T amin(int *error=NULL) {
    INIT_ERROR(error);
    T v;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      T *p = data_;
      v = fabs(*p);
      for (int i = 0; i < size_; i++, p++) {
	v = MIN(v, fabs(*p));
      }
#ifdef HAVE_CUDA
    }
    else {
      v = cu_amin(size_, data_, error);
      PASS_ERROR_WITH_RET(error, 0.0);
    }
#endif
    return v;
  }

  /*!
   * Return maximum value over all matrix elements
   */
  T max(int *error=NULL) {
    INIT_ERROR(error);
    T v;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      T *p = data_;
      v = abs(*p);
      for (int i = 0; i < size_; i++, p++) {
	v = MAX(v, *p);
      }
#ifdef HAVE_CUDA
    }
    else {
      v = cu_max(size_, data_, error);
      PASS_ERROR_WITH_RET(error, 0.0);
    }
#endif
    return v;
  }

  /*!
   * Return minium value over all matrix elements
   */
  T min(int *error=NULL) {
    INIT_ERROR(error);
    T v;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      T *p = data_;
      v = abs(*p);
      for (int i = 0; i < size_; i++, p++) {
	v = MIN(v, *p);
      }
#ifdef HAVE_CUDA
    }
    else {
      v = cu_min(size_, data_, error);
      PASS_ERROR_WITH_RET(error, 0.0);
    }
#endif
    return v;
  }

  double nrm2(int *error=NULL) {
    INIT_ERROR(error);
    double n;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      double acc = 0.0;
      for (int i = 0; i < size_; i++) {
	/* acc += creal(conj(data_[i])*data_[i]); */
	acc += data_[i]*data_[i];
      }
      n = sqrt(acc);
#ifdef HAVE_CUDA
    }
    else {
      PASS_CUBLAS_ERROR_WITH_RET( error,
				  0.0, cublasDnrm2(cublas_handle_, size_, data_,
						   1, &n) );
    }
#endif
    return n;
  }

  double nrm(int ord=2, int *error=NULL) {
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      double acc = 0.0;
      for (int i = 0; i < size_; i++) {
	acc += pow(data_[i], ord);
      }
      return pow(acc, 1.0/ord);
#ifdef HAVE_CUDA
    }
    else {
      RAISE_ERROR_WITH_RET(error, 0.0, "On-device *nrm* not yet supported.");
    }
#endif
  }

  /*!
   * Return sum over all matrix elements
   */
  T sum(int *error=NULL) {
    T v = 0.0;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      T *p = data_;
      for (int i = 0; i < size_; i++, p++) {
	v += *p;
      }
#ifdef HAVE_CUDA
    }
    else {
      v = cu_sum(size_, data_, error);
      PASS_ERROR_WITH_RET(error, 0.0);
    }
#endif
    return v;
  }

  /*!
   * Return traces of the matrix
   */
  T trace(int *error=NULL) {
    T tr = 0.0;
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      int m = 0;
      for (int i = 0; i < dim_; i++, m+=dim_+1) {
	tr += data_[m];
      }
#ifdef HAVE_CUDA
    }
    else {
      tr = cu_trace(dim_, data_, error);
      PASS_ERROR_WITH_RET(error, 0.0);
    }
#endif
    return tr;
  }

  /*!
   * In-place transpose operation
   */
  void transpose(int *error=NULL) {
#ifdef HAVE_CUDA
    if (on_host()) {
#endif
      for (int i = 1; i < dim_; i++){
	int m = i*dim_;
	int n = i;
	for (int j = 0; j < i; j++, m++, n+=dim_) {
	  T tmp = data_[m];
	  data_[m] = data_[n];
	  data_[n] = tmp;
	}
      }
#ifdef HAVE_CUDA
    }
    else {
      RAISE_ERROR(error, "On-device *transpose* not yet supported.");
    }
#endif
  }

#ifdef HAVE_CUDA
  /*!
   * Check whether this matrix resides on the device
   */
  bool on_device(int *error=NULL) const {
    return cu_on_device(data_, error);
  }

  /*!
   * Check whether this matrix resides on the device
   */
  bool on_host(int *error=NULL) const {
    return cu_on_host(data_, error);
  }
#endif

  /*!
   * Return the cublas handle
   */
  cublasHandle_t cublas_handle() const {
    return cublas_handle_;
  }


 protected:
  int dim_, size_;
  T *data_;
  bool own_data_;

  cublasHandle_t cublas_handle_;
};

#endif
