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
#ifndef __VEC_H
#define __VEC_H

/*
 * Simple vector class, including support for addition and multiplication
 */

template<typename T>
class vec {
 public:
  vec(int dim, T *data = NULL) {
    dim_ = dim;
    data_ = data;
    own_data_ = false;
    if (!data_) {
      data_ = new T[dim_];
      own_data_ = true;
      memset(data_, 0, dim_*sizeof(T));
    }
  }

  ~vec() {
    if (own_data_) {
      delete [] data_;
    }
  }

  void operator=(const vec<T> &mat) {
    memcpy(data_, mat.data_, dim_*sizeof(T));
  }

  void operator=(T *data) {
    memcpy(data_, data, dim_*sizeof(T));
  }

  T &operator[](int x) {
    return data_[x];
  }

  void operator+=(const vec<T> &A) {
    for (int i = 0; i < dim_; i++) {
      data_[i] += A.data_[i];
    }
  }

  void operator+=(T *A) {
    for (int i = 0; i < dim_; i++) {
      data_[i] += A[i];
    }
  }

  void operator-=(const vec<T> &A) {
    for (int i = 0; i < dim_; i++) {
      data_[i] -= A.data_[i];
    }
  }

  void operator-=(T *A) {
    for (int i = 0; i < dim_; i++) {
      data_[i] -= A[i];
    }
  }

  void fill_with(T value) {
    for (int i = 0; i < dim_; i++){
      data_[i] = value;
    }
  }

  T *data() {
    return data_;
  }

  double norm2() {
    double acc = 0.0;
    for (int i = 0; i < dim_; i++) {
      acc += creal(conj(data_[i])*data_[i]);
    }
    return sqrt(acc);
  }

  double norm(int ord=2) {
    double acc = 0.0;
    for (int i = 0; i < dim_; i++) {
      acc += pow(data_[i], ord);
    }
    return pow(acc, 1.0/ord);
  }

  int dim_;
  T *data_;
  bool own_data_;
};

#endif
