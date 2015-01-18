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

#include "cu_mat.h"

/*!
 * Multiply matrix with a scalar and add elements
 */
template<typename T>
__global__ void _cu_add_to_diagonal(int dim, T *data, T v)
{
  int i = (blockIdx.y*gridDim.x + blockIdx.x)*blockDim.x + threadIdx.x;

  if (i < dim) {
    data[i*(dim+1)] += v;
  }
}


/*!
 * Multiply matrix with a scalar and add elements
 */
void cu_add_to_diagonal(int dim, double *data, double v, int *error)
{
  INIT_ERROR;
  int num_blocks = dim / THREADS_PER_BLOCK;
  if (num_blocks*THREADS_PER_BLOCK < dim)  num_blocks++;

  _cu_add_to_diagonal<<<num_blocks, THREADS_PER_BLOCK>>>(dim, data, v);
  PASS_CUDA_ERROR;
}


