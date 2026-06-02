#pragma once

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"

#define THREADS_CUDA 128
#define cdouble_cuda cuDoubleComplex 
__global__ void kernel_update_heff(const double ex, const double ey, const double ez, cdouble_cuda* d_heff, cdouble_cuda* d_h0, cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc, const int num_points, const int num_orbitals);

void call_kernel_update_heff(const double ex, const double ey, const double ez, cdouble_cuda* d_heff, cdouble_cuda* d_h0, cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc, const int num_points, const int num_orbitals, cudaStream_t* stream);
#endif 
