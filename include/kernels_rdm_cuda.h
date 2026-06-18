#pragma once

#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"

#define THREADS_CUDA 256
#define cdouble_cuda cuDoubleComplex

__global__ void kernel_normalize_n(cdouble_cuda *data, const int num_points, const int num_orbitals);

void call_kernel_normalize_n(cdouble_cuda *data, const int num_points, const int num_orbitals, cudaStream_t* stream);

#endif
