#pragma once

#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"

#define THREADS_CUDA 256
#define cdouble_cuda cuDoubleComplex

__global__ void kernel_calc_observable(cdouble_cuda* op, cdouble_cuda* rho, cdouble_cuda* data, double* r_x, double* r_y, double* r_z, double ax, double ay, double az, const int ti, const int num_points, const int num_orbitals);

__global__ void kernel_calc_observable_shared(cdouble_cuda* op, cdouble_cuda* rho, cdouble_cuda* data, double* r_x, double* r_y, double* r_z, double ax, double ay, double az, const int ti, const int num_points, const int num_orbitals);

void call_kernel_calc_observable(cdouble_cuda *op, cdouble_cuda* rho, cdouble_cuda* data, double* r_x, double* r_y, double* r_z, double ax, double ay, double az, const int ti, const int num_points, const int num_orbitals, cudaStream_t* stream);

#endif
