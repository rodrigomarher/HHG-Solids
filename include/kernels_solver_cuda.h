#pragma once

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include "cuda_helpers/helper_cuda.h"

#define THREADS_CUDA 256
#define cdouble_cuda cuDoubleComplex 

__global__ void kernel_calculate_peierls_phase(double ax, double ay, double az, cdouble_cuda* d_peierls_phase, double* d_r_vec_x, double* d_r_vec_y, double* d_r_vec_z, const int num_points);

__global__ void kernel_reorder_and_peierls_k1(cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, 
                                           cdouble_cuda* d_heff, cdouble_cuda* d_rho,
                                           cdouble_cuda* d_peierls, const int num_points, const int num_orbitals); 

__global__ void kernel_reorder_and_peierls_kn(cdouble_cuda prefac, cdouble_cuda* kn,cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, 
                                              cdouble_cuda* d_heff, cdouble_cuda* d_rho,
                                              cdouble_cuda* d_peierls, const int num_points, const int num_orbitals); 
__global__ void kernel_fftshift2D(cdouble_cuda* data, const int nr1, const int nr2, const int num_orbitals);
__global__ void kernel_scale_inverse_fft(cdouble_cuda* data, const int num_points, const int num_orbitals);
__global__ void kernel_update_heff(const double ex, const double ey, const double ez, cdouble_cuda* d_heff, cdouble_cuda* d_h0, cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc, const int num_points, const int num_orbitals);
__global__ void kernel_commutator_k(cdouble_cuda* d_comm_k, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, const int num_points, const int num_orbitals);
__global__ void kernel_update_kn(cdouble_cuda* d_kn, cdouble_cuda* d_comm, const int num_points, const int num_orbitals);
__global__ void kernel_step_rho(cdouble_cuda* d_rho, cdouble_cuda* d_k1, cdouble_cuda* d_k2, cdouble_cuda* d_k3, cdouble_cuda* d_k4, const double dt, const int num_points, const int num_orbitals);

void call_kernel_calculate_peierls_phase(double ax, double ay, double az, cdouble_cuda* d_peierls_phase, double* d_r_vec_x, double* d_r_vec_y, double* d_r_vec_z, const int num_points, cudaStream_t* stream);

void call_kernel_update_heff(const double ex, const double ey, const double ez, cdouble_cuda* d_heff, cdouble_cuda* d_h0, cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc, const int num_points, const int num_orbitals, cudaStream_t* stream);

void call_pipeline_update_k1(cdouble_cuda* d_k1, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan);
void call_pipeline_update_k2(cdouble_cuda* d_k2, cdouble_cuda* d_k1, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan);
void call_pipeline_update_k3(cdouble_cuda* d_k3, cdouble_cuda* d_k2, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan);
void call_pipeline_update_k4(cdouble_cuda* d_k4, cdouble_cuda* d_k3, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan);
void call_kernel_step_rho(cdouble_cuda* d_rho, cdouble_cuda* d_k1, cdouble_cuda* d_k2, cdouble_cuda* d_k3, cdouble_cuda* d_k4, const double dt, const int num_points, const int num_orbitals, cudaStream_t* stream);
#endif 
