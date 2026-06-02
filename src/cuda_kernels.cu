#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"
#include "cuda_kernels.h"

__global__ void kernel_update_heff(const double ex, const double ey, const double ez,
                                   cdouble_cuda* d_heff, cdouble_cuda* d_h0, 
                                   cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc,
                                   const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx_r >= num_points){return;}
    
    cdouble_cuda value = make_cuDoubleComplex(0.0,0.0);
    for(int iorb = 0; iorb<num_orbitals*num_orbitals; iorb++){
        value = d_h0[idx_r*num_orbitals*num_orbitals + iorb];
        value = cuCadd(value, cuCmul(make_cuDoubleComplex(ex,0.0), d_xbc[idx_r*num_orbitals*num_orbitals + iorb]));
        value = cuCadd(value, cuCmul(make_cuDoubleComplex(ey,0.0), d_ybc[idx_r*num_orbitals*num_orbitals + iorb]));
        value = cuCadd(value, cuCmul(make_cuDoubleComplex(ez,0.0), d_zbc[idx_r*num_orbitals*num_orbitals + iorb]));
        d_heff[idx_r*num_orbitals*num_orbitals + iorb] = value;
    }
}

void call_kernel_update_heff(const double ex, const double ey, const double ez, cdouble_cuda* d_heff, cdouble_cuda* d_h0, cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc, const int num_points, const int num_orbitals, cudaStream_t* stream){
    int threadsPerBlock = THREADS_CUDA;
    int blocks = (num_points + threadsPerBlock - 1)/threadsPerBlock; 
    kernel_update_heff<<<blocks, threadsPerBlock, 0, *stream>>>(ex, ey, ez, d_heff, d_h0, d_xbc, d_ybc, d_zbc, num_points, num_orbitals);
}
