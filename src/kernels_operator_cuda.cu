#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"
#include "kernels_observable_cuda.h"

__global__ void kernel_calc_observable(cdouble_cuda* op, cdouble_cuda* rho, cdouble_cuda* data, const int ti, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x*blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;

    if(idx_r > num_points || iorb > num_orbitals || jorb > num_orbitals){return;}
    
    int idx = idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb;
    cdouble_cuda op_value = op[idx];
    cdouble_cuda rho_value = cuConj(rho[idx])
    cdouble_cuda value = cuCmul(op_value, rho_value);     
    atomicAdd(&data[ti], value.x);
    atomicAdd(&data[ti] + 1, value.y);
}

void call_kernel_calc_observable(cdouble_cuda* op, cdouble_cuda* rho, cdouble_cuda* data, const int ti, const int num_points, const int num_orbitals, cudaStream_t* stream){

    dim3 block(THREADS_CUDA);
    dim3 grid((num_points + THREADS_CUDA -1)/THREADS_CUDA, num_orbitals, num_orbitals);

    kernel_calc_observable<<<grid, bloc, 0, *stream>>>(op, rho, data, ti, num_points, num_orbitals);
}
#endif
