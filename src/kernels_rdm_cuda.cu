#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"
#include "kernels_rdm_cuda.h"

__global__ void kernel_rdm_normalize_n(cdouble_cuda* data, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x*blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if (idx_r > num_points || iorb>num_orbitals || jorb > num_points){return;}

    int idx = idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb;
    data[idx].x = data[idx].x/num_points;
    data[idx].y = data[idx].y/num_points;
}

void call_kernel_normalize_n(cdouble_cuda* data, const int num_points, const int num_orbitals, cudaStream_t* stream ){
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points+THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_rdm_normalize_n<<<grid, block, 0, *stream>>>(data, num_points, num_orbitals);
    }
}

#endif
