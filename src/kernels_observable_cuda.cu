#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"
#include "kernels_observable_cuda.h"

__global__ void kernel_calc_observable(cdouble_cuda* op, cdouble_cuda* rho, cdouble_cuda* data, double* r_x, double* r_y, double* r_z, double ax, double ay, double az, const int ti, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x*blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;

    if(idx_r > num_points || iorb > num_orbitals || jorb > num_orbitals){return;}
    
    int idx1 = idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb;
    double phase =   ax*r_x[idx_r] 
                   + ay*r_y[idx_r]
                   + az*r_z[idx_r];
    cdouble_cuda peierls_phase = make_cuDoubleComplex(cos(phase), sin(phase));
    cdouble_cuda op_value = op[idx1];
    cdouble_cuda rho_value = cuConj(cuCmul(cuConj(peierls_phase),rho[idx1]));
    cdouble_cuda value = cuCmul(op_value, rho_value);     
//    data[ti] = cuCadd(data[ti],value);
    atomicAdd(&reinterpret_cast<double*>(data)[2*ti], value.x);    
//    atomicAdd(&reinterpret_cast<double*>(data)[2*ti+1], value.y);
}

__global__ void kernel_calc_observable_shared(
    cdouble_cuda* op,
    cdouble_cuda* rho,
    cdouble_cuda* data,
    double* r_x, double* r_y, double* r_z,
    double ax, double ay, double az,
    const int ti,
    const int num_points,
    const int num_orbitals)
{
    extern __shared__ double smem[];
    double* s_re = smem;
    double* s_im = smem + blockDim.x;

    int idx_r = blockIdx.x * blockDim.x + threadIdx.x;
    int iorb  = blockIdx.y;
    int jorb  = blockIdx.z;
    int tid   = threadIdx.x;

    // --- each thread computes its contribution ---
    double re = 0.0, im = 0.0;
    if (idx_r < num_points && iorb < num_orbitals && jorb < num_orbitals) {
        int idx1 = idx_r * num_orbitals * num_orbitals
                 + iorb  * num_orbitals
                 + jorb;

        double phase = ax * r_x[idx_r]
                     + ay * r_y[idx_r]
                     + az * r_z[idx_r];

        cdouble_cuda peierls_phase = make_cuDoubleComplex(cos(phase), sin(phase));
        cdouble_cuda op_value      = op[idx1];
        cdouble_cuda rho_value     = cuConj(cuCmul(cuConj(peierls_phase), rho[idx1]));
        cdouble_cuda value         = cuCmul(op_value, rho_value);

        re = value.x;
        im = value.y;
    }
    s_re[tid] = re;
    s_im[tid] = im;
    __syncthreads();

    // --- tree reduction in shared memory ---
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            s_re[tid] += s_re[tid + stride];
            s_im[tid] += s_im[tid + stride];
        }
        __syncthreads();
    }

    // --- one atomicAdd per block (not per thread) ---
    if (tid == 0) {
        double* base = reinterpret_cast<double*>(data);
        atomicAdd(&base[ti * 2],     s_re[0]);
        atomicAdd(&base[ti * 2 + 1], s_im[0]);
    }
}

void call_kernel_calc_observable(cdouble_cuda* op, cdouble_cuda* rho, cdouble_cuda* data, double* r_x, double* r_y, double* r_z, double ax, double ay, double az, const int ti, const int num_points, const int num_orbitals, cudaStream_t* stream){

    dim3 block(THREADS_CUDA);
    dim3 grid((num_points + THREADS_CUDA -1)/THREADS_CUDA, num_orbitals, num_orbitals);

    int smem = 2 * THREADS_CUDA * sizeof(double);

    kernel_calc_observable_shared<<<grid, block, smem, *stream>>>(op, rho, data, r_x, r_y, r_z, ax, ay, az, ti, num_points, num_orbitals);
}
#endif
