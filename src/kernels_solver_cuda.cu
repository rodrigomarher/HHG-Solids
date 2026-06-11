#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include "cuda_helpers/helper_cuda.h"
#include "kernels_solver_cuda.h"

__global__ void kernel_calculate_peierls_phase(double ax, double ay, double az, cdouble_cuda* d_peierls_phase, double* d_r_vec_x, double* d_r_vec_y, double* d_r_vec_z, const int num_points){
    int idx_r = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx_r>=num_points){return;}

    double phase =   ax*d_r_vec_x[idx_r] 
                   + ay*d_r_vec_y[idx_r]
                   + az*d_r_vec_z[idx_r];
    
    d_peierls_phase[idx_r] = make_cuDoubleComplex(cos(phase), sin(phase));
}

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

__global__ void kernel_reorder_and_peierls_k1(cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, 
                                         cdouble_cuda* d_heff, cdouble_cuda* d_rho,
                                         cdouble_cuda* d_peierls, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x * blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if (idx_r >= num_points || iorb >= num_orbitals || jorb >= num_orbitals){return;}
    int idx_out = iorb*num_orbitals*num_points + jorb*num_points + idx_r;
    d_heff_k[idx_out] = cuCmul(d_peierls[idx_r], d_heff[idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb]);
    d_rho_k[idx_out] = d_rho[idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb];
}

__global__ void kernel_reorder_and_peierls_kn(cdouble_cuda prefac, cdouble_cuda* kn, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k,
                                         cdouble_cuda* d_heff, cdouble_cuda* d_rho,
                                         cdouble_cuda* d_peierls, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x * blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if (idx_r >= num_points || iorb >= num_orbitals || jorb >= num_orbitals){return;}
    int idx_out = iorb*num_orbitals*num_points + jorb*num_points + idx_r;
    d_heff_k[idx_out] = cuCmul(d_peierls[idx_r], d_heff[idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb]);
    d_rho_k[idx_out] = cuCadd(d_rho[idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb], cuCmul(prefac, kn[idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb]));
}

__global__ void kernel_fftshift2D(cdouble_cuda* data, const int nr1, const int nr2, const int num_orbitals){
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    int iorb = blockIdx.z;
    if (row >= nr1/2 || col >= nr2/2 || iorb>= num_orbitals*num_orbitals){return;}
    cdouble_cuda* s = data + iorb*nr1*nr2;
    int q1 = row * nr2 + col;
    int q2 = row * nr2 + (col + nr2/2);
    int q3 = (row + nr1/2) * nr2 + col;
    int q4 = (row + nr1/2) * nr2 + (col + nr2/2);
    cdouble_cuda t;
    t = s[q1]; s[q1] = s[q4]; s[q4] = t;
    t = s[q2]; s[q2] = s[q3]; s[q3] = t;
}

__global__ void kernel_scale_inverse_fft(cdouble_cuda* data, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x * blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if(idx_r >= num_points || iorb >= num_orbitals || jorb>= num_orbitals){return;}
    int idx = iorb*num_orbitals*num_points + jorb*num_points + idx_r;
    data[idx] = cuCmul(make_cuDoubleComplex(1.0/num_points,0.0), data[idx]);
}

__global__ void kernel_commutator_k(cdouble_cuda* d_comm_k, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, const int num_points, const int num_orbitals){
    int idx_k = blockIdx.x*blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if (idx_k >= num_points || iorb>=num_orbitals || jorb>=num_orbitals){return;}
    cdouble_cuda comm_k = make_cuDoubleComplex(0.0,0.0);
    for(int korb = 0; korb<num_orbitals; korb++){
        cdouble_cuda h0 = d_heff_k[iorb*num_points*num_orbitals + korb*num_points + idx_k]; 
        cdouble_cuda h1 = d_heff_k[korb*num_points*num_orbitals + jorb*num_points + idx_k]; 
        cdouble_cuda rho0 = d_rho_k[korb*num_points*num_orbitals + jorb*num_points + idx_k]; 
        cdouble_cuda rho1 = d_rho_k[iorb*num_points*num_orbitals + korb*num_points + idx_k]; 
        comm_k = cuCadd(comm_k, cuCsub(cuCmul(h0,rho0), cuCmul(rho1,h1)));
    }
    d_comm_k[iorb*num_points*num_orbitals + jorb*num_points + idx_k] = comm_k;
}

__global__ void kernel_update_kn(cdouble_cuda* d_kn, cdouble_cuda* d_comm, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x*blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if(idx_r >= num_points || iorb>=num_orbitals || jorb>=num_orbitals){return;}

    cdouble_cuda value = cuCmul(d_comm[iorb*num_points*num_orbitals + jorb*num_points + idx_r], make_cuDoubleComplex(1.0/(double)num_points,0.0));
    d_kn[idx_r*num_orbitals*num_orbitals + iorb*num_orbitals + jorb] = cuCmul(make_cuDoubleComplex(0.0,-1.0), value); 
}

__global__ void kernel_step_rho(cdouble_cuda* d_rho, cdouble_cuda* d_k1, cdouble_cuda* d_k2, cdouble_cuda* d_k3, cdouble_cuda* d_k4, const double dt, const int num_points, const int num_orbitals){
    int idx_r = blockIdx.x*blockDim.x + threadIdx.x;
    int iorb = blockIdx.y;
    int jorb = blockIdx.z;
    if(idx_r >= num_points || iorb >= num_orbitals || jorb >= num_orbitals){return;}

    int idx = idx_r*num_orbitals*num_orbitals + iorb * num_orbitals + jorb;
    cdouble_cuda prefac = make_cuDoubleComplex(0.1666666666666*dt, 0.0);
    cdouble_cuda tmp1 = cuCmul(make_cuDoubleComplex(1.0, 0.0), d_k1[idx]);
    cdouble_cuda tmp2 = cuCmul(make_cuDoubleComplex(2.0, 0.0), d_k2[idx]);
    cdouble_cuda tmp3 = cuCmul(make_cuDoubleComplex(2.0, 0.0), d_k3[idx]);
    cdouble_cuda tmp4 = cuCmul(make_cuDoubleComplex(1.0, 0.0), d_k4[idx]);

    cdouble_cuda value = d_rho[idx]; 
    cdouble_cuda tmp = tmp1;
    tmp = cuCadd(tmp, tmp2);
    tmp = cuCadd(tmp, tmp3);
    tmp = cuCadd(tmp, tmp4);
    tmp = cuCmul(prefac, tmp);
    value = cuCadd(value, tmp);
    d_rho[idx] = value;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void call_kernel_calculate_peierls_phase(double ax, double ay, double az, cdouble_cuda* d_peierls_phase, double* d_r_vec_x, double* d_r_vec_y, double* d_r_vec_z, const int num_points, cudaStream_t* stream){
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA - 1)/THREADS_CUDA);
        kernel_calculate_peierls_phase<<<grid, block, 0, *stream>>>(ax, ay, az, d_peierls_phase, d_r_vec_x, d_r_vec_y, d_r_vec_z, num_points);
        checkCudaErrors(cudaPeekAtLastError());
    }
}

void call_kernel_update_heff(const double ex, const double ey, const double ez, cdouble_cuda* d_heff, cdouble_cuda* d_h0, cdouble_cuda* d_xbc, cdouble_cuda* d_ybc, cdouble_cuda* d_zbc, const int num_points, const int num_orbitals, cudaStream_t* stream){
    int threadsPerBlock = THREADS_CUDA;
    int blocks = (num_points + threadsPerBlock - 1)/threadsPerBlock; 
    kernel_update_heff<<<blocks, threadsPerBlock, 0, *stream>>>(ex, ey, ez, d_heff, d_h0, d_xbc, d_ybc, d_zbc, num_points, num_orbitals);
    checkCudaErrors(cudaPeekAtLastError());
}

void call_pipeline_update_k1(cdouble_cuda* d_k1, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan){
    int num_points = nr1*nr2;
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_reorder_and_peierls_k1<<<grid, block,0, *stream>>>(d_heff_k, d_rho_k, d_heff, d_rho, d_peierls_phase, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_heff_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_rho_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        cufftExecZ2Z(fft_plan, d_heff_k, d_heff_k, CUFFT_FORWARD);   
        cufftExecZ2Z(fft_plan, d_rho_k, d_rho_k, CUFFT_FORWARD);   
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_commutator_k<<<grid, block, 0, *stream>>>(d_comm_k, d_heff_k, d_rho_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
    
    {
        cufftExecZ2Z(fft_plan, d_comm_k, d_comm_k, CUFFT_INVERSE);   
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_comm_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_update_kn<<<grid, block, 0, *stream>>>(d_k1, d_comm_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
}

void call_pipeline_update_k2(cdouble_cuda* d_k2, cdouble_cuda* d_k1, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan){
    int num_points = nr1*nr2;
    cdouble_cuda prefac = make_cuDoubleComplex(0.5*dt, 0.0);
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_reorder_and_peierls_kn<<<grid, block,0, *stream>>>(prefac, d_k1, d_heff_k, d_rho_k, d_heff, d_rho, d_peierls_phase, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_heff_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_rho_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        cufftExecZ2Z(fft_plan, d_heff_k, d_heff_k, CUFFT_FORWARD);   
        cufftExecZ2Z(fft_plan, d_rho_k, d_rho_k, CUFFT_FORWARD);   
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_commutator_k<<<grid, block, 0, *stream>>>(d_comm_k, d_heff_k, d_rho_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
    
    {
        cufftExecZ2Z(fft_plan, d_comm_k, d_comm_k, CUFFT_INVERSE);   
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_comm_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_update_kn<<<grid, block, 0, *stream>>>(d_k2, d_comm_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
}

void call_pipeline_update_k3(cdouble_cuda* d_k3, cdouble_cuda* d_k2, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan){
    int num_points = nr1*nr2;
    cdouble_cuda prefac = make_cuDoubleComplex(0.5*dt, 0.0);
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_reorder_and_peierls_kn<<<grid, block,0, *stream>>>(prefac, d_k2, d_heff_k, d_rho_k, d_heff, d_rho, d_peierls_phase, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_heff_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_rho_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        cufftExecZ2Z(fft_plan, d_heff_k, d_heff_k, CUFFT_FORWARD);   
        cufftExecZ2Z(fft_plan, d_rho_k, d_rho_k, CUFFT_FORWARD);   
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_commutator_k<<<grid, block, 0, *stream>>>(d_comm_k, d_heff_k, d_rho_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
    
    {
        cufftExecZ2Z(fft_plan, d_comm_k, d_comm_k, CUFFT_INVERSE);   
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_comm_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_update_kn<<<grid, block, 0, *stream>>>(d_k3, d_comm_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
}

void call_pipeline_update_k4(cdouble_cuda* d_k4, cdouble_cuda* d_k3, cdouble_cuda* d_heff, cdouble_cuda* d_rho, cdouble_cuda* d_heff_k, cdouble_cuda* d_rho_k, cdouble_cuda* d_comm_k, cdouble_cuda* d_peierls_phase, const double dt, const int nr1, const int nr2, const int num_orbitals, cudaStream_t* stream, cufftHandle fft_plan){
    int num_points = nr1*nr2;
    cdouble_cuda prefac = make_cuDoubleComplex(dt, 0.0);
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_reorder_and_peierls_kn<<<grid, block,0, *stream>>>(prefac, d_k3, d_heff_k, d_rho_k, d_heff, d_rho, d_peierls_phase, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_heff_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_rho_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
        cufftExecZ2Z(fft_plan, d_heff_k, d_heff_k, CUFFT_FORWARD);   
        cufftExecZ2Z(fft_plan, d_rho_k, d_rho_k, CUFFT_FORWARD);   
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_commutator_k<<<grid, block, 0, *stream>>>(d_comm_k, d_heff_k, d_rho_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
    
    {
        cufftExecZ2Z(fft_plan, d_comm_k, d_comm_k, CUFFT_INVERSE);   
        dim3 block(16,16);
        dim3 grid((nr2 + 16 - 1)/16, (nr1 + 16 -1)/16, num_orbitals*num_orbitals);
        kernel_fftshift2D<<<grid, block, 0, *stream>>>(d_comm_k, nr1, nr2, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }

    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_update_kn<<<grid, block, 0, *stream>>>(d_k4, d_comm_k, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }
}


void call_kernel_step_rho(cdouble_cuda* d_rho, cdouble_cuda* d_k1, cdouble_cuda* d_k2, cdouble_cuda* d_k3, cdouble_cuda* d_k4, const double dt, const int num_points, const int num_orbitals, cudaStream_t* stream){
    {
        dim3 block(THREADS_CUDA);
        dim3 grid((num_points + THREADS_CUDA-1)/THREADS_CUDA, num_orbitals, num_orbitals);
        kernel_step_rho<<<grid, block, 0, *stream>>>(d_rho, d_k1, d_k2, d_k3, d_k4, dt, num_points, num_orbitals);
        checkCudaErrors(cudaPeekAtLastError());
    }   
}
