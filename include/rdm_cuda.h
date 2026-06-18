#pragma once

#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include "rdm.h"
#include "settings.h"

#define cdouble_cuda cuDoubleComplex
class RDM_cuda{
    private:
        RDM* _rho_cpu;
        Settings* _settings;
        cdouble_cuda *_d_data;
        cudaStream_t* _stream;
        cufftHandle _cufft_plan;
        int _space_type;
        int _gauge;
        int _unit_system;
        void _allocate();
        void _deallocate();
        void _create_cufft_plan();

    public:
        RDM_cuda(Settings* settings, RDM* rho, cudaStream_t* stream);
        void copy_cpu_to_gpu();
        void copy_gpu_to_cpu();
        void convert_to_k();
        void convert_to_r();
        cdouble_cuda* get_ptr_cuda();
        ~RDM_cuda();
};

#endif
