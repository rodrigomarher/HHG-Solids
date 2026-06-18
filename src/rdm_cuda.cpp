#ifdef HAVE_CUDA

#include <iostream>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cuda_helpers/helper_cuda.h"
#include "settings.h"
#include "rdm.h"
#include "rdm_cuda.h"
#include "kernels_rdm_cuda.h"

#define cdouble std::complex<double>
#define cdouble_cuda cuDoubleComplex

RDM_cuda::RDM_cuda(Settings* settings, RDM* rho, cudaStream_t* stream){
    _rho_cpu = rho;
    _settings = settings;
    _stream = stream;
    _allocate();
    _create_cufft_plan();
    copy_cpu_to_gpu();
}

void RDM_cuda::_allocate(){
    int nr1 = _settings->nr1;
    int nr2 = _settings->nr2;
    int nr3 = _settings->nr3;
    int num_orbitals = _settings->num_orb;
    int num_points = nr1*nr2*nr3;
    int num_total = num_points*num_orbitals*num_orbitals;
    checkCudaErrors(cudaMallocAsync((void**)&_d_data, sizeof(cdouble_cuda)*num_total, *_stream)); 
    checkCudaErrors(cudaMemsetAsync(_d_data, 0, sizeof(cdouble_cuda)*num_total, *_stream));
    
    checkCudaErrors(cudaStreamSynchronize(*_stream));
}

void RDM_cuda::_deallocate(){
    checkCudaErrors(cudaFreeAsync(_d_data, *_stream));
    cufftDestroy(_cufft_plan);
}

void RDM_cuda::_create_cufft_plan(){
    int dims[2] = {_settings->nr1, _settings->nr2};
    int embed[2] = {_settings->nr1, _settings->nr2};
    cufftPlanMany(&_cufft_plan, 2, dims, 
                  embed, _settings->num_orb*_settings->num_orb, 1,
                  embed, _settings->num_orb*_settings->num_orb, 1,
                  CUFFT_Z2Z, _settings->num_orb*_settings->num_orb);
    cufftSetStream(_cufft_plan, *_stream);
}

void RDM_cuda::copy_cpu_to_gpu(){
    int num_points = _settings->nr1*_settings->nr2*_settings->nr3;
    int num_orbitals = _settings->num_orb;
    cdouble *tmp = new cdouble[num_points*num_orbitals*num_orbitals];

    for(int idx = 0; idx<num_points; idx++){
        memcpy(tmp + idx*num_orbitals*num_orbitals, _rho_cpu->data_ptr()[idx], sizeof(cdouble)*num_orbitals*num_orbitals);
    }
    
    checkCudaErrors(cudaMemcpyAsync(_d_data, tmp, sizeof(cdouble_cuda)*num_orbitals*num_orbitals*num_points, cudaMemcpyHostToDevice,*_stream));
    
    delete [] tmp;
    
    _gauge = _rho_cpu->gauge(); 
    _space_type = _rho_cpu->space_type();
    _unit_system = _rho_cpu->unit_system();
}


void RDM_cuda::copy_gpu_to_cpu(){
    if(_gauge != _rho_cpu->gauge() ||
       _space_type != _rho_cpu->space_type() ||
       _unit_system != _rho_cpu->unit_system()){
        std::cout<<"[RDM_cuda::copy_gpu_to_cpu] Device gauge, space_type or unit_system not consistent with cpu RDM."<<std::endl;
        std::exit(1);
    }
    int num_points = _settings->nr1*_settings->nr2*_settings->nr3;
    int num_orbitals = _settings->num_orb;
    cdouble *tmp = new cdouble[num_orbitals*num_orbitals*num_points];

    checkCudaErrors(cudaMemcpyAsync(tmp, _d_data, sizeof(cdouble_cuda)*num_points*num_orbitals*num_orbitals, cudaMemcpyDeviceToHost, *_stream ));

    for(int idx = 0; idx<num_points; idx++){
        memcpy(_rho_cpu->data_ptr()[idx], tmp + idx*num_orbitals*num_orbitals, sizeof(cdouble)*num_orbitals*num_orbitals);
    }
    delete[] tmp;
}

void RDM_cuda::convert_to_k(){
    if(_space_type == KSPACE){return;}
    int num_points = _settings->nr1*_settings->nr2;
    int num_orbitals = _settings->num_orb;
    cufftExecZ2Z(_cufft_plan, _d_data, _d_data, CUFFT_FORWARD);   
    _space_type = KSPACE;
}

void RDM_cuda::convert_to_r(){
    if(_space_type == RSPACE){return;}
    int num_points = _settings->nr1*_settings->nr2;
    int num_orbitals = _settings->num_orb;
    cufftExecZ2Z(_cufft_plan, _d_data, _d_data, CUFFT_INVERSE);   
    call_kernel_normalize_n(_d_data, num_points, num_orbitals, _stream);
    _space_type = RSPACE;
}

cdouble_cuda* RDM_cuda::get_ptr_cuda(){
    return _d_data;
}

RDM_cuda::~RDM_cuda(){
    _deallocate();
}
#endif
