#include <fstream>
#include <iostream>
#include "observable_cuda.h"
#include "kernels_observable_cuda.h"
#include "cuda_helpers/helper_cuda.h"

Observable_cuda::Observable_cuda(Settings* settings, Grid * grid, RDM_cuda* rho, Operator* op, cudaStream_t* stream){
    _unit_system = grid->unit_system();
    _settings = settings;
    _rho_cuda = rho;
    _operator = op;
    _stream = stream;
    _grid = grid;
    
    _num_points = settings->nr1*settings->nr2*settings->nr3;
    _num_orbitals = settings->num_orb;
    _data = new cdouble[_settings->nt];

    _allocate();
//    _init_device();
}

void Observable_cuda::_allocate(){
    checkCudaErrors(cudaMallocAsync((void**)&_d_data, sizeof(cdouble_cuda)*_settings->nt, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_data, 0, sizeof(cdouble_cuda)*_settings->nt, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_operator, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_operator, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_r_vec_x, sizeof(double)*_num_points, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_r_vec_x, 0, sizeof(double)*_num_points, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_r_vec_y, sizeof(double)*_num_points, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_r_vec_y, 0, sizeof(double)*_num_points, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_r_vec_z, sizeof(double)*_num_points, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_r_vec_z, 0, sizeof(double)*_num_points, *_stream));
}

void Observable_cuda::init_device(){
    cdouble* tmp = new cdouble[_num_points*_num_orbitals*_num_orbitals];
    double* tmp_r_vec_x = new double[_num_points]; 
    double* tmp_r_vec_y = new double[_num_points]; 
    double* tmp_r_vec_z = new double[_num_points]; 
    for(int idx_r = 0; idx_r<_num_points; idx_r++){
        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _operator->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
        tmp_r_vec_x[idx_r] = _grid->Rvecs(idx_r)[0];
        tmp_r_vec_y[idx_r] = _grid->Rvecs(idx_r)[1];
        tmp_r_vec_z[idx_r] = _grid->Rvecs(idx_r)[2];
    }
    checkCudaErrors(cudaMemcpyAsync(_d_operator, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

    checkCudaErrors(cudaMemcpyAsync(_d_r_vec_x, tmp_r_vec_x, sizeof(double)*_num_points, cudaMemcpyHostToDevice, *_stream));
    checkCudaErrors(cudaMemcpyAsync(_d_r_vec_y, tmp_r_vec_y, sizeof(double)*_num_points, cudaMemcpyHostToDevice, *_stream));
    checkCudaErrors(cudaMemcpyAsync(_d_r_vec_z, tmp_r_vec_z, sizeof(double)*_num_points, cudaMemcpyHostToDevice, *_stream));
    
    checkCudaErrors(cudaStreamSynchronize(*_stream));
    delete[] tmp;
    delete[] tmp_r_vec_x;
    delete[] tmp_r_vec_y;
    delete[] tmp_r_vec_z;
}

void Observable_cuda::_deallocate(){
    checkCudaErrors(cudaFreeAsync(_d_data, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_r_vec_x, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_r_vec_y, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_r_vec_z, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_operator, *_stream));
}

void Observable_cuda::calculate(int ti, double ax, double ay, double az){
    cdouble_cuda* d_operator = _d_operator;
    cdouble_cuda* d_rho = _rho_cuda->get_ptr_cuda();
    cdouble_cuda* d_data = _d_data;
    double* d_r_x = _d_r_vec_x;
    double* d_r_y = _d_r_vec_y;
    double* d_r_z = _d_r_vec_z;
    
    call_kernel_calc_observable(d_operator, d_rho, d_data, d_r_x, d_r_y, d_r_z, ax, ay, az, ti, _num_points, _num_orbitals, _stream);
}


cdouble* Observable_cuda::get_ptr(){
    return _data;
}

cdouble_cuda* Observable_cuda::get_d_ptr(){
    return _d_data;
}

void Observable_cuda::copy_device_to_host(){
    checkCudaErrors(cudaMemcpyAsync(_data, _d_data, sizeof(cdouble_cuda)*_settings->nt, cudaMemcpyDeviceToHost, *_stream));
}

void Observable_cuda::write(std::string filename){
    std::ofstream file(filename);
    if(file.is_open()){
        file << _unit_system << std::endl;
        for(int i=0; i <_settings->nt; i++){
            file << std::abs(_data[i]) << " " <<std::atan2(_data[i].imag(), _data[i].real()) << std::endl;
        }
        file.close();
    }
    else{
        std::cout<< "[Observable_cuda] Error opening file"<< filename<<std::endl;
    }
}

Observable_cuda::~Observable_cuda(){
    _deallocate();
    delete[] _data;
}

