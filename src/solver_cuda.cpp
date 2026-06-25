#ifdef HAVE_CUDA
#include <iostream>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include "cblas.h"
#include "solver.h"
#include "vec3_util.h"
#include "fftw_helper.h"
#include "cuda_helpers/helper_cuda.h"
#include "kernels_solver_cuda.h"
#include "solver_cuda.h"

#define MOD(a, b) (((a) % (b) + (b)) % (b))

//double _interpolate(const int n, const double theta, double* arr, const int nmax){
//    if(n>nmax-1){return arr[nmax-1];}
//    return (1.0-theta)*arr[n] + theta*arr[n-1];
//}

Solver_cuda::Solver_cuda(){

}

void Solver_cuda::init(){
    _init_device_arrays();
    _create_cufft_plan(); 
    return;
}

Solver_cuda::Solver_cuda(Settings_swe *settings_swe, Grid* grid, Hamiltonian* hamiltonian, RDM_cuda* rho, BerryConnection** r_bc, Efield* efield, WannierTB* wannier, cudaStream_t* stream){
    _settings_swe = settings_swe; 
    _grid = grid;
    _hamiltonian = hamiltonian;
    _rho = rho;
    _efield = efield;
    _wannier = wannier;
    _r_bc = r_bc;

    _num_points = _settings_swe->nr1 * _settings_swe->nr2 * _settings_swe->nr3;
    _num_orbitals = _settings_swe->num_orb;

    _stream = stream;
    _d_rho = _rho->get_ptr_cuda();
    
    _allocate();
}

void Solver_cuda::_allocate(){
    checkCudaErrors(cudaMallocAsync((void**)&_d_k1, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_k2, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_k3, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_k4, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_h0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_xbc, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_ybc, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_zbc, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
//    checkCudaErrors(cudaMallocAsync((void**)&_d_rho, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_heff, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_rho_k, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_heff_k, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_comm_k, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_peierls_phase, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_r_vec_x, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_r_vec_y, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_r_vec_z, sizeof(cdouble_cuda)*_num_points, *_stream));

    checkCudaErrors(cudaMemsetAsync(_d_k1, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k2, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k3, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k4, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_h0, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_xbc, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_ybc, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_zbc, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
//    checkCudaErrors(cudaMemsetAsync(_d_rho, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_heff, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_rho_k, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_heff_k, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_comm_k, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_peierls_phase, 0, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_r_vec_x, 0, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_r_vec_y, 0, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_r_vec_z, 0, sizeof(cdouble_cuda)*_num_points, *_stream));
    checkCudaErrors(cudaStreamSynchronize(*_stream));
}

void Solver_cuda::_deallocate(){
    checkCudaErrors(cudaFreeAsync(_d_k1, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_k2, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_k3, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_k4, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_h0, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_xbc, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_ybc, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_zbc, *_stream));
//    checkCudaErrors(cudaFreeAsync(_d_rho, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_heff,*_stream));
    checkCudaErrors(cudaFreeAsync(_d_rho_k, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_heff_k,*_stream));
    checkCudaErrors(cudaFreeAsync(_d_comm_k,*_stream));
    checkCudaErrors(cudaFreeAsync(_d_peierls_phase,*_stream));
    checkCudaErrors(cudaFreeAsync(_d_r_vec_x,*_stream));
    checkCudaErrors(cudaFreeAsync(_d_r_vec_y,*_stream));
    checkCudaErrors(cudaFreeAsync(_d_r_vec_z,*_stream));
    checkCudaErrors(cudaStreamSynchronize(*_stream));
    cufftDestroy(_cufft_plan);
}

Solver_cuda::~Solver_cuda(){
    _deallocate();

}

void Solver_cuda::_create_cufft_plan(){
    int dims[2] = {_settings_swe->nr1, _settings_swe->nr2};
    int embed[2] = {_settings_swe->nr1, _settings_swe->nr2};
    cufftPlanMany(&_cufft_plan, 2, dims,
                  embed, 1, _settings_swe->nr1*_settings_swe->nr2,
                  embed, 1, _settings_swe->nr1*_settings_swe->nr2,
                  CUFFT_Z2Z, _num_orbitals*_num_orbitals);
    cufftSetStream(_cufft_plan, *_stream);
}

void Solver_cuda::step_rk4(const int ti){
    double t = _grid->t(ti);
    _dt = _grid->t()[1] - _grid->t()[0];
    double ax = 0;
    double ax_dt2 = 0;
    double ax_dt = 0;
    double ex = 0;
    double ex_dt2 = 0;
    double ex_dt = 0;
    double ay = 0;
    double ay_dt2 = 0;
    double ay_dt = 0;
    double ey = 0;
    double ey_dt2 = 0;
    double ey_dt = 0;
    double az = 0;
    double az_dt2 = 0;
    double az_dt = 0;
    double ez = 0;
    double ez_dt2 = 0;
    double ez_dt = 0;
    
    if(ti<_settings_swe->nt-1){
        ax = _efield->A_x[ti];
        ax_dt2 = 0.5*_efield->A_x[ti] + 0.5*_efield->A_x[ti+1];
        ax_dt = _efield->A_x[ti+1];
        ex = _efield->E_x[ti];
        ex_dt2 = 0.5*_efield->E_x[ti] + 0.5*_efield->E_x[ti+1];
        ex_dt = _efield->E_x[ti+1];
        ay = _efield->A_y[ti];
        ay_dt2 = 0.5*_efield->A_y[ti] + 0.5*_efield->A_y[ti+1];
        ay_dt = _efield->A_y[ti+1];
        ey = _efield->E_y[ti];
        ey_dt2 = 0.5*_efield->E_y[ti] + 0.5*_efield->E_y[ti+1];
        ey_dt = _efield->E_y[ti+1];
        az = _efield->A_z[ti];
        az_dt2 = 0.5*_efield->A_z[ti] + 0.5*_efield->A_z[ti+1];
        az_dt = _efield->A_z[ti+1];
        ez = _efield->E_z[ti];
        ez_dt2 = 0.5*_efield->E_z[ti] + 0.5*_efield->E_z[ti+1];
        ez_dt = _efield->E_z[ti+1];
    }
    else {
        ax = _efield->A_x[_settings_swe->nt-1];
        ax_dt2 = _efield->A_x[_settings_swe->nt-1];
        ax_dt = _efield->A_x[_settings_swe->nt-1];
        ex = _efield->E_x[_settings_swe->nt-1];
        ex_dt2 = _efield->E_x[_settings_swe->nt-1];
        ex_dt = _efield->E_x[_settings_swe->nt-1];
        ay = _efield->A_y[_settings_swe->nt-1];
        ay_dt2 = _efield->A_y[_settings_swe->nt-1];
        ay_dt = _efield->A_y[_settings_swe->nt-1];
        ey = _efield->E_y[_settings_swe->nt-1];
        ey_dt2 = _efield->E_y[_settings_swe->nt-1];
        ey_dt = _efield->E_y[_settings_swe->nt-1];
        az = _efield->A_z[_settings_swe->nt-1];
        az_dt2 = _efield->A_z[_settings_swe->nt-1];
        az_dt = _efield->A_z[_settings_swe->nt-1];
        ez = _efield->E_z[_settings_swe->nt-1];
        ez_dt2 = _efield->E_z[_settings_swe->nt-1];
        ez_dt = _efield->E_z[_settings_swe->nt-1];
    }
    _clear_kn();

    _update_heff(ex, ey, ez, ax, ay, az);
    _calculate_peierls_phase(ax, ay, az); 
    _update_k1_conv_cuda(ex,  ey, ez, ax, ay, az); 

    _update_heff(ex_dt2, ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2);
    _calculate_peierls_phase(ax_dt2, ay_dt2, az_dt2); 
    _update_k2_conv_cuda(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    _update_k3_conv_cuda(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 

    _update_heff(ex_dt, ey_dt, ez_dt, ax_dt, ay_dt, az_dt);
    _calculate_peierls_phase(ax_dt, ay_dt, az_dt); 
    _update_k4_conv_cuda(ex_dt,  ey_dt, ez_dt, ax_dt, ay_dt, az_dt); 
     
    _step_rho();
    //_copy_rho_device_to_host();
}

void Solver_cuda::_init_device_arrays(){
    cdouble *tmp = new cdouble[_num_points*_num_orbitals*_num_orbitals];
    double* tmp_r_vec_x = new double[_num_points];
    double* tmp_r_vec_y = new double[_num_points];
    double* tmp_r_vec_z = new double[_num_points];

//    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
//        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _rho->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
//    }
//    checkCudaErrors(cudaMemcpyAsync(_d_rho, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _hamiltonian->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
    }
    checkCudaErrors(cudaMemcpyAsync(_d_h0, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _r_bc[0]->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
    }
    checkCudaErrors(cudaMemcpyAsync(_d_xbc, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _r_bc[1]->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
    }
    checkCudaErrors(cudaMemcpyAsync(_d_ybc, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _r_bc[2]->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
    }
    checkCudaErrors(cudaMemcpyAsync(_d_zbc, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

    for(int idx_r = 0; idx_r<_num_points; idx_r++){
        tmp_r_vec_x[idx_r] = _grid->Rvecs(idx_r)[0];
        tmp_r_vec_y[idx_r] = _grid->Rvecs(idx_r)[1];
        tmp_r_vec_z[idx_r] = _grid->Rvecs(idx_r)[2];
    }
    checkCudaErrors(cudaMemcpyAsync(_d_r_vec_x, tmp_r_vec_x, sizeof(double)*_num_points, cudaMemcpyHostToDevice,*_stream));
    checkCudaErrors(cudaMemcpyAsync(_d_r_vec_y, tmp_r_vec_y, sizeof(double)*_num_points, cudaMemcpyHostToDevice,*_stream));
    checkCudaErrors(cudaMemcpyAsync(_d_r_vec_z, tmp_r_vec_z, sizeof(double)*_num_points, cudaMemcpyHostToDevice,*_stream));
    checkCudaErrors(cudaStreamSynchronize(*_stream));
    delete[] tmp;
    delete[] tmp_r_vec_x;
    delete[] tmp_r_vec_y;
    delete[] tmp_r_vec_z;
}

void Solver_cuda::_calculate_peierls_phase(const double ax, const double ay, const double az){
    call_kernel_calculate_peierls_phase(ax, ay, az, _d_peierls_phase, _d_r_vec_x, _d_r_vec_y, _d_r_vec_z, _num_points, _stream);
}

void Solver_cuda::_update_heff(const double ex, const double ey, const double ez,
                          const double ax, const double ay, const double az){
    call_kernel_update_heff(ex, ey, ez, _d_heff, _d_h0, _d_xbc, _d_ybc, _d_zbc, _num_points, _num_orbitals, _stream);

}

void Solver_cuda::_update_k1_conv_cuda(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    double dt = _grid->t()[1] - _grid->t()[0];
    call_pipeline_update_k1(_d_k1, _d_heff, _d_rho, _d_heff_k, _d_rho_k, _d_comm_k, _d_peierls_phase, dt, _settings_swe->nr1, _settings_swe->nr2, _num_orbitals, _stream, _cufft_plan); 
}

void Solver_cuda::_update_k2_conv_cuda(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    double dt = _grid->t()[1] - _grid->t()[0];
    call_pipeline_update_k2(_d_k2, _d_k1,  _d_heff, _d_rho, _d_heff_k, _d_rho_k, _d_comm_k, _d_peierls_phase, dt, _settings_swe->nr1, _settings_swe->nr2, _num_orbitals, _stream, _cufft_plan); 
}

void Solver_cuda::_update_k3_conv_cuda(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    double dt = _grid->t()[1] - _grid->t()[0];
    call_pipeline_update_k3(_d_k3, _d_k2,  _d_heff, _d_rho, _d_heff_k, _d_rho_k, _d_comm_k, _d_peierls_phase, dt, _settings_swe->nr1, _settings_swe->nr2, _num_orbitals, _stream, _cufft_plan); 
}

void Solver_cuda::_update_k4_conv_cuda(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    double dt = _grid->t()[1] - _grid->t()[0];
    call_pipeline_update_k4(_d_k4, _d_k3,  _d_heff, _d_rho, _d_heff_k, _d_rho_k, _d_comm_k, _d_peierls_phase, dt, _settings_swe->nr1, _settings_swe->nr2, _num_orbitals, _stream, _cufft_plan); 
}

void Solver_cuda::_step_rho(){
    double dt = _grid->t()[1] - _grid->t()[0];;
    call_kernel_step_rho(_d_rho, _d_k1, _d_k2, _d_k3, _d_k4, dt, _num_points, _num_orbitals, _stream);
}

void Solver_cuda::_copy_rho_device_to_host(){
//    cdouble* tmp = new cdouble[_num_points*_num_orbitals*_num_orbitals];
//    checkCudaErrors(cudaMemcpyAsync(tmp, _d_rho, sizeof(cdouble)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyDeviceToHost, *_stream));
//    checkCudaErrors(cudaStreamSynchronize(*_stream));
//    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
//            memcpy(_rho->data_ptr()[idx_r],tmp + idx_r*_num_orbitals*_num_orbitals, sizeof(cdouble)*_num_orbitals*_num_orbitals);
//    }
//    delete [] tmp;
}

void Solver_cuda::_clear_kn(){
    checkCudaErrors(cudaMemsetAsync(_d_k1, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k2, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k3, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k4, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
}

#endif
