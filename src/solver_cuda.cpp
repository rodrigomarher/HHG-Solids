#ifdef HAVE_CUDA
#include <iostream>
#include <cuda_runtime.h>
#include <cuComplex.h>
#include "cblas.h"
#include "solver.h"
#include "vec3_util.h"
#include "fftw_helper.h"
#include "cuda_helpers/helper_cuda.h"
#include "cuda_kernels.h"
#include "solver_cuda.h"

#define MOD(a, b) (((a) % (b) + (b)) % (b))

//double _interpolate(const int n, const double theta, double* arr, const int nmax){
//    if(n>nmax-1){return arr[nmax-1];}
//    return (1.0-theta)*arr[n] + theta*arr[n-1];
//}

Solver_cuda::Solver_cuda(){

}

void Solver_cuda::init(){
    return;
}

Solver_cuda::Solver_cuda(Settings *settings, Grid* grid, Hamiltonian* hamiltonian, RDM* rho, BerryConnection** r_bc, Efield* efield, WannierTB* wannier){
    _settings = settings; 
    _grid = grid;
    _hamiltonian = hamiltonian;
    _rho = rho;
    _efield = efield;
    _wannier = wannier;
    _r_bc = r_bc;

    _num_points = _settings->nr1 * _settings->nr2 * _settings->nr3;
    _num_orbitals = _settings->num_orb;

    _allocate();
    _init_device_h0_rbc_rho();
    
}

void Solver_cuda::_allocate(){
    checkCudaErrors(cudaMallocAsync((void**)&_d_k1, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_k2, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_k3, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_k4, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_heff, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_h0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_xbc, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_ybc, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_zbc, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMallocAsync((void**)&_d_rho, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));

    checkCudaErrors(cudaMemsetAsync(_d_k1, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k2, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k3, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k4, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_heff, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_h0, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_xbc, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_ybc, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaMemsetAsync(_d_zbc, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, *_stream));
    checkCudaErrors(cudaStreamSynchronize(*_stream));
}

void Solver_cuda::_deallocate(){
    checkCudaErrors(cudaFreeAsync(_d_k1, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_k2, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_k3, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_k4, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_heff, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_h0, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_xbc, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_ybc, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_zbc, *_stream));
    checkCudaErrors(cudaFreeAsync(_d_rho, *_stream));
    checkCudaErrors(cudaStreamSynchronize(*_stream));
}

Solver_cuda::~Solver_cuda(){
    _deallocate();
}



void Solver_cuda::step_rk4(const int ti, cdouble* peierls_phase){
    _peierls_phase = peierls_phase;
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
    
    if(ti<_settings->nt-1){
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
        ax = _efield->A_x[_settings->nt-1];
        ax_dt2 = _efield->A_x[_settings->nt-1];
        ax_dt = _efield->A_x[_settings->nt-1];
        ex = _efield->E_x[_settings->nt-1];
        ex_dt2 = _efield->E_x[_settings->nt-1];
        ex_dt = _efield->E_x[_settings->nt-1];
        ay = _efield->A_y[_settings->nt-1];
        ay_dt2 = _efield->A_y[_settings->nt-1];
        ay_dt = _efield->A_y[_settings->nt-1];
        ey = _efield->E_y[_settings->nt-1];
        ey_dt2 = _efield->E_y[_settings->nt-1];
        ey_dt = _efield->E_y[_settings->nt-1];
        az = _efield->A_z[_settings->nt-1];
        az_dt2 = _efield->A_z[_settings->nt-1];
        az_dt = _efield->A_z[_settings->nt-1];
        ez = _efield->E_z[_settings->nt-1];
        ez_dt2 = _efield->E_z[_settings->nt-1];
        ez_dt = _efield->E_z[_settings->nt-1];
    }
    _clear_kn();

    _update_heff(ex, ey, ez, ax, ay, az);
    //_update_k1_conv_fftw(ex,  ey, ez, ax, ay, az); 
    //_update_heff(ex_dt2, ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2);
    //_update_k2_conv_fftw(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    //_update_k3_conv_fftw(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    //_update_heff(ex_dt, ey_dt, ez_dt, ax_dt, ay_dt, az_dt);
    //_update_k4_conv_fftw(ex_dt,  ey_dt, ez_dt, ax_dt, ay_dt, az_dt); 

    //cdouble prefac = {0.1666666666666*_dt, 0.0};
    //for(int idx_r = 0; idx_r < _num_points; idx_r++){
    //    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
    //        cdouble value = _rho->get(idx_r, iorb);
    //        value = value + prefac*(_k1[idx_r][iorb] 
    //                                + 2.0*_k2[idx_r][iorb]
    //                                + 2.0*_k3[idx_r][iorb]
    //                                + _k4[idx_r][iorb]);
    //        _rho->set(value, idx_r, iorb);
    //    }
    //}
}

 Solver_cuda::_init_device_h0_rbc_rho(){
    cdouble *tmp = new cdouble[_num_points*_num_orbitals*_num_orbitals];

    for(int idx_r = 0; idx_r<_num_points; idx_r ++){
        memcpy(tmp + idx_r*_num_orbitals*_num_orbitals, _rho->data_ptr()[idx_r], sizeof(cdouble)*_num_orbitals*_num_orbitals);
    }
    checkCudaErrors(cudaMemcpyAsync(_d_rho, tmp, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals, cudaMemcpyHostToDevice, *_stream));

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
    checkCudaErrors(cudaStreamSynchronize(*_stream));
    delete[] tmp;
}

void Solver_cuda::_calc_commutator(cdouble *A, cdouble *B, cdouble *C, const int n){
    for ( int iorb = 0; iorb < n; iorb++){
                for (int jorb = 0; jorb < n; jorb++){
                    cdouble value_comm = 0;
                    for (int korb = 0; korb<n; korb++){
                        cdouble a_0 =  A[iorb*_num_orbitals + korb];
                        cdouble a_1 =  A[korb*_num_orbitals + jorb];
                        cdouble b_0 = B[korb*_num_orbitals + jorb];
                        cdouble b_1 = B[iorb*_num_orbitals + korb];
                        value_comm += a_0 * b_0 - a_1*b_1;
                    }
                    C[iorb*_num_orbitals + jorb] += value_comm; 
                }
            }
}

void Solver_cuda::_update_heff(const double ex, const double ey, const double ez,
                          const double ax, const double ay, const double az){
    call_kernel_update_heff(ex, ey, ez, _d_heff, _d_h0, _d_xbc, _d_ybc, _d_zbc, _num_points, _num_orbitals, _stream);

    //cdouble **h0 = _hamiltonian->data_ptr();
    //cdouble **xbc = _r_bc[0]->data_ptr();
    //cdouble **ybc = _r_bc[1]->data_ptr();
    //cdouble **zbc = _r_bc[2]->data_ptr();

 
   //cdouble heff_value = {0.0,0.0};
    //for(int idx_r=0; idx_r<_num_points; idx_r++){
    //    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
    //        heff_value = h0[idx_r][iorb] + ex*xbc[idx_r][iorb] + ey* ybc[idx_r][iorb] + ez*zbc[idx_r][iorb];
    //        _heff[idx_r][iorb] = heff_value;
    //    }
    //}
}


void Solver_cuda::_update_k1_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
   // cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *tmp_1 = new cdouble[_num_points];
   // cdouble *tmp_2 = new cdouble[_num_points];
   // cdouble **rho_ptr = _rho->data_ptr();

   // fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
   // fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for(int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
   //             tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
   //             tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb];
   //         }
   //         fftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         fftshift(tmp_2, _settings->nr1, _settings->nr2);
   //         fft3(tmp_1, in, out, _num_points, forward);
   //         fft3(tmp_2, in, out, _num_points, forward);

   //         for(int idx_k = 0; idx_k < _num_points; idx_k++){
   //             heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
   //             rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
   //         }
   //     }
   // }
 
   // for(int idx_k=0; idx_k <_num_points; idx_k++){
   //     for(int iorb = 0; iorb<_num_orbitals; iorb++){
   //         for(int jorb = 0; jorb <_num_orbitals; jorb++){
   //             cdouble comm_value = {0.0, 0.0};
   //             for (int korb = 0; korb<_num_orbitals; korb++){
   //                 cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 
   //                 cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 comm_value += heff0*rho0 - rho1*heff1;
   //             }
   //             comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
   //         }
   //     }
   // } 

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for (int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_k = 0; idx_k <_num_points; idx_k++){
   //             tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //         }

   //         ifft3(tmp_1, in, out, _num_points, backward);
   //         ifftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             _k1[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
   //         }
   //     }
   // }

   // fftw_free(in);
   // fftw_free(out);
   // fftw_destroy_plan(forward);
   // fftw_destroy_plan(backward);
   // delete[] heff_k;
   // delete[] rho_k;    
   // delete[] comm_k;
   // delete[] tmp_1;
   // delete[] tmp_2;
}

void Solver_cuda::_update_k2_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
   // cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *tmp_1 = new cdouble[_num_points];
   // cdouble *tmp_2 = new cdouble[_num_points];
   // cdouble **rho_ptr = _rho->data_ptr();

   // fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
   // fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for(int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
   //             tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
   //             tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb] + 0.5*_dt*_k1[idx_r][iorb*_num_orbitals + jorb];
   //         }

   //         fftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         fftshift(tmp_2, _settings->nr1, _settings->nr2);
   //         fft3(tmp_1, in, out, _num_points, forward);
   //         fft3(tmp_2, in, out, _num_points, forward);

   //         for(int idx_k = 0; idx_k < _num_points; idx_k++){
   //             heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
   //             rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
   //         }
   //     }
   // }
 
   // for(int idx_k=0; idx_k <_num_points; idx_k++){
   //     for(int iorb = 0; iorb<_num_orbitals; iorb++){
   //         for(int jorb = 0; jorb <_num_orbitals; jorb++){
   //             cdouble comm_value = {0.0, 0.0};
   //             for (int korb = 0; korb<_num_orbitals; korb++){
   //                 cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 
   //                 cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 comm_value += heff0*rho0 - rho1*heff1;
   //             }
   //             comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
   //         }
   //     }
   // } 

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for (int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_k = 0; idx_k <_num_points; idx_k++){
   //             tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //         }

   //         ifft3(tmp_1, in, out, _num_points, backward);
   //         ifftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             _k2[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
   //         }
   //     }
   // }

   // fftw_free(in);
   // fftw_free(out);
   // fftw_destroy_plan(forward);
   // fftw_destroy_plan(backward);
   // delete[] heff_k;
   // delete[] rho_k;    
   // delete[] comm_k;
   // delete[] tmp_1;
   // delete[] tmp_2;
}

void Solver_cuda::_update_k3_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
   // cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *tmp_1 = new cdouble[_num_points];
   // cdouble *tmp_2 = new cdouble[_num_points];
   // cdouble **rho_ptr = _rho->data_ptr();

   // fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
   // fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for(int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
   //             //cdouble peierls_phase = std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
   //             tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
   //             tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb] + 0.5*_dt*_k2[idx_r][iorb*_num_orbitals + jorb];
   //         }

   //         fftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         fftshift(tmp_2, _settings->nr1, _settings->nr2);
   //         fft3(tmp_1, in, out, _num_points, forward);
   //         fft3(tmp_2, in, out, _num_points, forward);

   //         for(int idx_k = 0; idx_k < _num_points; idx_k++){
   //             heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
   //             rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
   //         }
   //     }
   // }
 
   // for(int idx_k=0; idx_k <_num_points; idx_k++){
   //     for(int iorb = 0; iorb<_num_orbitals; iorb++){
   //         for(int jorb = 0; jorb <_num_orbitals; jorb++){
   //             cdouble comm_value = {0.0, 0.0};
   //             for (int korb = 0; korb<_num_orbitals; korb++){
   //                 cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 
   //                 cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 comm_value += heff0*rho0 - rho1*heff1;
   //             }
   //             comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
   //         }
   //     }
   // } 

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for (int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_k = 0; idx_k <_num_points; idx_k++){
   //             tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //         }

   //         ifft3(tmp_1, in, out, _num_points, backward);
   //         ifftshift(tmp_1, _settings->nr1, _settings->nr2);

   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             _k3[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
   //         }
   //     }
   // }

   // fftw_free(in);
   // fftw_free(out);
   // fftw_destroy_plan(forward);
   // fftw_destroy_plan(backward);
   // delete[] heff_k;
   // delete[] rho_k;    
   // delete[] comm_k;
   // delete[] tmp_1;
   // delete[] tmp_2;
}

void Solver_cuda::_update_k4_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
   // cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
   // cdouble *tmp_1 = new cdouble[_num_points];
   // cdouble *tmp_2 = new cdouble[_num_points];
   // cdouble **rho_ptr = _rho->data_ptr();

   // fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
   // fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
   // fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for(int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
   //             //cdouble peierls_phase = std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
   //             tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
   //             tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb] + _dt*_k3[idx_r][iorb*_num_orbitals + jorb];
   //         }
   //         fftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         fftshift(tmp_2, _settings->nr1, _settings->nr2);
   //         fft3(tmp_1, in, out, _num_points, forward);
   //         fft3(tmp_2, in, out, _num_points, forward);

   //         for(int idx_k = 0; idx_k < _num_points; idx_k++){
   //             heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
   //             rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
   //         }
   //     }
   // }
 
   // for(int idx_k=0; idx_k <_num_points; idx_k++){
   //     for(int iorb = 0; iorb<_num_orbitals; iorb++){
   //         for(int jorb = 0; jorb <_num_orbitals; jorb++){
   //             cdouble comm_value = {0.0, 0.0};
   //             for (int korb = 0; korb<_num_orbitals; korb++){
   //                 cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 
   //                 cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //                 cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
   //                 comm_value += heff0*rho0 - rho1*heff1;
   //             }
   //             comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
   //         }
   //     }
   // } 

   // for(int iorb = 0; iorb < _num_orbitals; iorb++){
   //     for (int jorb = 0; jorb < _num_orbitals; jorb++){
   //         for(int idx_k = 0; idx_k <_num_points; idx_k++){
   //             tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
   //         }

   //         ifft3(tmp_1, in, out, _num_points, backward);
   //         ifftshift(tmp_1, _settings->nr1, _settings->nr2);
   //         
   //         for(int idx_r = 0; idx_r < _num_points; idx_r++){
   //             _k4[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
   //         }
   //     }
   // }

   // fftw_free(in);
   // fftw_free(out);
   // fftw_destroy_plan(forward);
   // fftw_destroy_plan(backward);
   // delete[] heff_k;
   // delete[] rho_k;    
   // delete[] comm_k;
   // delete[] tmp_1;
   // delete[] tmp_2;
}

void Solver_cuda::_clear_kn(){
    checkCudaErrors(cudaMemsetAsync(_d_k1, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k2, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k3, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
    checkCudaErrors(cudaMemsetAsync(_d_k4, 0, sizeof(cdouble_cuda)*_num_points*_num_orbitals*_num_orbitals,*_stream));
}

#endif
