#include <iostream>
#include "cblas.h"
#include "solver_kspace.h"
#include "vec3_util.h"
#include "fftw_helper.h"
#include "interpolation.h"

#define MOD(a, b) (((a) % (b) + (b)) % (b))

Solver_kspace::Solver_kspace(){

}

Solver_kspace::Solver_kspace(Settings *settings, Grid* grid, Hamiltonian* hamiltonian, RDM* rho, BerryConnection** r_bc, Efield* efield, WannierTB* wannier){
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
 
}

void Solver_kspace::init(){
    _fill_aux_data();   
    _convert_to_kspace(_h0_k);
    _convert_to_kspace(_xbc_k);
    _convert_to_kspace(_ybc_k);
    _convert_to_kspace(_zbc_k);
    //for(int idx_k = 0; idx_k<_num_points; idx_k++){
    //    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
    //        std::cout<<"idx_k ("<<idx_k<<") iorb ("<<iorb<<") "<<_h0_k[idx_k][iorb]<<std::endl;
    //    }
    //}
}

void Solver_kspace::_allocate(){
    _k1 = new cdouble*[_num_points];
    _k2 = new cdouble*[_num_points];
    _k3 = new cdouble*[_num_points];
    _k4 = new cdouble*[_num_points];
    _heff_k = new cdouble*[_num_points];
    _rho_k = new cdouble*[_num_points];
    _h0_k = new cdouble*[_num_points];
    _xbc_k = new cdouble*[_num_points];
    _ybc_k = new cdouble*[_num_points];
    _zbc_k = new cdouble*[_num_points];

    for(int idx_k = 0; idx_k < _num_points; idx_k++){
        _k1[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _k2[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _k3[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _k4[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _heff_k[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _rho_k[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _h0_k[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _xbc_k[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _ybc_k[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        _zbc_k[idx_k] = new cdouble[_num_orbitals*_num_orbitals];
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            _k1[idx_k][iorb] = cdouble(0.0,0.0);
            _k2[idx_k][iorb] = cdouble(0.0,0.0);
            _k3[idx_k][iorb] = cdouble(0.0,0.0);
            _k4[idx_k][iorb] = cdouble(0.0,0.0);
            _rho_k[idx_k][iorb] = cdouble(0.0,0.0);
            _heff_k[idx_k][iorb] = cdouble(0.0,0.0);
            _h0_k[idx_k][iorb] = cdouble(0.0,0.0);
            _xbc_k[idx_k][iorb] = cdouble(0.0,0.0);
            _ybc_k[idx_k][iorb] = cdouble(0.0,0.0);
            _zbc_k[idx_k][iorb] = cdouble(0.0,0.0);
        }
    }
}

void Solver_kspace::_fill_aux_data(){
    for(int idx_r = 0; idx_r<_num_points; idx_r++){
        for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
            _rho_k[idx_r][iorb] = _rho->data_ptr()[idx_r][iorb];
            _h0_k[idx_r][iorb] = _hamiltonian->data_ptr()[idx_r][iorb];
            _xbc_k[idx_r][iorb] = _r_bc[0]->data_ptr()[idx_r][iorb];
            _ybc_k[idx_r][iorb] = _r_bc[1]->data_ptr()[idx_r][iorb];
            _zbc_k[idx_r][iorb] = _r_bc[2]->data_ptr()[idx_r][iorb];
        }
    }
}

void Solver_kspace::_deallocate(){
    for(int idx_k = 0; idx_k < _num_points; idx_k++){
        delete[] _k1[idx_k];
        delete[] _k2[idx_k];
        delete[] _k3[idx_k];
        delete[] _k4[idx_k];
        delete[] _heff_k[idx_k];
        delete[] _rho_k[idx_k];
        delete[] _h0_k[idx_k];
        delete[] _xbc_k[idx_k];
        delete[] _ybc_k[idx_k];
        delete[] _zbc_k[idx_k];
    }
    delete[] _k1;
    delete[] _k2;
    delete[] _k3;
    delete[] _k4;
    delete[] _heff_k;
    delete[] _rho_k;
    delete[] _h0_k;
    delete[] _xbc_k;
    delete[] _ybc_k;
    delete[] _zbc_k;
}

Solver_kspace::~Solver_kspace(){
    _deallocate();
}

void Solver_kspace::_convert_to_kspace(cdouble** data){
    cdouble* tmp1 = new cdouble[_num_points];
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points); 
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points); 
    fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    
    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
        for(int idx_r = 0; idx_r < _num_points; idx_r++){
            tmp1[idx_r] = data[idx_r][iorb];
        }
        ifftshift(tmp1, _settings->nr1, _settings->nr2);
        fft3(tmp1, in, out, _num_points, forward);
        for(int idx_k = 0; idx_k < _num_points; idx_k++){
            data[idx_k][iorb] = tmp1[idx_k];
        }
    } 
    
    delete[] tmp1;
    fftw_free(in);
    fftw_free(out);
} 

void Solver_kspace::_convert_to_rspace(cdouble** data){
    cdouble* tmp1 = new cdouble[_num_points];
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points); 
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points); 
    fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    
    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
        for(int idx_k = 0; idx_k < _num_points; idx_k++){
            tmp1[idx_k] = data[idx_k][iorb];
        }
        //fftshift(tmp1, _settings->nr1, _settings->nr2);
        ifft3(tmp1, in, out, _num_points, backward);
        ifftshift(tmp1, _settings->nr1, _settings->nr2);
        for(int idx_r = 0; idx_r < _num_points; idx_r++){
            data[idx_r][iorb] = tmp1[idx_r];
        }
    } 
    
    delete[] tmp1;
    fftw_free(in);
    fftw_free(out);
}

void Solver_kspace::step_rk4(const int ti){
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
    //for(int idx_r = 0; idx_r<_num_points; idx_r++){
    //    for (int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
    //        std::cout<<"idx_r ("<<idx_r<<") iorb ("<<iorb<<"): "<<_rho_k[idx_r][iorb]<<std::endl;
    //    }
    //} 
    _clear_kn();
    _convert_to_kspace(_rho_k);
     
    _update_heff(ex, ey, ez, ax, ay, az);
    _update_k1(ex,  ey, ez, ax, ay, az); 
    _update_heff(ex_dt2, ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2);
    _update_k2(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    _update_k3(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    _update_heff(ex_dt, ey_dt, ez_dt, ax_dt, ay_dt, az_dt);
    _update_k4(ex_dt,  ey_dt, ez_dt, ax_dt, ay_dt, az_dt); 
    
    //for(int idx_k = 0; idx_k<_num_points; idx_k++){
    //    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
    //        std::cout<<"idx_k ("<<idx_k<<") iorb ("<<iorb<<") "<<_rho_k[idx_k][iorb]<<std::endl;
    //    }
    //}
    cdouble prefac = {0.1666666666666*_dt, 0.0};
    for(int idx_k = 0; idx_k < _num_points; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
            cdouble value = _rho_k[idx_k][iorb];
            value = value + prefac*(_k1[idx_k][iorb] 
                                    + 2.0*_k2[idx_k][iorb]
                                    + 2.0*_k3[idx_k][iorb]
                                    + _k4[idx_k][iorb]);
            //std::cout<<"idx_k ("<<idx_k<<") iorb ("<<iorb<<"): "<<_k4[idx_k][iorb]<<std::endl;
            _rho_k[idx_k][iorb] = value;
        }
    }
   
    _convert_to_rspace(_rho_k); 
    for(int idx_r = 0; idx_r<_num_points; idx_r++){
        for (int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            _rho->set(_rho_k[idx_r][iorb], idx_r, iorb);
           // std::cout<<"idx_r ("<<idx_r<<") iorb ("<<iorb<<"): "<<_rho_k[idx_r][iorb]<<std::endl;
        }
    } 
}

void Solver_kspace::_update_k1(const double ex, const double ey, const double ez,
                               const double ax, const double ay, const double az){
    int N1 = _settings->nr1;
    int N2 = _settings->nr2;
    int N3 = _settings->nr3;
    cdouble* heff_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* rho_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* comm_k = new cdouble[_num_orbitals*_num_orbitals];

    double **reciprocal_lattice_vector;
    reciprocal_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){reciprocal_lattice_vector[i] = new double[3];}
    _grid->reciprocal_vector(reciprocal_lattice_vector);

    for(int idx_k_1 = 0; idx_k_1 < _settings->nr1; idx_k_1++){
        for(int idx_k_2 = 0; idx_k_2 < _settings->nr2; idx_k_2++){
            double k_vec[3];
            k_vec[0] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[0] - ax;
            k_vec[1] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[1] - ay;
            k_vec[2] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[2] - az;
            double coefs_d[3];
            solve_3by3(reciprocal_lattice_vector, k_vec, coefs_d);
            //coefs_d[0] *= (double)N1;
            //coefs_d[1] *= (double)N2;
            //coefs_d[2] *= (double)N3;
            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                heff_comm_k[iorb] = bicubicInterpolate(_heff_k, N1, N2, coefs_d[0], coefs_d[1], iorb); 
                //heff_comm_k[iorb] =_heff_k[idx_k_1*N2 + idx_k_2][iorb]; 
                rho_comm_k[iorb] = _rho_k[idx_k_1*N2 + idx_k_2][iorb]; 
                //if(std::abs(heff_comm_k[iorb] - _heff_k[idx_k_1*N2 + idx_k_2][iorb])>1e-14){
                //    std::cout<<"idx_k ("<<idx_k_1<<") iorb ("<<iorb<<"): "<<std::abs(heff_comm_k[iorb])
                //             <<" "<< std::abs(_heff_k[idx_k_1*N2 + idx_k_2][iorb])<<std::endl;
                //}
            } 
            _calc_commutator(heff_comm_k, rho_comm_k, comm_k, _num_orbitals);
            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                _k1[idx_k_1*N2 + idx_k_2][iorb] = cdouble(0.0,-1.0)*comm_k[iorb];
            }
        }
    } 
    delete[] heff_comm_k;
    delete[] rho_comm_k;
    delete[] comm_k;
}

void Solver_kspace::_update_k2(const double ex, const double ey, const double ez,
                               const double ax, const double ay, const double az){
    int N1 = _settings->nr1;
    int N2 = _settings->nr2;
    int N3 = _settings->nr3;
    cdouble* heff_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* rho_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* comm_k = new cdouble[_num_orbitals*_num_orbitals];
    double **reciprocal_lattice_vector;
    reciprocal_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){reciprocal_lattice_vector[i] = new double[3];}
    _grid->reciprocal_vector(reciprocal_lattice_vector);
    for(int idx_k_1 = 0; idx_k_1 < _settings->nr1; idx_k_1++){
        for(int idx_k_2 = 0; idx_k_2 < _settings->nr2; idx_k_2++){
            double k_vec[3];
            k_vec[0] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[0] - ax;
            k_vec[1] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[1] - ay;
            k_vec[2] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[2] - az;
            double coefs_d[3];
            solve_3by3(reciprocal_lattice_vector, k_vec, coefs_d);
            //coefs_d[0] *= N1;
            //coefs_d[1] *= N2;
            //coefs_d[2] *= N3;

            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                //heff_comm_k[iorb] =_heff_k[idx_k_1*N2 + idx_k_2][iorb]; 
                heff_comm_k[iorb] = bicubicInterpolate(_heff_k, N1, N2, coefs_d[0], coefs_d[1], iorb); 
                rho_comm_k[iorb] = _rho_k[idx_k_1*N2 + idx_k_2][iorb] + 0.5*_dt*_k1[idx_k_1*N2 + idx_k_2][iorb]; 
            } 
            _calc_commutator(heff_comm_k, rho_comm_k, comm_k, _num_orbitals);
            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                _k2[idx_k_1*N2 + idx_k_2][iorb] = cdouble(0.0,-1.0)*comm_k[iorb];
            }
        }
    } 
    delete[] heff_comm_k;
    delete[] rho_comm_k;
    delete[] comm_k;
}

void Solver_kspace::_update_k3(const double ex, const double ey, const double ez,
                               const double ax, const double ay, const double az){
    int N1 = _settings->nr1;
    int N2 = _settings->nr2;
    int N3 = _settings->nr3;
    cdouble* heff_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* rho_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* comm_k = new cdouble[_num_orbitals*_num_orbitals];
    double **reciprocal_lattice_vector;
    reciprocal_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){reciprocal_lattice_vector[i] = new double[3];}
    _grid->reciprocal_vector(reciprocal_lattice_vector);
    for(int idx_k_1 = 0; idx_k_1 < _settings->nr1; idx_k_1++){
        for(int idx_k_2 = 0; idx_k_2 < _settings->nr2; idx_k_2++){
            double k_vec[3];
            k_vec[0] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[0] - ax;
            k_vec[1] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[1] - ay;
            k_vec[2] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[2] - az;
            double coefs_d[3];
            solve_3by3(reciprocal_lattice_vector, k_vec, coefs_d);

            //coefs_d[0] *= N1;
            //coefs_d[1] *= N2;
            //coefs_d[2] *= N3;
            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                //heff_comm_k[iorb] =_heff_k[idx_k_1*N2 + idx_k_2][iorb]; 
                heff_comm_k[iorb] = bicubicInterpolate(_heff_k, N1, N2, coefs_d[0], coefs_d[1], iorb); 
                rho_comm_k[iorb] = _rho_k[idx_k_1*N2 + idx_k_2][iorb] + 0.5*_dt*_k2[idx_k_1*N2 + idx_k_2][iorb]; 
            } 
            _calc_commutator(heff_comm_k, rho_comm_k, comm_k, _num_orbitals);
            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                _k3[idx_k_1*N2 + idx_k_2][iorb] = cdouble(0.0,-1.0)*comm_k[iorb];
            }
        }
    } 
    delete[] heff_comm_k;
    delete[] rho_comm_k;
    delete[] comm_k;
}

void Solver_kspace::_update_k4(const double ex, const double ey, const double ez,
                               const double ax, const double ay, const double az){
    int N1 = _settings->nr1;
    int N2 = _settings->nr2;
    int N3 = _settings->nr3;
    cdouble* heff_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* rho_comm_k = new cdouble[_num_orbitals*_num_orbitals];
    cdouble* comm_k = new cdouble[_num_orbitals*_num_orbitals];
    double **reciprocal_lattice_vector;
    reciprocal_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){reciprocal_lattice_vector[i] = new double[3];}
    _grid->reciprocal_vector(reciprocal_lattice_vector);
    for(int idx_k_1 = 0; idx_k_1 < _settings->nr1; idx_k_1++){
        for(int idx_k_2 = 0; idx_k_2 < _settings->nr2; idx_k_2++){
            double k_vec[3];
            k_vec[0] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[0] - ax;
            k_vec[1] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[1] - ay;
            k_vec[2] = _grid->Kvecs(idx_k_1*N2 + idx_k_2)[2] - az;
            double coefs_d[3];
            solve_3by3(reciprocal_lattice_vector, k_vec, coefs_d);
            //coefs_d[0] *= N1;
            //coefs_d[1] *= N2;
            //coefs_d[2] *= N3;

            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                //heff_comm_k[iorb] =_heff_k[idx_k_1*N2 + idx_k_2][iorb]; 
                heff_comm_k[iorb] = bicubicInterpolate(_heff_k, N1, N2, coefs_d[0], coefs_d[1], iorb); 
                rho_comm_k[iorb] = _rho_k[idx_k_1*N2 + idx_k_2][iorb] + _dt*_k3[idx_k_1*N2 + idx_k_2][iorb]; 
            } 
            _calc_commutator(heff_comm_k, rho_comm_k, comm_k, _num_orbitals);
            for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
                _k4[idx_k_1*N2 + idx_k_2][iorb] = cdouble(0.0,-1.0)*comm_k[iorb];
            }
        }
    } 
    delete[] heff_comm_k;
    delete[] rho_comm_k;
    delete[] comm_k;
}

void Solver_kspace::_calc_commutator(cdouble *A, cdouble *B, cdouble *C, const int n){
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
                    C[iorb*_num_orbitals + jorb] = value_comm; 
                }
            }
}

void Solver_kspace::_update_heff(const double ex, const double ey, const double ez,
                          const double ax, const double ay, const double az){
    cdouble heff_value = {0.0,0.0};
    for(int idx_k=0; idx_k<_num_points; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
            heff_value =  _h0_k[idx_k][iorb] 
                        + ex*_xbc_k[idx_k][iorb] 
                        + ey*_ybc_k[idx_k][iorb] 
                        + ez*_zbc_k[idx_k][iorb];
            _heff_k[idx_k][iorb] = heff_value;
        }
    }
}


void Solver_kspace::_clear_kn(){
    for(int idx_k=0; idx_k < _num_points; idx_k++){
        for (int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            _k1[idx_k][iorb] = {0.0,0.0};
            _k2[idx_k][iorb] = {0.0,0.0};
            _k3[idx_k][iorb] = {0.0,0.0};
            _k4[idx_k][iorb] = {0.0,0.0};
        }
    }
}
