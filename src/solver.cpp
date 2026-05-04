#include <iostream>
#include "cblas.h"
#include "solver.h"
#include "vec3_util.h"
#include "fftw_helper.h"

#define MOD(a, b) (((a) % (b) + (b)) % (b))

double _interpolate(const int n, const double theta, double* arr, const int nmax){
    if(n>nmax-1){return arr[nmax-1];}
    return (1.0-theta)*arr[n] + theta*arr[n-1];
}

Solver::Solver(){

}

void Solver::init(){
    return;
}

Solver::Solver(Settings *settings, Grid* grid, Hamiltonian* hamiltonian, RDM* rho, BerryConnection** r_bc, Efield* efield, WannierTB* wannier){
    _settings = settings; 
    _grid = grid;
    _hamiltonian = hamiltonian;
    _rho = rho;
    _efield = efield;
    _wannier = wannier;
    _rho_aux = new RDM(settings, grid, wannier, rho->gauge(), rho->space_type());
    _rho->copy_data(_rho_aux);
    _r_bc = r_bc;

    _num_points = _settings->nr1 * _settings->nr2 * _settings->nr3;
    _num_orbitals = _settings->num_orb;

    _allocate();
    
}

void Solver::_allocate(){
    _k1 = new cdouble*[_num_points];
    _k2 = new cdouble*[_num_points];
    _k3 = new cdouble*[_num_points];
    _k4 = new cdouble*[_num_points];
    _heff = new cdouble*[_num_points];

    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        _k1[idx_r] = new cdouble[_num_orbitals*_num_orbitals];
        _k2[idx_r] = new cdouble[_num_orbitals*_num_orbitals];
        _k3[idx_r] = new cdouble[_num_orbitals*_num_orbitals];
        _k4[idx_r] = new cdouble[_num_orbitals*_num_orbitals];
        _heff[idx_r] = new cdouble[_num_orbitals*_num_orbitals];
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            _k1[idx_r][iorb] = cdouble(0.0,0.0);
            _k2[idx_r][iorb] = cdouble(0.0,0.0);
            _k3[idx_r][iorb] = cdouble(0.0,0.0);
            _k4[idx_r][iorb] = cdouble(0.0,0.0);
            _heff[idx_r][iorb] = cdouble(0.0,0.0);
        }
    }
}

void Solver::_deallocate(){
    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        delete[] _k1[idx_r];
        delete[] _k2[idx_r];
        delete[] _k3[idx_r];
        delete[] _k4[idx_r];
        delete[] _heff[idx_r];
    }
    delete[] _k1;
    delete[] _k2;
    delete[] _k3;
    delete[] _k4;
    delete[] _heff;
}

Solver::~Solver(){
    delete _rho_aux;

    _deallocate();
}



void Solver::step_rk4(const int ti, cdouble* peierls_phase){
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
    //_update_k1_conv_fftw(ex,  ey, ez, ax, ay, az); 
   // _update_heff(ex_dt2, ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2);
    //_update_k2_conv_fftw(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    //_update_k3_conv_fftw(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
   // _update_heff(ex_dt, ey_dt, ez_dt, ax_dt, ay_dt, az_dt);
    //_update_k4_conv_fftw(ex_dt,  ey_dt, ez_dt, ax_dt, ay_dt, az_dt); 
    //_update_heff(ex, ey, ez, ax, ay, az);
    _update_heff(ex, ey, ez, ax, ay, az);
    _update_k1_conv_fftw(ex,  ey, ez, ax, ay, az); 
    _update_heff(ex_dt2, ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2);
    _update_k2_conv_fftw(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    _update_k3_conv_fftw(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    _update_heff(ex_dt, ey_dt, ez_dt, ax_dt, ay_dt, az_dt);
    _update_k4_conv_fftw(ex_dt,  ey_dt, ez_dt, ax_dt, ay_dt, az_dt); 

    //_update_k2_no_blas(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    //_update_k3_no_blas(ex_dt2,  ey_dt2, ez_dt2, ax_dt2, ay_dt2, az_dt2); 
    //_update_k4_no_blas(ex_dt,  ey_dt, ez_dt, ax_dt, ay_dt, az_dt); 
    
    
    cdouble prefac = {0.1666666666666*_dt, 0.0};
    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
            cdouble value = _rho->get(idx_r, iorb);
            value = value + prefac*(_k1[idx_r][iorb] 
                                    + 2.0*_k2[idx_r][iorb]
                                    + 2.0*_k3[idx_r][iorb]
                                    + _k4[idx_r][iorb]);
            _rho->set(value, idx_r, iorb);
        }
    }
}

void Solver::_update_k1(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble *tmp_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_h0_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_xbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_ybc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_zbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_sum = new cdouble[_num_orbitals*_num_orbitals];
    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            tmp_sum[iorb] = cdouble(0.0,0.0);
        }
        for(int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
            //double r_diff[3];
            //double coefs_d[3];
            //int coefs_i[3];
            //r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
            //r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
            //r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
            //solve_3by3(temp_lattice_vector, r_diff, coefs_d);
            //coefs_i[0] = (int)coefs_d[0];
            //coefs_i[1] = (int)coefs_d[1];
            //coefs_i[2] = (int)coefs_d[2];
            //coefs_i[0] = MOD(((coefs_i[0] + _settings->nr1/2) ),_settings->nr1);
            //coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
            //coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
            //int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
            //               + coefs_i[1]*_settings->nr3
            //               + coefs_i[2]; 
            int idx_diff = idx_r_p;
            cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                              -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                              -1.0*az*_grid->Rvecs(idx_r_p)[2]));
            //cdouble *r_bc[3];
            cdouble *heff = _heff[idx_r_p];//_hamiltonian->data_ptr()[idx_r_p];
            cdouble *rho0 = _rho->data_ptr()[idx_diff];
            cblas_zcopy(_num_orbitals*_num_orbitals,
                    rho0,
                    1,
                    tmp_rho,
                    1);

            //r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
            //r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
            //r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

            cdouble alpha = {1.0, 0.0};
            alpha *= peierls_phase;
            cdouble beta = {0.0, 0.0};
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        heff,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);
            alpha = {-1.0, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasColMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        heff,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);

            // X - Electric field component
            //alpha = {ex, 0.0};
            //alpha *= peierls_phase;
            //beta = {0.0, 0.0};
            //cblas_zgemm(CblasRowMajor,
            //            CblasNoTrans,
            //            CblasNoTrans,
            //            _num_orbitals,
            //            _num_orbitals,
            //            _num_orbitals,
            //            &alpha,
            //            r_bc[0],
            //            _num_orbitals,
            //            tmp_rho,
            //            _num_orbitals,
            //            &beta,
            //            tmp_xbc_rho,
            //            _num_orbitals);
            //alpha = {-1.0*ex, 0.0};
            //alpha *= peierls_phase;
            //beta = {1.0, 0.0};
            //cblas_zgemm(CblasRowMajor,
            //            CblasNoTrans,
            //            CblasNoTrans,
            //            _num_orbitals,
            //            _num_orbitals,
            //            _num_orbitals,
            //            &alpha,
            //            tmp_rho,
            //            _num_orbitals,
            //            r_bc[0],
            //            _num_orbitals,
            //            &beta,
            //            tmp_xbc_rho,
            //            _num_orbitals);
            //// Y - Electric field component
            //alpha = {ey, 0.0};
            //alpha *= peierls_phase;
            //beta = {0.0, 0.0};
            //cblas_zgemm(CblasRowMajor,
            //            CblasNoTrans,
            //            CblasNoTrans,
            //            _num_orbitals,
            //            _num_orbitals,
            //            _num_orbitals,
            //            &alpha,
            //            r_bc[1],
            //            _num_orbitals,
            //            tmp_rho,
            //            _num_orbitals,
            //            &beta,
            //            tmp_ybc_rho,
            //            _num_orbitals);
            //alpha = {-1.0*ey, 0.0};
            //alpha *= peierls_phase;
            //beta = {1.0, 0.0};
            //cblas_zgemm(CblasRowMajor,
            //            CblasNoTrans,
            //            CblasNoTrans,
            //            _num_orbitals,
            //            _num_orbitals,
            //            _num_orbitals,
            //            &alpha,
            //            tmp_rho,
            //            _num_orbitals,
            //            r_bc[1],
            //            _num_orbitals,
            //            &beta,
            //            tmp_ybc_rho,
            //            _num_orbitals);
            //// Z - Electric field component
            //alpha = {ez, 0.0};
            //alpha *= peierls_phase;
            //beta = {0.0, 0.0};
            //cblas_zgemm(CblasRowMajor,
            //            CblasNoTrans,
            //            CblasNoTrans,
            //            _num_orbitals,
            //            _num_orbitals,
            //            _num_orbitals,
            //            &alpha,
            //            r_bc[2],
            //            _num_orbitals,
            //            tmp_rho,
            //            _num_orbitals,
            //            &beta,
            //            tmp_zbc_rho,
            //            _num_orbitals);
            //alpha = {-1.0*ez, 0.0};
            //alpha *= peierls_phase;
            //beta = {1.0, 0.0};
            //cblas_zgemm(CblasRowMajor,
            //            CblasNoTrans,
            //            CblasNoTrans,
            //            _num_orbitals,
            //            _num_orbitals,
            //            _num_orbitals,
            //            &alpha,
            //            tmp_rho,
            //            _num_orbitals,
            //            r_bc[2],
            //            _num_orbitals,
            //            &beta,
            //            tmp_zbc_rho,
            //            _num_orbitals);

            ////Add all intermediatsums
            //alpha = {1.0,0.0};
            //cblas_zaxpy(_num_orbitals*_num_orbitals,
            //            &alpha,
            //            tmp_h0_rho,
            //            1,
            //            tmp_sum,
            //            1);
            //cblas_zaxpy(_num_orbitals*_num_orbitals,
            //            &alpha,
            //            tmp_xbc_rho,
            //            1,
            //            tmp_sum,
            //            1);
            //cblas_zaxpy(_num_orbitals*_num_orbitals,
            //            &alpha,
            //            tmp_ybc_rho,
            //            1,
            //            tmp_sum,
            //            1);
            //cblas_zaxpy(_num_orbitals*_num_orbitals,
            //            &alpha,
            //            tmp_zbc_rho,
            //            1,
            //            tmp_sum,
            //            1);
        }
        cblas_zcopy(_num_orbitals*_num_orbitals,
                    tmp_h0_rho,
                    1,
                    _k1[idx_r],
                    1);
    }

    delete[] tmp_rho;
    delete[] tmp_h0_rho;
    delete[] tmp_xbc_rho;
    delete[] tmp_ybc_rho;
    delete[] tmp_zbc_rho;
    delete[] tmp_sum;
    for (int i=0; i<3; i++){
        delete[] temp_lattice_vector[i];
    }
    delete[] temp_lattice_vector;
}

void Solver::_update_k2(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble *tmp_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_h0_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_xbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_ybc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_zbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_sum = new cdouble[_num_orbitals*_num_orbitals];
    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            tmp_sum[iorb] = cdouble(0.0,0.0);
        }
        for(int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
            double r_diff[3];
            double coefs_d[3];
            int coefs_i[3];
            r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
            r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
            r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
            solve_3by3(temp_lattice_vector, r_diff, coefs_d);
            coefs_i[0] = (int)coefs_d[0];
            coefs_i[1] = (int)coefs_d[1];
            coefs_i[2] = (int)coefs_d[2];
            coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
            coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
            coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
            int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                           + coefs_i[1]*_settings->nr3
                           + coefs_i[2]; 

            cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                              -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                              -1.0*az*_grid->Rvecs(idx_r_p)[2]));
            cdouble *r_bc[3];
            cdouble *h0 = _hamiltonian->data_ptr()[idx_r_p];
            cdouble *rho0 = _rho->data_ptr()[idx_diff];

            cblas_zcopy(_num_orbitals*_num_orbitals,
                        rho0,
                        1,
                        tmp_rho,
                        1);
            cdouble alpha = {0.5*_dt, 0.0};
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                    &alpha,
                    _k1[idx_diff],
                    1,
                    tmp_rho,
                    1);

            r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
            r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
            r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

            alpha = {1.0, 0.0};
            alpha *= peierls_phase;
            cdouble beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        h0,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);
            alpha = {-1.0, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        h0,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);

            // X - Electric field component
            alpha = {ex, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[0],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_xbc_rho,
                        _num_orbitals);
            alpha = {-1.0*ex, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[0],
                        _num_orbitals,
                        &beta,
                        tmp_xbc_rho,
                        _num_orbitals);
            // Y - Electric field component
            alpha = {ey, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[1],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_ybc_rho,
                        _num_orbitals);
            alpha = {-1.0*ey, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[1],
                        _num_orbitals,
                        &beta,
                        tmp_ybc_rho,
                        _num_orbitals);
            // Z - Electric field component
            alpha = {ez, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[2],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_zbc_rho,
                        _num_orbitals);
            alpha = {-1.0*ez, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[2],
                        _num_orbitals,
                        &beta,
                        tmp_zbc_rho,
                        _num_orbitals);

            //Add all intermediatsums
            alpha = {1.0,0.0};
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_h0_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_xbc_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_ybc_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_zbc_rho,
                        1,
                        tmp_sum,
                        1);
        }
        cblas_zcopy(_num_orbitals*_num_orbitals,
                    tmp_sum,
                    1,
                    _k1[idx_r],
                    1);
    }

    delete[] tmp_rho;
    delete[] tmp_h0_rho;
    delete[] tmp_xbc_rho;
    delete[] tmp_ybc_rho;
    delete[] tmp_zbc_rho;
    delete[] tmp_sum;
    for (int i=0; i<3; i++){
        delete[] temp_lattice_vector[i];
    }
    delete[] temp_lattice_vector;
}

void Solver::_update_k3(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble *tmp_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_h0_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_xbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_ybc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_zbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_sum = new cdouble[_num_orbitals*_num_orbitals];
    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            tmp_sum[iorb] = cdouble(0.0,0.0);
        }
        for(int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
            double r_diff[3];
            double coefs_d[3];
            int coefs_i[3];
            r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
            r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
            r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
            solve_3by3(temp_lattice_vector, r_diff, coefs_d);
            coefs_i[0] = (int)coefs_d[0];
            coefs_i[1] = (int)coefs_d[1];
            coefs_i[2] = (int)coefs_d[2];
            coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
            coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
            coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
            int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                           + coefs_i[1]*_settings->nr3
                           + coefs_i[2]; 

            cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                              -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                              -1.0*az*_grid->Rvecs(idx_r_p)[2]));
            cdouble *r_bc[3];
            cdouble *h0 = _hamiltonian->data_ptr()[idx_r_p];
            cdouble *rho0 = _rho->data_ptr()[idx_diff];

            cblas_zcopy(_num_orbitals*_num_orbitals,
                        rho0,
                        1,
                        tmp_rho,
                        1);
            cdouble alpha = {0.5*_dt, 0.0};
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                    &alpha,
                    _k2[idx_diff],
                    1,
                    tmp_rho,
                    1);

            r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
            r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
            r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

            alpha = {1.0, 0.0};
            alpha *= peierls_phase;
            cdouble beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        h0,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);
            alpha = {-1.0, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        h0,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);

            // X - Electric field component
            alpha = {ex, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[0],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_xbc_rho,
                        _num_orbitals);
            alpha = {-1.0*ex, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[0],
                        _num_orbitals,
                        &beta,
                        tmp_xbc_rho,
                        _num_orbitals);
            // Y - Electric field component
            alpha = {ey, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[1],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_ybc_rho,
                        _num_orbitals);
            alpha = {-1.0*ey, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[1],
                        _num_orbitals,
                        &beta,
                        tmp_ybc_rho,
                        _num_orbitals);
            // Z - Electric field component
            alpha = {ez, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[2],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_zbc_rho,
                        _num_orbitals);
            alpha = {-1.0*ez, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[2],
                        _num_orbitals,
                        &beta,
                        tmp_zbc_rho,
                        _num_orbitals);

            //Add all intermediatsums
            alpha = {1.0,0.0};
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_h0_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_xbc_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_ybc_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_zbc_rho,
                        1,
                        tmp_sum,
                        1);
        }
        cblas_zcopy(_num_orbitals*_num_orbitals,
                    tmp_sum,
                    1,
                    _k1[idx_r],
                    1);
    }

    delete[] tmp_rho;
    delete[] tmp_h0_rho;
    delete[] tmp_xbc_rho;
    delete[] tmp_ybc_rho;
    delete[] tmp_zbc_rho;
    delete[] tmp_sum;
    for (int i=0; i<3; i++){
        delete[] temp_lattice_vector[i];
    }
    delete[] temp_lattice_vector;
}

void Solver::_update_k4(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble *tmp_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_h0_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_xbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_ybc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_zbc_rho = new cdouble[_num_orbitals*_num_orbitals];
    cdouble *tmp_sum = new cdouble[_num_orbitals*_num_orbitals];
    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for(int idx_r = 0; idx_r < _num_points; idx_r++){
        for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            tmp_sum[iorb] = cdouble(0.0,0.0);
        }
        for(int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
            double r_diff[3];
            double coefs_d[3];
            int coefs_i[3];
            r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
            r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
            r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
            solve_3by3(temp_lattice_vector, r_diff, coefs_d);
            coefs_i[0] = (int)coefs_d[0];
            coefs_i[1] = (int)coefs_d[1];
            coefs_i[2] = (int)coefs_d[2];
            coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
            coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
            coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
            int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                           + coefs_i[1]*_settings->nr3
                           + coefs_i[2]; 

            cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                              -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                              -1.0*az*_grid->Rvecs(idx_r_p)[2]));
            cdouble *r_bc[3];
            cdouble *h0 = _hamiltonian->data_ptr()[idx_r_p];
            cdouble *rho0 = _rho->data_ptr()[idx_diff];

            cblas_zcopy(_num_orbitals*_num_orbitals,
                        rho0,
                        1,
                        tmp_rho,
                        1);
            cdouble alpha = {1.0*_dt, 0.0};
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                    &alpha,
                    _k3[idx_diff],
                    1,
                    tmp_rho,
                    1);

            r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
            r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
            r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

            alpha = {1.0, 0.0};
            alpha *= peierls_phase;
            cdouble beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        h0,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);
            alpha = {-1.0, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        h0,
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_h0_rho,
                        _num_orbitals);

            // X - Electric field component
            alpha = {ex, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[0],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_xbc_rho,
                        _num_orbitals);
            alpha = {-1.0*ex, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[0],
                        _num_orbitals,
                        &beta,
                        tmp_xbc_rho,
                        _num_orbitals);
            // Y - Electric field component
            alpha = {ey, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[1],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_ybc_rho,
                        _num_orbitals);
            alpha = {-1.0*ey, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[1],
                        _num_orbitals,
                        &beta,
                        tmp_ybc_rho,
                        _num_orbitals);
            // Z - Electric field component
            alpha = {ez, 0.0};
            alpha *= peierls_phase;
            beta = {0.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        r_bc[2],
                        _num_orbitals,
                        tmp_rho,
                        _num_orbitals,
                        &beta,
                        tmp_zbc_rho,
                        _num_orbitals);
            alpha = {-1.0*ez, 0.0};
            alpha *= peierls_phase;
            beta = {1.0, 0.0};
            cblas_zgemm(CblasRowMajor,
                        CblasNoTrans,
                        CblasNoTrans,
                        _num_orbitals,
                        _num_orbitals,
                        _num_orbitals,
                        &alpha,
                        tmp_rho,
                        _num_orbitals,
                        r_bc[2],
                        _num_orbitals,
                        &beta,
                        tmp_zbc_rho,
                        _num_orbitals);

            //Add all intermediatsums
            alpha = {1.0,0.0};
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_h0_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_xbc_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_ybc_rho,
                        1,
                        tmp_sum,
                        1);
            cblas_zaxpy(_num_orbitals*_num_orbitals,
                        &alpha,
                        tmp_zbc_rho,
                        1,
                        tmp_sum,
                        1);
        }
        cblas_zcopy(_num_orbitals*_num_orbitals,
                    tmp_sum,
                    1,
                    _k1[idx_r],
                    1);
    }

    delete[] tmp_rho;
    delete[] tmp_h0_rho;
    delete[] tmp_xbc_rho;
    delete[] tmp_ybc_rho;
    delete[] tmp_zbc_rho;
    delete[] tmp_sum;
    for (int i=0; i<3; i++){
        delete[] temp_lattice_vector[i];
    }
    delete[] temp_lattice_vector;
}


void Solver::_update_k1_no_blas(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble** ham_op = _hamiltonian->data_ptr();
    cdouble** rho_op = _rho->data_ptr();
    cdouble** xbc_op = _r_bc[0]->data_ptr();
    cdouble** ybc_op = _r_bc[1]->data_ptr();
    cdouble** zbc_op = _r_bc[2]->data_ptr();
    cdouble* h0;
    cdouble* rho; 

    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);
    
    int counter = 0; 
    //#pragma omp parallel for schedule(static)
    for (int idx_r = 0; idx_r < _num_points; idx_r ++){
        for ( int iorb = 0; iorb < _num_orbitals; iorb++){
            for (int jorb = 0; jorb < _num_orbitals; jorb++){
                cdouble value_sum_r_p = {0.0,0.0};
                for (int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
                    double r_diff[3];
                    double coefs_d[3];
                    int coefs_i[3];
                    r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
                    r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
                    r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
                    solve_3by3(temp_lattice_vector, r_diff, coefs_d);
                    coefs_i[0] = (int)coefs_d[0];
                    coefs_i[1] = (int)coefs_d[1];
                    coefs_i[2] = (int)coefs_d[2];
                    coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
                    coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
                    coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
                    int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                                   + coefs_i[1]*_settings->nr3
                                   + coefs_i[2]; 


                    h0 = _heff[idx_r_p];
                    rho = _rho->data_ptr()[idx_r_p];

                    cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                                      -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                                      -1.0*az*_grid->Rvecs(idx_r_p)[2]));
                    cdouble value_comm = 0;
                    for (int korb = 0; korb<_num_orbitals; korb++){
                        value_comm += ( h0[iorb*_num_orbitals + korb]*(rho[korb*_num_orbitals + jorb])
                                     -(rho[iorb*_num_orbitals + korb])*h0[korb*_num_orbitals + jorb]);
                    }
                    value_sum_r_p += peierls_phase*value_comm; 
                }
                _k1[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*value_sum_r_p; 
            }
        }
    }
    
    std::cout<<"counter: "<<counter<<std::endl;
}

void Solver::_update_k2_no_blas(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble* h0;
    cdouble* rho;
    cdouble* kn;
    cdouble* r_bc[3];

    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for (int idx_r = 0; idx_r < _num_points; idx_r ++){
        for ( int iorb = 0; iorb < _num_orbitals; iorb++){
            for (int jorb = 0; jorb < _num_orbitals; jorb++){
                cdouble value_sum_r_p = {0.0,0.0};
                for (int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
                    double r_diff[3];
                    double coefs_d[3];
                    int coefs_i[3];
                    r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
                    r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
                    r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
                    solve_3by3(temp_lattice_vector, r_diff, coefs_d);
                    coefs_i[0] = (int)coefs_d[0];
                    coefs_i[1] = (int)coefs_d[1];
                    coefs_i[2] = (int)coefs_d[2];
                    coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
                    coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
                    coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
                    int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                                   + coefs_i[1]*_settings->nr3
                                   + coefs_i[2]; 


                    h0 = _heff[idx_r_p];
                    rho = _rho->data_ptr()[idx_r_p];
                    kn = _k1[idx_diff];
                    r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
                    r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
                    r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

                    cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                                      -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                                      -1.0*az*_grid->Rvecs(idx_r_p)[2]));
                    cdouble value_comm = 0;
                    for (int korb = 0; korb<_num_orbitals; korb++){
                        value_comm += ( h0[iorb*_num_orbitals + korb]*(rho[korb*_num_orbitals + jorb]+0.5*_dt*kn[korb*_num_orbitals + jorb])
                                     -(rho[iorb*_num_orbitals + korb]+0.5*_dt*kn[iorb*_num_orbitals + korb])*h0[korb*_num_orbitals + jorb]);
                    }
                    value_sum_r_p += peierls_phase*value_comm; 
                }
                _k2[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*value_sum_r_p; 
            }
        }
    }
}

void Solver::_update_k3_no_blas(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble* h0;
    cdouble* rho;
    cdouble* kn;
    cdouble* r_bc[3];

    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for (int idx_r = 0; idx_r < _num_points; idx_r ++){
        for ( int iorb = 0; iorb < _num_orbitals; iorb++){
            for (int jorb = 0; jorb < _num_orbitals; jorb++){
                cdouble value_sum_r_p = {0.0,0.0};
                for (int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
                    double r_diff[3];
                    double coefs_d[3];
                    int coefs_i[3];
                    r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
                    r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
                    r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
                    solve_3by3(temp_lattice_vector, r_diff, coefs_d);
                    coefs_i[0] = (int)coefs_d[0];
                    coefs_i[1] = (int)coefs_d[1];
                    coefs_i[2] = (int)coefs_d[2];
                    coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
                    coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
                    coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
                    int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                                   + coefs_i[1]*_settings->nr3
                                   + coefs_i[2]; 


                    h0 = _heff[idx_r_p];
                    rho = _rho->data_ptr()[idx_r_p];
                    kn = _k2[idx_diff];
                    r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
                    r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
                    r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

                    cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                                      -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                                      -1.0*az*_grid->Rvecs(idx_r_p)[2]));
                    cdouble value_comm = 0;
                    for (int korb = 0; korb<_num_orbitals; korb++){
                        value_comm +=    ( h0[iorb*_num_orbitals + korb]*(rho[korb*_num_orbitals + jorb]+0.5*_dt*kn[korb*_num_orbitals + jorb])
                                     -(rho[iorb*_num_orbitals + korb]+0.5*_dt*kn[iorb*_num_orbitals + korb])*h0[korb*_num_orbitals + jorb]);
                    }
                    value_sum_r_p += peierls_phase*value_comm; 
                }
                _k3[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*value_sum_r_p; 
            }
        }
    }
}

void Solver::_update_k4_no_blas(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az){
    cdouble* h0;
    cdouble* rho;
    cdouble* kn;
    cdouble* r_bc[3];

    double **temp_lattice_vector;
    temp_lattice_vector = new double*[3];
    for (int i=0; i<3; i++){temp_lattice_vector[i] = new double[3];}
    _grid->lattice_vector(temp_lattice_vector);

    for (int idx_r = 0; idx_r < _num_points; idx_r ++){
        for ( int iorb = 0; iorb < _num_orbitals; iorb++){
            for (int jorb = 0; jorb < _num_orbitals; jorb++){
                cdouble value_sum_r_p = {0.0,0.0};
                for (int idx_r_p = 0; idx_r_p < _num_points; idx_r_p++){
                    double r_diff[3];
                    double coefs_d[3];
                    int coefs_i[3];
                    r_diff[0] = _grid->Rvecs(idx_r)[0] - _grid->Rvecs(idx_r_p)[0];
                    r_diff[1] = _grid->Rvecs(idx_r)[1] - _grid->Rvecs(idx_r_p)[1];
                    r_diff[2] = _grid->Rvecs(idx_r)[2] - _grid->Rvecs(idx_r_p)[2];
                    solve_3by3(temp_lattice_vector, r_diff, coefs_d);
                    coefs_i[0] = (int)coefs_d[0];
                    coefs_i[1] = (int)coefs_d[1];
                    coefs_i[2] = (int)coefs_d[2];
                    coefs_i[0] = MOD((coefs_i[0] + _settings->nr1/2),_settings->nr1);
                    coefs_i[1] = MOD((coefs_i[1] + _settings->nr2/2),_settings->nr2);
                    coefs_i[2] = MOD((coefs_i[2] + _settings->nr3/2),_settings->nr3);
                    int idx_diff =   coefs_i[0]*_settings->nr2*_settings->nr3
                                   + coefs_i[1]*_settings->nr3
                                   + coefs_i[2]; 


                    h0 = _heff[idx_r_p];
                    rho = _rho->data_ptr()[idx_r_p];
                    kn = _k3[idx_diff];
                    r_bc[0] = _r_bc[0]->data_ptr()[idx_r_p]; 
                    r_bc[1] = _r_bc[1]->data_ptr()[idx_r_p]; 
                    r_bc[2] = _r_bc[2]->data_ptr()[idx_r_p]; 

                    cdouble peierls_phase = std::exp(-1.0*cdouble(0.0,-1.0*ax*_grid->Rvecs(idx_r_p)[0]
                                                                      -1.0*ay*_grid->Rvecs(idx_r_p)[1]
                                                                      -1.0*az*_grid->Rvecs(idx_r_p)[2]));
                    cdouble value_comm = 0;
                    for (int korb = 0; korb<_num_orbitals; korb++){
                        value_comm +=    ( h0[iorb*_num_orbitals + korb]*(rho[korb*_num_orbitals + jorb]+1.0*_dt*kn[korb*_num_orbitals + jorb])
                                     -(rho[iorb*_num_orbitals + korb]+1.0*_dt*kn[iorb*_num_orbitals + korb])*h0[korb*_num_orbitals + jorb]);
                    }
                    value_sum_r_p += peierls_phase*value_comm; 
                }
                _k4[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*value_sum_r_p; 
            }
        }
    }
}

void Solver::_calc_commutator(cdouble *A, cdouble *B, cdouble *C, const int n){
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

void Solver::_update_heff(const double ex, const double ey, const double ez,
                          const double ax, const double ay, const double az){
    cdouble **h0 = _hamiltonian->data_ptr();
    cdouble **xbc = _r_bc[0]->data_ptr();
    cdouble **ybc = _r_bc[1]->data_ptr();
    cdouble **zbc = _r_bc[2]->data_ptr();

    cdouble heff_value = {0.0,0.0};
    for(int idx_r=0; idx_r<_num_points; idx_r++){
        for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
            heff_value = h0[idx_r][iorb] + ex*xbc[idx_r][iorb] + ey* ybc[idx_r][iorb] + ez*zbc[idx_r][iorb];
            _heff[idx_r][iorb] = heff_value;
        }
    }
}


void Solver::_update_k1_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *tmp_1 = new cdouble[_num_points];
    cdouble *tmp_2 = new cdouble[_num_points];
    cdouble **rho_ptr = _rho->data_ptr();

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for(int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
                tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
                tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb];
            }
            fftshift(tmp_1, _settings->nr1, _settings->nr2);
            fftshift(tmp_2, _settings->nr1, _settings->nr2);
            fft3(tmp_1, in, out, _num_points, forward);
            fft3(tmp_2, in, out, _num_points, forward);

            for(int idx_k = 0; idx_k < _num_points; idx_k++){
                heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
                rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
            }
        }
    }
 
    for(int idx_k=0; idx_k <_num_points; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                cdouble comm_value = {0.0, 0.0};
                for (int korb = 0; korb<_num_orbitals; korb++){
                    cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    
                    cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    comm_value += heff0*rho0 - rho1*heff1;
                }
                comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
            }
        }
    } 

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for (int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_k = 0; idx_k <_num_points; idx_k++){
                tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
            }

            ifft3(tmp_1, in, out, _num_points, backward);
            ifftshift(tmp_1, _settings->nr1, _settings->nr2);
            
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                _k1[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
            }
        }
    }

    fftw_free(in);
    fftw_free(out);
    delete[] heff_k;
    delete[] rho_k;    
    delete[] comm_k;
    delete[] tmp_1;
    delete[] tmp_2;
}

void Solver::_update_k2_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *tmp_1 = new cdouble[_num_points];
    cdouble *tmp_2 = new cdouble[_num_points];
    cdouble **rho_ptr = _rho->data_ptr();

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for(int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
                tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
                tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb] + 0.5*_dt*_k1[idx_r][iorb*_num_orbitals + jorb];
            }

            fftshift(tmp_1, _settings->nr1, _settings->nr2);
            fftshift(tmp_2, _settings->nr1, _settings->nr2);
            fft3(tmp_1, in, out, _num_points, forward);
            fft3(tmp_2, in, out, _num_points, forward);

            for(int idx_k = 0; idx_k < _num_points; idx_k++){
                heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
                rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
            }
        }
    }
 
    for(int idx_k=0; idx_k <_num_points; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                cdouble comm_value = {0.0, 0.0};
                for (int korb = 0; korb<_num_orbitals; korb++){
                    cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    
                    cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    comm_value += heff0*rho0 - rho1*heff1;
                }
                comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
            }
        }
    } 

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for (int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_k = 0; idx_k <_num_points; idx_k++){
                tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
            }

            ifft3(tmp_1, in, out, _num_points, backward);
            ifftshift(tmp_1, _settings->nr1, _settings->nr2);
            
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                _k2[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
            }
        }
    }

    fftw_free(in);
    fftw_free(out);
    delete[] heff_k;
    delete[] rho_k;    
    delete[] comm_k;
    delete[] tmp_1;
    delete[] tmp_2;
}

void Solver::_update_k3_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *tmp_1 = new cdouble[_num_points];
    cdouble *tmp_2 = new cdouble[_num_points];
    cdouble **rho_ptr = _rho->data_ptr();

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for(int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
                //cdouble peierls_phase = std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
                tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
                tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb] + 0.5*_dt*_k2[idx_r][iorb*_num_orbitals + jorb];
            }

            fftshift(tmp_1, _settings->nr1, _settings->nr2);
            fftshift(tmp_2, _settings->nr1, _settings->nr2);
            fft3(tmp_1, in, out, _num_points, forward);
            fft3(tmp_2, in, out, _num_points, forward);

            for(int idx_k = 0; idx_k < _num_points; idx_k++){
                heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
                rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
            }
        }
    }
 
    for(int idx_k=0; idx_k <_num_points; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                cdouble comm_value = {0.0, 0.0};
                for (int korb = 0; korb<_num_orbitals; korb++){
                    cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    
                    cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    comm_value += heff0*rho0 - rho1*heff1;
                }
                comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
            }
        }
    } 

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for (int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_k = 0; idx_k <_num_points; idx_k++){
                tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
            }

            ifft3(tmp_1, in, out, _num_points, backward);
            ifftshift(tmp_1, _settings->nr1, _settings->nr2);

            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                _k3[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
            }
        }
    }

    fftw_free(in);
    fftw_free(out);
    delete[] heff_k;
    delete[] rho_k;    
    delete[] comm_k;
    delete[] tmp_1;
    delete[] tmp_2;
}

void Solver::_update_k4_conv_fftw(const double ex, const double ey, const double ez,
                                  const double ax, const double ay, const double az){
    cdouble *heff_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *rho_k  = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_num_points];
    cdouble *tmp_1 = new cdouble[_num_points];
    cdouble *tmp_2 = new cdouble[_num_points];
    cdouble **rho_ptr = _rho->data_ptr();

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_num_points);
    fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for(int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                //cdouble peierls_phase = 1.0;//std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
                //cdouble peierls_phase = std::exp(cdouble(0.0,_grid->Rvecs(idx_r)[0]*ax + _grid->Rvecs(idx_r)[1]*ay + _grid->Rvecs(idx_r)[2]*az));
                tmp_1[idx_r] = _peierls_phase[idx_r]*_heff[idx_r][iorb*_num_orbitals + jorb];
                tmp_2[idx_r] = rho_ptr[idx_r][iorb*_num_orbitals + jorb] + _dt*_k3[idx_r][iorb*_num_orbitals + jorb];
            }
            fftshift(tmp_1, _settings->nr1, _settings->nr2);
            fftshift(tmp_2, _settings->nr1, _settings->nr2);
            fft3(tmp_1, in, out, _num_points, forward);
            fft3(tmp_2, in, out, _num_points, forward);

            for(int idx_k = 0; idx_k < _num_points; idx_k++){
                heff_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_1[idx_k];
                rho_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = tmp_2[idx_k];
            }
        }
    }
 
    for(int idx_k=0; idx_k <_num_points; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                cdouble comm_value = {0.0, 0.0};
                for (int korb = 0; korb<_num_orbitals; korb++){
                    cdouble heff0 = heff_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    cdouble heff1 = heff_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    
                    cdouble rho0 = rho_k[korb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
                    cdouble rho1 = rho_k[iorb*_num_orbitals*_num_points + korb*_num_points + idx_k];
                    comm_value += heff0*rho0 - rho1*heff1;
                }
                comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k] = comm_value;
            }
        }
    } 

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for (int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_k = 0; idx_k <_num_points; idx_k++){
                tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_num_points + jorb*_num_points + idx_k];
            }

            ifft3(tmp_1, in, out, _num_points, backward);
            ifftshift(tmp_1, _settings->nr1, _settings->nr2);
            
            for(int idx_r = 0; idx_r < _num_points; idx_r++){
                _k4[idx_r][iorb*_num_orbitals + jorb] = cdouble(0.0,-1.0)*tmp_1[idx_r];
            }
        }
    }

    fftw_free(in);
    fftw_free(out);
    delete[] heff_k;
    delete[] rho_k;    
    delete[] comm_k;
    delete[] tmp_1;
    delete[] tmp_2;
}

void Solver::_clear_kn(){
    for(int idx_r=0; idx_r < _num_points; idx_r++){
        for (int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            _k1[idx_r][iorb] = {0.0,0.0};
            _k2[idx_r][iorb] = {0.0,0.0};
            _k3[idx_r][iorb] = {0.0,0.0};
            _k4[idx_r][iorb] = {0.0,0.0};
        }
    }
}
