#include <iostream>
#include <iomanip>
#include <lapacke.h>
#include <string>
#include "swe.h"
#include "observable.h"
#include "fftw_helper.h"

SWESim::SWESim(){}

SWESim::SWESim(const std::string &path_tb, Settings* settings){
    _path_tb = path_tb;
    _unit_system = SI;
    _settings = settings;
    _init();        
}

void SWESim::restart(){
    delete _rho;
    delete _r_bc[0];
    delete _r_bc[1];
    delete _r_bc[2];
    delete _v[0];
    delete _v[1];
    delete _v[2];
    delete _diagonalization;
    delete _solver;
    delete _wannier;
    delete _grid;
    delete _efield;
    delete _jx;
    delete _jy;
    delete _jz;
    _init();
}

void SWESim::init(){
    _init();
}

void SWESim::set_path_tb(const std::string &path_tb){
    _path_tb = path_tb;
}

void SWESim::set_settings(Settings* settings){
    _settings = settings;
}

void SWESim::_init(){

    _wannier = new WannierTB(_path_tb);
    _settings->nt = _settings->tmax/_settings->dt;
    _settings->path_tb = _path_tb;
    _settings->num_orb = _wannier->num_orb;
    _settings->num_sites = _wannier->num_sites;

    _grid = new Grid(_settings, _wannier);
    _hamiltonian = new Hamiltonian(_settings, _grid, _wannier, WGAUGE, RSPACE);
    _rho = new RDM(_settings, _grid, _wannier, WGAUGE, RSPACE);
    _r_bc[0] = new BerryConnection(_settings, _grid, _wannier, WGAUGE, RSPACE , XAXIS);
    _r_bc[1] = new BerryConnection(_settings, _grid, _wannier, WGAUGE, RSPACE , YAXIS);
    _r_bc[2] = new BerryConnection(_settings, _grid, _wannier, WGAUGE, RSPACE , ZAXIS);
    _v[0] = new Velocity(_settings, _grid, _wannier, WGAUGE, RSPACE);
    _v[1] = new Velocity(_settings, _grid, _wannier, WGAUGE, RSPACE);
    _v[2] = new Velocity(_settings, _grid, _wannier, WGAUGE, RSPACE);
    _v[0]->setup(_hamiltonian, _r_bc[0], XAXIS);
    _v[1]->setup(_hamiltonian, _r_bc[1], YAXIS);
    _v[2]->setup(_hamiltonian, _r_bc[2], ZAXIS);
    _diagonalization = new MatrixField(_settings, _wannier);
    _efield = new Efield(_settings, _grid, _unit_system);

    _jx = new Observable(_settings, _grid, _rho, _v[0]);
    _jy = new Observable(_settings, _grid, _rho, _v[1]);
    _jz = new Observable(_settings, _grid, _rho, _v[2]);

    _solver = new Solver(_settings, _grid, _hamiltonian, _rho, _r_bc, _efield, _wannier);

    _convert_to_au(); 
    
    _efield->init_fields();
    _hamiltonian->convert_to_k();
    _hamiltonian->diagonalization_W_KSPACE(_diagonalization);
    _rho->convert_to_k();
    _hamiltonian->convert_w_to_h(_diagonalization);
    _rho->convert_w_to_h(_diagonalization);
    _rho->calculate_equilibrium(_hamiltonian);
    _hamiltonian->convert_h_to_w(_diagonalization);
    _hamiltonian->convert_to_r();
    _rho->convert_h_to_w(_diagonalization);
    _rho->convert_to_r();
    
    _solver->init();
}

void SWESim::_convert_to_au(){
    if (_unit_system == AU){
        return;
    }
    _wannier->convert_to_au();
    _grid->convert_to_au();
    _hamiltonian->convert_to_au();
    _rho->convert_to_au();
    _r_bc[0]->convert_to_au();
    _r_bc[1]->convert_to_au();
    _r_bc[2]->convert_to_au();
    _efield->convert_to_au();
    _v[0]->convert_to_au();
    _v[1]->convert_to_au();
    _v[2]->convert_to_au();
    _unit_system = AU;
}

void SWESim::_convert_to_si(){
    if (_unit_system == SI){
        return;
    }
    _wannier->convert_to_si();
    _grid->convert_to_si();
    _hamiltonian->convert_to_si();
    _rho->convert_to_si();
    _r_bc[0]->convert_to_si();
    _r_bc[1]->convert_to_si();
    _r_bc[2]->convert_to_si();
    _efield->convert_to_si();
    _v[0]->convert_to_si();
    _v[1]->convert_to_si();
    _v[2]->convert_to_si();
    _unit_system = SI;
}

void SWESim::test_files(){
    std::cout<<"Writing test files"<<std::endl;
    _convert_to_si();
    _grid->supercell_to_file("test/Rvecs_SI.dat");
    _grid->reciprocal_to_file("test/Kvecs_SI.dat");
    _hamiltonian->write("test/hamiltonian_SI.dat");

    _convert_to_au();
    _grid->supercell_to_file("test/Rvecs_AU.dat");
    _grid->reciprocal_to_file("test/Kvecs_AU.dat");
    _hamiltonian->write("test/hamiltonian_AU.dat");

    std::cout<<"Checking direct to reciprocal transformation"<<std::endl;
    _hamiltonian->convert_to_k();
    _hamiltonian->write("test/hamiltonian_k_AU.dat");
    
    std::cout<<"Checking reciprocal to direct transformation"<<std::endl;
    _hamiltonian->convert_to_r();
    _hamiltonian->write("test/hamiltonian_r_AU.dat");

    std::cout<<"Calculating Hamiltonian U matrix"<<std::endl;
    _hamiltonian->convert_to_k();
    _hamiltonian->diagonalization_W_KSPACE(_diagonalization);
    
    _diagonalization->write_to_file("test/Uk_AU.dat");

    _hamiltonian->convert_w_to_h(_diagonalization);
    _hamiltonian->write("test/hamiltonian_k_H_AU.dat");

    _hamiltonian->convert_h_to_w(_diagonalization);
    _hamiltonian->write("test/hamiltonian_k_W_AU.dat");

    _hamiltonian->convert_w_to_h(_diagonalization);
    _rho->convert_to_k();
    _rho->convert_w_to_h(_diagonalization);
    
    _rho->calculate_equilibrium(_hamiltonian);
    _rho->write("test/rho_k_H_AU.dat");
    _rho->convert_h_to_w(_diagonalization);
    _rho->write("test/rho_k_W_AU.dat");
    _rho->convert_to_r();
    _rho->write("test/rho_r_W_AU.dat");

    _efield->init_fields();
    _efield->write();

    _r_bc[0]->write("test/x_bc.dat");
    _r_bc[1]->write("test/y_bc.dat");
    _r_bc[2]->write("test/z_bc.dat");

}

void SWESim::run_simulation(){
    int _num_points  =_settings->nr1*_settings->nr2*_settings->nr3;
    int _num_orbitals = _settings->num_orb;
    cdouble* peierls_phase = new cdouble[_num_points];
    //Observable rho(_settings, _grid, _rho, _rho);
    for(int ti = 0; ti<_settings->nt; ti++){
        if(ti%10 == 0){
            std::cout<<"\t Progress: " << ((double)ti/(double)_settings->nt)*100<<std::setprecision(3)<<" %\r"<<std::flush;
        }
        double ax = _efield->A_x[ti];
        double ay = _efield->A_y[ti];
        double az = _efield->A_z[ti];
        _calc_peierls_phase(ax, ay, az, peierls_phase);
        _solver->step_rk4(ti, peierls_phase);
        for(int idx_r = 0; idx_r<_num_points; idx_r ++){
            for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
                cdouble rho_value = _rho->data_ptr()[idx_r][iorb];
                _rho->set(rho_value*std::conj(peierls_phase[idx_r]), idx_r, iorb);
            }
        }
        _jx->calculate();
        _jy->calculate();
        _jz->calculate();
        //rho.calculate();
        //if(ti%100 == 0){
        //    std::ostringstream oss;
        //    oss << "test/rho_r_W_AU_" << ti << ".dat"; 
         //    _rho->write(oss.str());
        //}
        for(int idx_r = 0; idx_r<_num_points; idx_r ++){
            for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
                cdouble rho_value = _rho->data_ptr()[idx_r][iorb];
                _rho->set(rho_value*peierls_phase[idx_r], idx_r, iorb);
            }
        }
    }
    //oss.str("");
    //oss << _settings->path_results<<"/rho.dat";
    //rho.write(oss.str());
    delete[] peierls_phase;
}

void SWESim::_calc_peierls_phase(double ax, double ay, double az, cdouble* peierls_phase){
    int _num_points  =_settings->nr1*_settings->nr2*_settings->nr3;
    for(int idx_r = 0; idx_r<_num_points; idx_r++){
        cdouble phase = std::exp(cdouble(0.0,+ax*_grid->Rvecs(idx_r)[0]
                                             +ay*_grid->Rvecs(idx_r)[1]
                                             +az*_grid->Rvecs(idx_r)[2]));
        peierls_phase[idx_r] = phase;
    }
}

void SWESim::get_current(double* time, cdouble* jx, cdouble* jy, cdouble* jz){
    for(int ti = 0; ti < _settings->nt; ti++){
        time[ti] = _grid->t()[ti];
        jx[ti] = _jx->get_ptr()[ti];
        jy[ti] = _jy->get_ptr()[ti];
        jz[ti] = _jz->get_ptr()[ti];
    }    
}

void SWESim::save_current(const std::string &path){
    std::ostringstream oss;
    oss << path<<"/jx.dat";
    _jx->write(oss.str());
    oss.str("");
    oss << path<<"/jy.dat";
    _jy->write(oss.str());
    oss.str("");
    oss << path<<"/jz.dat";
    _jz->write(oss.str());

}

SWESim::~SWESim(){
    delete _rho;
    delete _r_bc[0];
    delete _r_bc[1];
    delete _r_bc[2];
    delete _v[0];
    delete _v[1];
    delete _v[2];
    delete _diagonalization;
    delete _solver;
    delete _wannier;
    delete _grid;    
    delete _efield;
    delete _jx;
    delete _jy;
    delete _jz;
}
