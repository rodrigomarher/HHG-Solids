#include <iostream>
#include <lapacke.h>
#include <string>
#include "swe.h"
#include "observable.h"
#include "fftw_helper.h"

SWESim::SWESim(const std::string &path_tb){
    _path_tb = path_tb;
    _init();        
}

void SWESim::_init(){
    _unit_system = SI;
    _wannier = new WannierTB(_path_tb);


    _settings = new Settings();
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

    _solver = new Solver(_settings, _grid, _hamiltonian, _rho, _r_bc, _efield, _wannier);

    _convert_to_au(); 
    
    _efield->init_fields();
    _hamiltonian->write("test/hamiltonian_r_W_AU.dat");
    _hamiltonian->convert_to_k();
    _hamiltonian->write("test/hamiltonian_k_W_AU.dat");
    _hamiltonian->diagonalization_W_KSPACE(_diagonalization);
    _rho->convert_to_k();
    _hamiltonian->convert_w_to_h(_diagonalization);
    _rho->convert_w_to_h(_diagonalization);
    _rho->calculate_equilibrium(_hamiltonian);
    _hamiltonian->convert_h_to_w(_diagonalization);
    _hamiltonian->write("test/hamiltonian_k_W_AU.dat");
    _hamiltonian->convert_to_r();
    _rho->convert_h_to_w(_diagonalization);
    _rho->write("test/rho_k_W_AU.dat");
    _rho->convert_to_r();
    
    _grid->supercell_to_file("test/Rvecs_AU.dat");
    _grid->reciprocal_to_file("test/Kvecs_AU.dat");
    _rho->write("test/rho_r_W_AU.dat");
    _r_bc[0]->write("test/rx_W_AU.dat");
    _r_bc[1]->write("test/ry_W_AU.dat");
    _r_bc[2]->write("test/rz_W_AU.dat");
    _v[0]->write("test/vx_W_AU.dat");
    _v[1]->write("test/vy_W_AU.dat");
    _v[2]->write("test/vz_W_AU.dat");
    _efield->write();
    _diagonalization->write_to_file("test/U_k_AU.dat");
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
    Observable rho(_settings, _grid, _rho, _rho);
    Observable jx(_settings, _grid, _rho, _v[0]);
    Observable jy(_settings, _grid, _rho, _v[1]);
    Observable jz(_settings, _grid, _rho, _v[2]);
    for(int ti = 0; ti<_settings->nt; ti++){
        std::cout<<"\t ti: " << ti<<std::endl;
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
        jx.calculate();
        jy.calculate();
        jz.calculate();
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
    std::ostringstream oss;
    oss << _settings->path_results<<"/jx.dat";
    jx.write(oss.str());
    oss.str("");
    oss << _settings->path_results<<"/jy.dat";
    jy.write(oss.str());
    oss.str("");
    oss << _settings->path_results<<"/jz.dat";
    jz.write(oss.str());
    oss.str("");
    oss << _settings->path_results<<"/rho.dat";
    rho.write(oss.str());
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

SWESim::~SWESim(){
    delete _wannier;
    delete _settings;
    delete _grid;
    delete _hamiltonian;
    delete _diagonalization;
    delete _r_bc[0];
    delete _r_bc[1];
    delete _r_bc[2];
}
