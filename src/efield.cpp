#include <filesystem>
#include <string>
#include <fstream>
#include <iostream>
#include "efield.h"

Efield::Efield(Settings *settings, Grid* grid, int unit_system){
    _grid = grid; _settings = settings;
    _t = grid->t();
    _lambda = settings->lambda_nm*10.0;
    _omega = 2.0*M_PI*c_si/_lambda;
    _tmax_field = settings->tmax_field;
    _unit_system = unit_system; 
    _svec[0] = settings->sx;
    _svec[1] = settings->sy;
    _svec[2] = settings->sz;
    _allocate();

}

void Efield::_allocate(){
    E_x = new double[_settings->nt];
    E_y = new double[_settings->nt];
    E_z = new double[_settings->nt];
    A_x = new double[_settings->nt];
    A_y = new double[_settings->nt];
    A_z = new double[_settings->nt];
}

void Efield::init_fields(){
    if (_unit_system != AU){
        std::cout<<"[Efield::_init_fields] Error, unit system must be AU."<<std::endl;
        std::exit(1);
    }
    double efield_peak = sqrtf(_settings->intensity_wcm2*intensity_Wcm2au);
    double dt = _t[1]-_t[0];
    for(int i=0; i<_settings->nt; i++){
        E_x[i] = _settings->sx*efield_peak*_env_sin2(_t[i])*sinf(_omega*_t[i] + _settings->phi_x);
        E_y[i] = _settings->sy*efield_peak*_env_sin2(_t[i])*sinf(_omega*_t[i] + _settings->phi_y);
        E_z[i] = _settings->sz*efield_peak*_env_sin2(_t[i])*sinf(_omega*_t[i] + _settings->phi_z);
        A_x[i] = 0.0;
        A_y[i] = 0.0;
        A_z[i] = 0.0;
    }
    for(int i=1; i<_settings->nt; i++){
        A_x[i] = A_x[i-1] - E_x[i-1]*dt;
        A_y[i] = A_y[i-1] - E_y[i-1]*dt;
        A_z[i] = A_z[i-1] - E_z[i-1]*dt;
    }
}

double Efield::_env_sin2(double ti){
     if (ti<_tmax_field){
        return pow(sin(M_PI*ti/_tmax_field),2);
    }
    else {
        return 0.0;
    }

}

int Efield::unit_system(){
    return _unit_system;
}

void Efield::convert_to_au(){
    if(_unit_system == AU){
        return;
    }
    _lambda *= length_A2au;
    _tmax_field *= time_fs2au;
    _omega *= 1.0/(time_fs2au);
    _unit_system = AU;
    if(_unit_system != _grid->unit_system()){
        std::cout<<"[Warning] Efield units different from grid"<<std::endl;
    }
}
void Efield::convert_to_si(){
    if(_unit_system == SI){
        return;
    }
    _lambda *= length_au2A;
    _tmax_field *= time_au2fs;
    _omega *= 1.0/time_au2fs;
    _unit_system = SI;
    if(_unit_system != _grid->unit_system()){
        std::cout<<"[Warning] Efield units different from grid"<<std::endl;
    }
}

void Efield::write(){
    std::string path = "fields/";
    if(!std::filesystem::is_directory(path) || !std::filesystem::exists(path)){
        std::filesystem::create_directory(path);
    }
    
    std::ofstream file_efield( path + "efield.dat" );
    std::ofstream file_afield( path + "afield.dat" );

    if(file_efield.is_open() && file_afield.is_open()){
        file_efield << _unit_system<<std::endl;
        file_afield << _unit_system<<std::endl;
        for (int i=0; i<_settings->nt; i++){
            file_efield<<_t[i]<<" "<< E_x[i] << " "<< E_y[i]<<" "<< E_z[i]<<std::endl;
            file_afield<<_t[i]<<" "<< A_x[i] << " "<< A_y[i]<<" "<< A_z[i]<<std::endl;
        }
        file_efield.close();
        file_afield.close();
    }
    else{
        std::cout<< "[Efiel::write] Error opening files"<<std::endl;
    }
}


Efield::~Efield(){
    delete[] E_x;
    delete[] E_y;
    delete[] E_z;
    delete[] A_x;
    delete[] A_y;
    delete[] A_z;
}
