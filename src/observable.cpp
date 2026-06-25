#include <fstream>
#include <iostream>
#include "observable.h"

Observable::Observable(Settings_swe* settings_swe, Grid* grid, RDM* rho, Operator* op){
    _unit_system = grid->unit_system();
    _settings_swe = settings_swe;
    _grid = grid;
    _rho = rho;
    _operator = op;
    _num_points = settings_swe->nr1 * settings_swe->nr2 * settings_swe->nr3;
    _num_orbitals = settings_swe->num_orb;
    _data = new std::vector<cdouble>();
}

void Observable::calculate(){
    cdouble sum = {0.0,0.0};
    for (int idx_r_1 = 0; idx_r_1 < _settings_swe->nr1; idx_r_1++){
        for (int idx_r_2 = 0; idx_r_2 < _settings_swe->nr2; idx_r_2++){
            for(int iorb = 0; iorb < _num_orbitals; iorb++){
                for(int jorb = 0; jorb< _num_orbitals; jorb++){
                    cdouble value = _operator->data_ptr()[idx_r_1*_settings_swe->nr2 + idx_r_2][iorb*_num_orbitals + jorb];
                    value *= std::conj(_rho->data_ptr()[idx_r_1*_settings_swe->nr2 + idx_r_2][iorb*_num_orbitals + jorb]);
                    sum += value;
                }   
            }
        }
    }
    _data->push_back(sum);
}

void Observable::write(std::string filename){
    std::ofstream file(filename);
    if(file.is_open()){
        file << _unit_system << std::endl;
        for(int i = 0; i < _data->size(); i++){
            file << std::abs(_data->at(i)) << " "<< std::atan2(_data->at(i).imag(), _data->at(i).real()) <<std::endl;
        }
        file.close();
    }
    else{
        std::cout<< "[Observable] Error opening file: " << filename<< std::endl;
    }
}

cdouble* Observable::get_ptr(){
    return _data->data();
}

Observable::~Observable(){
    delete _data;
}
