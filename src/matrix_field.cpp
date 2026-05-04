#include <string>
#include <iostream>
#include <fstream>

#include "matrix_field.h"

MatrixField::MatrixField (Settings *settings, WannierTB *wannier){
    _settings = settings;  _wannier = wannier;
    _numpoints = settings->nr1*_settings->nr2*_settings->nr3;
    _num_orbitals = settings->num_orb;
    _alloc_matrix();
    _unit_system = _wannier->unit_system();
    
}

void MatrixField::_alloc_matrix(){
    _data = new cdouble*[_numpoints];
    for(int i=0; i<_numpoints; i++){
        _data[i] = new cdouble[_num_orbitals*_num_orbitals];
    }
    for(int i=0; i<_numpoints; i++){
        for(int j=0; j<_num_orbitals*_num_orbitals; j++){
            _data[i][j] = cdouble(0.0,0.0);
        }
    }
}

void MatrixField::write_to_file(const std::string &filename){
    std::ofstream file (filename);
    std::string unit_system_str;
    switch(_unit_system){
        case 0:
            unit_system_str = "_SI_";
            break;
        case 1:
            unit_system_str = "_AU_";
            break;
        default:
            unit_system_str = "_ADIM_";
    }
    if (file.is_open()){
        file<<unit_system_str<<std::endl;
        for(int i = 0; i < _numpoints; i++){
            for (int iorb=0; iorb < _num_orbitals*_num_orbitals; iorb++){
                file<<i<<" " << std::abs(_data[i][iorb]) << " " << std::atan2(_data[i][iorb].imag(), _data[i][iorb].real())<<std::endl;
            }
        }
        file.close();
    }
    else{
        std::cout<< "[MatrixField] Error opening file: "<< filename<<std::endl;
    }
}

void MatrixField::set(cdouble value, const int idx, const int iorb){
    _data[idx][iorb] = value;
}

cdouble MatrixField::get(const int idx, const int iorb){
    return _data[idx][iorb];
}

cdouble& MatrixField::operator()(const int idx, const int iorb){
    return _data[idx][iorb];
}

void MatrixField::set_unit_system(const int value){
    _unit_system = value;
}

int MatrixField::get_unit_system(){
    return _unit_system;
}

cdouble** MatrixField::data_ptr(){
    return _data;
}
MatrixField::~MatrixField(){
    for(int i=0; i<_numpoints; i++){
        delete[] _data[i];
    }
    delete[] _data;
}
