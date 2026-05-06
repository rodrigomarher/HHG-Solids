#include <iostream>
#include <fstream>
#include "settings.h"
#include "wannier_tb.h"
#include "grid.h"
#include "vec3_util.h"

Grid::Grid(){
}

Grid::Grid(Settings *settings, WannierTB *wannier){
    _settings = settings; _wannier = wannier;
    _unit_system = wannier->unit_system();
    for (int i=0; i<3; i++){
        _lattice_vector[i][0] = _wannier->lattice_vector[i][0];
        _lattice_vector[i][1] = _wannier->lattice_vector[i][1];
        _lattice_vector[i][2] = _wannier->lattice_vector[i][2];
    }
    _calculate_reciprocal_vectors();
    _allocate_grid();
    _setup_grid();
}
int Grid::unit_system(){
    return _unit_system;
}

void Grid::_calculate_reciprocal_vectors(){
    double tmp[3];
    cross(_lattice_vector[1], _lattice_vector[2], tmp);
    double V = dot(_lattice_vector[0], tmp);

    cross(_lattice_vector[1], _lattice_vector[2], tmp);
    _reciprocal_vector[0][0] = 2.0*M_PI/V*tmp[0];
    _reciprocal_vector[0][1] = 2.0*M_PI/V*tmp[1];
    _reciprocal_vector[0][2] = 2.0*M_PI/V*tmp[2];

    cross(_lattice_vector[2], _lattice_vector[0], tmp);
    _reciprocal_vector[1][0] = 2.0*M_PI/V*tmp[0];
    _reciprocal_vector[1][1] = 2.0*M_PI/V*tmp[1];
    _reciprocal_vector[1][2] = 2.0*M_PI/V*tmp[2];

    cross(_lattice_vector[0], _lattice_vector[1], tmp);
    _reciprocal_vector[2][0] = 2.0*M_PI/V*tmp[0];
    _reciprocal_vector[2][1] = 2.0*M_PI/V*tmp[1];
    _reciprocal_vector[2][2] = 2.0*M_PI/V*tmp[2];
}

void Grid::_allocate_grid(){
    _n1 = new int[_settings->nr1];
    _n2 = new int[_settings->nr2];
    _n3 = new int[_settings->nr3];
    _m1 = new int[_settings->nr1];
    _m2 = new int[_settings->nr2];
    _m3 = new int[_settings->nr3];
    _t  = new double[_settings->nt];
    _Rvecs = new double*[_settings->nr1*_settings->nr2*_settings->nr3];
    _Kvecs = new double*[_settings->nr1*_settings->nr2*_settings->nr3];
    for (int i=0; i<_settings->nr1*_settings->nr2*_settings->nr3; i++){
        _Rvecs[i] = new double[3];
        _Kvecs[i] = new double[3];
    }
}

void Grid::_setup_grid(){
    for(int i=0; i<_settings->nt; i++){
        _t[i] = (double)i*_settings->dt;
    }
    for(int i=0; i<_settings->nr1; i++){
        _n1[i] = -_settings->nr1/2 + (i);
        _m1[i] = i;
    }    
    for(int i=0; i<_settings->nr2; i++){
        _n2[i] = -_settings->nr2/2 + (i);
        _m2[i] = i;
    }    
    for(int i=0; i<_settings->nr3; i++){
        _n3[i] = -_settings->nr3/2 + (i);
        _m3[i] = i;
    }   

    for(int i=0; i<_settings->nr1; i++){
        for (int j=0; j<_settings->nr2; j++){
            for (int k=0 ; k<_settings->nr3; k++){
                double r_x = _n1[i]*_lattice_vector[0][0] +
                             _n2[j]*_lattice_vector[1][0] +
                             _n3[k]*_lattice_vector[2][0];
                double r_y = _n1[i]*_lattice_vector[0][1] +
                             _n2[j]*_lattice_vector[1][1] +
                             _n3[k]*_lattice_vector[2][1];
                double r_z = _n1[i]*_lattice_vector[0][2] +
                             _n2[j]*_lattice_vector[1][2] +
                             _n3[k]*_lattice_vector[2][2];
                double k_x = (double)_m1[i]/(double)_settings->nr1*_reciprocal_vector[0][0] +
                             (double)_m2[j]/(double)_settings->nr2*_reciprocal_vector[1][0] +
                             (double)_m3[k]/(double)_settings->nr3*_reciprocal_vector[2][0];
                double k_y = (double)_m1[i]/(double)_settings->nr1*_reciprocal_vector[0][1] +
                             (double)_m2[j]/(double)_settings->nr2*_reciprocal_vector[1][1] +
                             (double)_m3[k]/(double)_settings->nr3*_reciprocal_vector[2][1];
                double k_z = (double)_m1[i]/(double)_settings->nr1*_reciprocal_vector[0][2] +
                             (double)_m2[j]/(double)_settings->nr2*_reciprocal_vector[1][2] +
                             (double)_m3[k]/(double)_settings->nr3*_reciprocal_vector[2][2];
                _Rvecs[i*_settings->nr2*_settings->nr3 + j*_settings->nr3 + k][0] = r_x;
                _Rvecs[i*_settings->nr2*_settings->nr3 + j*_settings->nr3 + k][1] = r_y;
                _Rvecs[i*_settings->nr2*_settings->nr3 + j*_settings->nr3 + k][2] = r_z;
                _Kvecs[i*_settings->nr2*_settings->nr3 + j*_settings->nr3 + k][0] = k_x;
                _Kvecs[i*_settings->nr2*_settings->nr3 + j*_settings->nr3 + k][1] = k_y;
                _Kvecs[i*_settings->nr2*_settings->nr3 + j*_settings->nr3 + k][2] = k_z;
            }
        }
    } 
}

void Grid::convert_to_au(){
    if(_unit_system == AU){
        return;
    }
    for(int i=0; i<3; i++){
        _lattice_vector[i][0] *= length_A2au;
        _lattice_vector[i][1] *= length_A2au;
        _lattice_vector[i][2] *= length_A2au;
        _reciprocal_vector[i][0] /= length_A2au;
        _reciprocal_vector[i][1] /= length_A2au;
        _reciprocal_vector[i][2] /= length_A2au;
    }
    for (int i=0; i<_settings->nr1*_settings->nr2*_settings->nr3; i++){
        _Rvecs[i][0] *= length_A2au;
        _Rvecs[i][1] *= length_A2au;
        _Rvecs[i][2] *= length_A2au;
        _Kvecs[i][0] /= length_A2au;
        _Kvecs[i][1] /= length_A2au;
        _Kvecs[i][2] /= length_A2au;
    }
    for (int i=0; i<_settings->nt; i++){
        _t[i] *= time_fs2au;
    }
    _unit_system = AU;
}

void Grid::convert_to_si(){
    if(_unit_system == SI){
        return;
    }
    for(int i=0; i<3; i++){
        _lattice_vector[i][0] *= length_au2A;
        _lattice_vector[i][1] *= length_au2A;
        _lattice_vector[i][2] *= length_au2A;
        _reciprocal_vector[i][0] /= length_au2A;
        _reciprocal_vector[i][1] /= length_au2A;
        _reciprocal_vector[i][2] /= length_au2A;
    }
    for (int i=0; i<_settings->nr1*_settings->nr2*_settings->nr3; i++){
        _Rvecs[i][0] /= length_A2au;
        _Rvecs[i][1] /= length_A2au;
        _Rvecs[i][2] /= length_A2au;
        _Kvecs[i][0] *= length_A2au;
        _Kvecs[i][1] *= length_A2au;
        _Kvecs[i][2] *= length_A2au;
    }

    for (int i=0; i<_settings->nt; i++){
        _t[i] *= time_au2fs;
    }
    _unit_system = SI;
}
void Grid::print_lattice(){
    std::string unit_system_str = _unit_system == SI ? "(SI units)" : "(AU units)";
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<" Lattice vectors "<< unit_system_str << " : "<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"\t a1: {"<<_lattice_vector[0][0]<<"; "<<_lattice_vector[0][1] << "; "<<_lattice_vector[0][2]<<"}"<<std::endl;
    std::cout<<"\t a2: {"<<_lattice_vector[1][0]<<"; "<<_lattice_vector[1][1] << "; "<<_lattice_vector[1][2]<<"}"<<std::endl;
    std::cout<<"\t a2: {"<<_lattice_vector[2][0]<<"; "<<_lattice_vector[2][1] << "; "<<_lattice_vector[2][2]<<"}"<<std::endl;
    std::cout<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<" Reciprocal vectors "<< unit_system_str << " : "<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"\t b1: {"<<_reciprocal_vector[0][0]<<"; "<<_reciprocal_vector[0][1] << "; "<<_reciprocal_vector[0][2]<<"}"<<std::endl;
    std::cout<<"\t b2: {"<<_reciprocal_vector[1][0]<<"; "<<_reciprocal_vector[1][1] << "; "<<_reciprocal_vector[1][2]<<"}"<<std::endl;
    std::cout<<"\t b2: {"<<_reciprocal_vector[2][0]<<"; "<<_reciprocal_vector[2][1] << "; "<<_reciprocal_vector[2][2]<<"}"<<std::endl;

}

void Grid::supercell_to_file(const std::string &filename){
    std::string unit_system_str = _unit_system == SI ? "_SI_" : "_AU_";
    std::ofstream file (filename);
    int size = _settings->nr1*_settings->nr2*_settings->nr3;
    if (file.is_open()){
        file<<unit_system_str<<std::endl;
        for(int i = 0; i < size; i++){
            file << _Rvecs[i][0] << " " << _Rvecs[i][1] << " "<< _Rvecs[i][2] <<std::endl;;
        }
        file.close();
    }
    else{
        std::cout<< "[Grid] Error opening file: "<< filename<<std::endl;
    }
}

void Grid::reciprocal_to_file(const std::string &filename){
    std::string unit_system_str = _unit_system == SI ? "_SI_" : "_AU_";
    std::ofstream file (filename);
    int size = _settings->nr1*_settings->nr2*_settings->nr3;
    if (file.is_open()){
        file<<unit_system_str<<std::endl;
        for(int i = 0; i < size; i++){
            file << _Kvecs[i][0] << " " << _Kvecs[i][1] << " "<< _Kvecs[i][2] <<std::endl;;
        }
        file.close();
    }
    else{
        std::cout<< "[Grid] Error opening file: "<< filename<<std::endl;
    }
    
}

double* Grid::Rvecs(const int idx){ return _Rvecs[idx];}
double* Grid::Kvecs(const int idx){ return _Kvecs[idx];}
double  Grid::t(const int idx){return _t[idx];}
double** Grid::Rvecs(){return _Rvecs;}
double** Grid::Kvecs(){return _Kvecs;}
double* Grid::t(){return _t;}
void Grid::lattice_vector(double **output){
    output[0][0] = _lattice_vector[0][0]; 
    output[0][1] = _lattice_vector[0][1]; 
    output[0][2] = _lattice_vector[0][2]; 
    output[1][0] = _lattice_vector[1][0]; 
    output[1][1] = _lattice_vector[1][1]; 
    output[1][2] = _lattice_vector[1][2]; 
    output[2][0] = _lattice_vector[2][0]; 
    output[2][1] = _lattice_vector[2][1]; 
    output[2][2] = _lattice_vector[2][2]; 
}
void Grid::reciprocal_vector(double **output){
    output[0][0] = _reciprocal_vector[0][0]; 
    output[0][1] = _reciprocal_vector[0][1]; 
    output[0][2] = _reciprocal_vector[0][2]; 
    output[1][0] = _reciprocal_vector[1][0]; 
    output[1][1] = _reciprocal_vector[1][1]; 
    output[1][2] = _reciprocal_vector[1][2]; 
    output[2][0] = _reciprocal_vector[2][0]; 
    output[2][1] = _reciprocal_vector[2][1]; 
    output[2][2] = _reciprocal_vector[2][2]; 
}

Grid::~Grid(){
    delete[] _n1;
    delete[] _n2;
    delete[] _n3;
    delete[] _m1;
    delete[] _m2;
    delete[] _m3;
    for (int i = 0; i<_settings->nr1*_settings->nr2*_settings->nr3; i++){
        delete[] _Rvecs[i];
        delete[] _Kvecs[i];
    }
    delete[] _Rvecs;
    delete[] _Kvecs;
    delete[] _t;
}
