#include <iostream>
#include "berry_connection.h"
#include "vec3_util.h"
#include "constants.h"


void BerryConnection::_setup(){
    //std::cout<<"Setting Hamiltonian"<<std::endl;
    double **temp_lattice_vector;
    temp_lattice_vector = new double* [3];

    for (int i=0; i<3; i++){
        temp_lattice_vector[i] = new double[3];
    }

    _grid->lattice_vector(temp_lattice_vector);

    cdouble *wann_berry_connection;
    switch(_axis){
        case XAXIS:
            wann_berry_connection = _wannier->x_bc;
            break;
        case YAXIS:
            wann_berry_connection = _wannier->y_bc;
            break;
        case ZAXIS:
            wann_berry_connection = _wannier->z_bc;
            break;
    }
    for(int idx_r=0; idx_r < _numpoints; idx_r++){
        double coefs_d[3];
        int coefs_i[3];
        solve_3by3(temp_lattice_vector, _grid->Rvecs(idx_r), coefs_d);
        coefs_i[0] = (int)coefs_d[0];
        coefs_i[1] = (int)coefs_d[1];
        coefs_i[2] = (int)coefs_d[2];
        int in_wannier = _wannier->check_vec_in_wannier(coefs_i);
        if (in_wannier != -1){
            for(int idx_orb = 0; idx_orb < _num_orbitals * _num_orbitals; idx_orb++){
                _matrix->set(-wann_berry_connection[in_wannier*_num_orbitals*_num_orbitals + idx_orb], idx_r, idx_orb);
            }
        }
    }
    _gauge = WGAUGE; _space_type = RSPACE;
    
    for(int i=0; i<3; i++){
        delete[] temp_lattice_vector[i];
    }
    delete[] temp_lattice_vector;
}


void BerryConnection::convert_to_au(){
    if (_unit_system == AU){
        return;
    }
    for (int idx_r = 0; idx_r < _numpoints; idx_r++){    
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            cdouble value = _matrix->get(idx_r, iorb);
            _matrix->set( value*length_A2au, idx_r, iorb);
        }
    }    
    _unit_system = AU;
    _matrix->set_unit_system(_unit_system);
}

void BerryConnection::convert_to_si(){
    if (_unit_system == SI){
        return;
    }
    for (int idx_r = 0; idx_r < _numpoints; idx_r++){    
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            cdouble value = _matrix->get(idx_r, iorb);
            _matrix->set(value*length_au2A, idx_r, iorb);
        }
    }    
    _unit_system = SI;
    _matrix->set_unit_system(_unit_system);
}

int BerryConnection::axis(){
    return _axis;
}
