#include <iostream>
#include "lapacke.h"
#include "hamiltonian.h"
#include "vec3_util.h"
#include "constants.h"


void Hamiltonian::_setup(){
    //std::cout<<"Setting Hamiltonian"<<std::endl;
    double **temp_lattice_vector;
    temp_lattice_vector = new double* [3];

    for (int i=0; i<3; i++){
        temp_lattice_vector[i] = new double[3];
    }

    _grid->lattice_vector(temp_lattice_vector);

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
                _matrix->set(_wannier->hopping[in_wannier*_num_orbitals*_num_orbitals + idx_orb], idx_r, idx_orb);
            }
        }
    }
    _gauge = WGAUGE; _space_type = RSPACE;
    
    for(int i=0; i<3; i++){
        delete[] temp_lattice_vector[i];
    }
    delete[] temp_lattice_vector;
}


void Hamiltonian::convert_to_au(){
    if (_unit_system == AU){
        return;
    }
    for (int idx_r = 0; idx_r < _numpoints; idx_r++){    
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            cdouble value = _matrix->get(idx_r, iorb);
            _matrix->set( value*energy_ev2au, idx_r, iorb);
        }
    }    
    _unit_system = AU;
    _matrix->set_unit_system(_unit_system);
}

void Hamiltonian::convert_to_si(){
    if (_unit_system == SI){
        return;
    }
    for (int idx_r = 0; idx_r < _numpoints; idx_r++){    
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            cdouble value = _matrix->get(idx_r, iorb);
            _matrix->set(value*energy_au2ev, idx_r, iorb);
        }
    }    
    _unit_system = SI;
    _matrix->set_unit_system(_unit_system);
}

void Hamiltonian::diagonalization_W_KSPACE(MatrixField *Umatrix){
    if (_space_type != KSPACE || _gauge != WGAUGE){
        std::cout<<"Error calculating diagonal transformation. Space type not KSPACE or gauge not W. "<<std::endl;
        std::exit(-1);
    }
    double *eigenvalues = new double[_num_orbitals];
    cdouble *matrix_temp = new cdouble[_num_orbitals*_num_orbitals];
    for(int idx_k = 0; idx_k<_numpoints; idx_k ++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb=0; jorb<_num_orbitals; jorb++){
                matrix_temp[jorb*_num_orbitals + iorb] = _matrix->get(idx_k, iorb*_num_orbitals +jorb);
            }
        }
        int info = LAPACKE_zheevd(
            LAPACK_COL_MAJOR,
            'V',                 // compute eigenvectors
            'L',                 // use upper triangle
            _num_orbitals,
            reinterpret_cast<lapack_complex_double*>(matrix_temp),
            _num_orbitals,
            eigenvalues
        ); 
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                Umatrix->set(matrix_temp[jorb*_num_orbitals + iorb], idx_k, iorb*_num_orbitals + jorb);
            }
        }
        //std::cout<<std::endl;
    }
    delete[] eigenvalues;
    delete[] matrix_temp;
}


