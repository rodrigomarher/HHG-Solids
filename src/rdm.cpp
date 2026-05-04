#include <iostream>
#include "rdm.h"
#include "constants.h"

void RDM::_setup(){
    return;
}

void RDM::convert_to_au(){
    if(_unit_system == AU){
        return;
    }
    _unit_system = AU;
}

void RDM::convert_to_si(){
    if(_unit_system == SI){
        return;
    }
    _unit_system = SI;
}

void RDM::calculate_equilibrium(Hamiltonian* ham){
    if(ham->space_type() != KSPACE || ham->gauge() !=HGAUGE 
       || _space_type != KSPACE || _gauge !=HGAUGE){
        std::cout<<"Calculation of RDM equilibirum population error. Hamiltonian and RDM must be in space type KSPACE and  gauge H"<<std::endl;
    }
    
    for(int idx_k = 0; idx_k < _numpoints; idx_k++){
        for(int iorb = 0; iorb < _num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                cdouble energy_value = ham->get(idx_k, iorb*_num_orbitals + jorb);
                cdouble value = energy_value.real() < 0.0 && iorb == jorb ? cdouble(1.0, 0.0) : cdouble(0.0, 0.0);
                _matrix->set(value, idx_k, iorb*_num_orbitals + jorb);
            }
        }
    }
}

