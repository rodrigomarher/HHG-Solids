#pragma once

#include "grid.h"
#include "operator.h"
#include "wannier_tb.h"


class Hamiltonian : public Operator{
    private: 
       void _setup() override; 
    public:
        Hamiltonian(Settings *settings, Grid *grid, WannierTB *wannier, int gauge, int space_type)
            : Operator(settings, grid, wannier , gauge, space_type){_setup();};
         
        void convert_to_au() override;
        void convert_to_si() override;
        void diagonalization_W_KSPACE(MatrixField *matrix);
};
