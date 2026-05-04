#pragma once

#include "grid.h"
#include "operator.h"
#include "hamiltonian.h"
#include "berry_connection.h"
#include "wannier_tb.h"

#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

class Velocity : public Operator{
    private: 
       void _setup() override; 
       int _axis;
    public:
        Velocity(Settings *settings, Grid *grid, WannierTB *wannier, int gauge, int space_type)
            : Operator(settings, grid, wannier , gauge, space_type){_setup();};
         
        void convert_to_au() override;
        void convert_to_si() override;
        void setup(Hamiltonian* ham, BerryConnection* rbc, const int axis);
        void diagonalization_W_KSPACE(MatrixField *matrix);
};
