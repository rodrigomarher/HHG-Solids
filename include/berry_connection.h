#pragma once

#include "grid.h"
#include "operator.h"
#include "hamiltonian.h"
#include "wannier_tb.h"

#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

class BerryConnection : public Operator{
    private:
        int _axis;
        void _setup() override;
    public:
        BerryConnection(Settings *settings, Grid *grid, WannierTB *wannier, int gauge, int space_type, int axis)
            : Operator(settings, grid, wannier, gauge, space_type), _axis{axis}{_setup();}
        void convert_to_au() override;
        void convert_to_si() override;
        int axis();
};
