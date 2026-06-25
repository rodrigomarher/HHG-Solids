#pragma once

#include <complex>
#include <vector>
#include <string>
#include "settings_swe.h"
#include "grid.h"
#include "operator.h"
#include "rdm.h"

#define cdouble std::complex <double>

class Observable{
    private:
        int _unit_system;
        std::vector<cdouble>* _data;
        int _num_points;
        int _num_orbitals;
        
        Settings_swe* _settings_swe;
        Grid* _grid;
        RDM* _rho;
        Operator* _operator;
    public:
        Observable(Settings_swe *settings_swe, Grid* grid, RDM* rho, Operator* op);
        void calculate();
        void write(std::string filename);
        cdouble* get_ptr();
        ~Observable();
};
