#pragma once

#include <complex>
#include "grid.h"
#include "matrix_field.h" 
#include "wannier_tb.h"

#define WGAUGE 0
#define HGAUGE 1
#define PGAUGE 2
#define LGAUGE 3
#define RSPACE 0
#define KSPACE 1

class Operator{
    protected:
        int _gauge;
        int _unit_system;
        int _space_type;
        int _numpoints;
        int _num_orbitals;
        MatrixField* _matrix;
        Settings* _settings;
        Grid* _grid;
        WannierTB* _wannier;
        virtual void _setup() = 0;
        
    public:
        Operator(Settings* settings, Grid *grid, WannierTB *wannier, int gauge, int space_type); 
        virtual void convert_to_au() = 0;
        virtual void convert_to_si() = 0;
        int unit_system();
        int gauge();
        int space_type();
        cdouble** data_ptr();
        void convert_to_k();
        void convert_to_r();
        void convert_w_to_h(MatrixField* U);
        void convert_h_to_w(MatrixField* U);
        void convert_l_to_p(double vecpot_x, double vecpot_y, double vecpot_z);
        void convert_p_to_l(double vecpot_x, double vecpot_y, double vecpot_z);
        void write(const std::string &filename);
        cdouble get(const int idx, const int iorb);      
        void set(cdouble value, const int idx, const int iorb);      
        void copy_data(Operator* op_dest);
        ~Operator();
};
