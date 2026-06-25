#pragma once

#include <string>
#include "grid.h"
#include "settings_swe.h"
#include "wannier_tb.h"

#define cdouble std::complex<double>
class MatrixField{
    protected:
        int _unit_system;
        Settings_swe *_settings_swe;
        WannierTB *_wannier;

        cdouble **_data;
        int _numpoints, _num_orbitals;
        void _alloc_matrix();
    public:
        MatrixField(Settings_swe *settings_swe, WannierTB *wannier); 
        void write_to_file(const std::string &filename);
        void set_unit_system(const int value);
        int get_unit_system();
        void set(cdouble value, const int idx, const int iorb);
        cdouble get(const int idx, const int iorb);
        cdouble& operator() (const int idx, const int iorb);
        cdouble** data_ptr();
        ~MatrixField();
};

