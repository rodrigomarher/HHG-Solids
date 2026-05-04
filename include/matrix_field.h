#pragma once
#include <complex>
#include <string>
#include "grid.h"
#include "settings.h"
#include "wannier_tb.h"

#define cdouble std::complex<double>
class MatrixField{
    protected:
        int _unit_system;
        Settings *_settings;
        WannierTB *_wannier;

        cdouble **_data;
        int _numpoints, _num_orbitals;
        void _alloc_matrix();
    public:
        MatrixField(Settings *settings, WannierTB *wannier); 
        void write_to_file(const std::string &filename);
        void set_unit_system(const int value);
        int get_unit_system();
        void set(cdouble value, const int idx, const int iorb);
        cdouble get(const int idx, const int iorb);
        cdouble& operator() (const int idx, const int iorb);
        cdouble** data_ptr();
        ~MatrixField();
};
