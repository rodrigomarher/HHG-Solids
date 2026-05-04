#pragma once

#include <string>
#include <complex>
#include "constants.h"

#define cdouble std::complex<double>

class WannierTB{
    private:
        std::string _path_tb;
        int _unit_system = SI;

    public:
        std::string comment;
        double lattice_vector[3][3];
        int num_orb;
        int num_sites;
        int *degeneracy;
        int **R_idx;
        cdouble *hopping;
        cdouble *x_bc;
        cdouble *y_bc;
        cdouble *z_bc;

        WannierTB(const std::string &filename);
        void parse_tb();
        void print_lattice();
        void print_hamiltonian();
        void print_berryconnection();
        void print_tb();
        void convert_to_au();
        void convert_to_si();
        int unit_system();
        int check_vec_in_wannier(int *vec);
        ~WannierTB();
};
