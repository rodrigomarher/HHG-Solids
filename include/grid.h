#pragma once

#include <cmath>
#include "settings.h"
#include "wannier_tb.h"
#include "constants.h"

class Grid{
    private:
        int _unit_system;
        const Settings *_settings;
        const WannierTB *_wannier;
        int *_n1, *_n2, *_n3;
        int *_m1, *_m2, *_m3;
        double *_t;
        double **_Rvecs;
        double **_Kvecs;
        double _lattice_vector[3][3];
        double _reciprocal_vector[3][3];
        void _calculate_reciprocal_vectors();
        void _allocate_grid();
        void _setup_grid();
    public:
        Grid();
        Grid(Settings *settings, WannierTB *wannier); 
        int unit_system();
        double get_r(int idx);
        double get_k(int idx);
        void print_lattice();
        void convert_to_au();
        void convert_to_si();
        void supercell_to_file(const std::string &filename);
        void reciprocal_to_file(const std::string &filename);
        double* Rvecs(const int idx); 
        double* Kvecs(const int idx); 
        double t(const int idx);
        double** Rvecs();
        double** Kvecs();
        double* t();
        void lattice_vector(double **output);
        void reciprocal_vector(double **output);
        ~Grid();
};
