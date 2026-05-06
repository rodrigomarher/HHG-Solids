#pragma once
#include <string>
#include "settings.h"
#include "grid.h"
#include "wannier_tb.h"
#include "operator.h"
#include "hamiltonian.h"
#include "rdm.h"
#include "berry_connection.h"
#include "velocity.h"
#include "matrix_field.h"
#include "efield.h"
#include "solver.h"
#include "solver_kspace.h"
#include "observable.h"

class SWESim{
    private:
        int _unit_system;
        Settings *_settings;
        Grid *_grid;
        WannierTB *_wannier;
        Hamiltonian *_hamiltonian;
        BerryConnection *_r_bc[3];
        Velocity *_v[3];
        RDM *_rho;
        MatrixField *_diagonalization;
        Efield* _efield;
        Solver* _solver;
        Observable* _jx;
        Observable* _jy;
        Observable* _jz;

        std::string _path_tb;
        
        void _init();
        void _convert_to_au();
        void _convert_to_si();
        void _calc_peierls_phase(double ax, double ay, double az, cdouble* peierls_phase);
    public:
        SWESim();
        SWESim(const std::string &path_tb, Settings* settings);
        void run_simulation();
        void test_files(); 
        void set_path_tb(const std::string &path_tb);
        void set_settings(Settings* settings);
        void init();
        void restart(); 
        void get_current(double* t, cdouble* jx, cdouble* jy, cdouble* jz);
        void save_current(const std::string &path);
        ~SWESim();
};
