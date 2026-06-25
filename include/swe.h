#pragma once
#include <string>
#include "settings_swe.h"
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
#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include "solver_cuda.h"
#include "rdm_cuda.h"
#include "observable_cuda.h"
#endif
#include "observable.h"

class SWESim{
    private:
        int _unit_system;
        Settings_swe *_settings_swe;
        Grid *_grid;
        WannierTB *_wannier;
        Hamiltonian *_hamiltonian;
        BerryConnection *_r_bc[3];
        Velocity *_v[3];
        RDM *_rho;
        MatrixField *_diagonalization;
        Efield* _efield;
        Observable* _jx;
        Observable* _jy;
        Observable* _jz;
        #ifdef HAVE_CUDA
        cudaStream_t* _stream;
        int _device;
        Solver_cuda* _solver;
        RDM_cuda* _rho_cuda;
        Observable_cuda* _jx_cuda;
        Observable_cuda* _jy_cuda;
        Observable_cuda* _jz_cuda;
        #else
        Solver* _solver;
        #endif

        std::string _path_tb;
        
        void _init();
        void _convert_to_au();
        void _convert_to_si();
        void _calc_peierls_phase(double ax, double ay, double az, cdouble* peierls_phase);
        void _run_simulation_cpu();
        #ifdef HAVE_CUDA
        void _run_simulation_cuda();
        #endif
    public:
        SWESim();
        SWESim(const std::string &path_tb, Settings_swe* settings_swe);
        void update_field(double *ex, double *ey, double *ez);
        void run_simulation();
        void test_files(); 
        void set_path_tb(const std::string &path_tb);
        void set_settings_swe(Settings_swe* settings_swe);
        void init();
        void restart(); 
        void get_current(double* t, cdouble* jx, cdouble* jy, cdouble* jz);
        void save_current(const std::string &path);
        ~SWESim();
};
