#pragma once

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include <cuComplex.h>
#include <cufft.h>
#include "settings.h"
#include "grid.h"
#include "hamiltonian.h"
#include "rdm.h"
#include "efield.h"
#include "berry_connection.h"
#include "wannier_tb.h"
#include "rdm_cuda.h"

#define cdouble_cuda cuDoubleComplex


class Solver_cuda{
    private:
        Settings* _settings;
        Grid* _grid;
        Hamiltonian* _hamiltonian;
        //RDM* _rho;
	RDM_cuda* _rho;
        Efield* _efield;
        BerryConnection** _r_bc;
        WannierTB* _wannier;

	cudaStream_t* _stream;
	cufftHandle _cufft_plan;

        int _num_points;
        int _num_orbitals;
        double _dt;

	cdouble_cuda* _d_rho;
	cdouble_cuda* _d_rho_k;
	cdouble_cuda* _d_h0;
	cdouble_cuda* _d_xbc;
	cdouble_cuda* _d_ybc;
	cdouble_cuda* _d_zbc;
        cdouble_cuda* _d_heff;
        cdouble_cuda* _d_heff_k;
        cdouble_cuda* _d_k1;
        cdouble_cuda* _d_k2;
        cdouble_cuda* _d_k3;
        cdouble_cuda* _d_k4;
	cdouble_cuda* _d_comm_k;
	cdouble_cuda* _d_peierls_phase;
	double* _d_r_vec_x;
	double* _d_r_vec_y;
	double* _d_r_vec_z;
	double *_d_Ex;
	double *_d_Ey;
	double *_d_Ez;
	double *_d_Ax;
	double *_d_Ay;
	double *_d_Az;
        cdouble* _peierls_phase;
        void _clear_kn();

        void _update_k1_conv_cuda(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_k2_conv_cuda(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_k3_conv_cuda(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_k4_conv_cuda(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_heff(const double ex, const double ey, const double ez,
                          const double ax, const double ay, const double az);
        
        void _allocate();
        void _deallocate();
	void _init_device_arrays();
	void _create_cufft_plan();
	void _calculate_peierls_phase(const double ax, const double ay, const double az);
	void _step_rho();
	void _copy_rho_device_to_host();
    public:
        Solver_cuda();
        Solver_cuda(Settings* settings,
               Grid* grid, 
               Hamiltonian* hamiltonian,
               RDM_cuda* rdm,
               BerryConnection** r_bc,
               Efield* efield,
               WannierTB* wannier,
               cudaStream_t* stream);
        
        void step_rk4(const int ti);
        void init();
        ~Solver_cuda();
};
#endif
