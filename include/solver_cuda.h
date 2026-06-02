#pragma once

#ifdef HAVE_CUDA
#include <cuda_runtime.h>
#include <cuComplex.h>
#include "settings.h"
#include "grid.h"
#include "hamiltonian.h"
#include "rdm.h"
#include "efield.h"
#include "berry_connection.h"
#include "wannier_tb.h"

#define cdouble_cuda cuDoubleComplex
class Solver_cuda{
    private:
        Settings* _settings;
        Grid* _grid;
        Hamiltonian* _hamiltonian;
        RDM* _rho;
        Efield* _efield;
        BerryConnection** _r_bc;
        WannierTB* _wannier;
	cudaStream_t* _stream;

        int _num_points;
        int _num_orbitals;
        double _dt;

	cdouble_cuda* _d_rho;
	cdouble_cuda* _d_h0;
	cdouble_cuda* _d_xbc;
	cdouble_cuda* _d_ybc;
	cdouble_cuda* _d_zbc;
        cdouble_cuda* _d_heff;
        cdouble_cuda* _d_k1;
        cdouble_cuda* _d_k2;
        cdouble_cuda* _d_k3;
        cdouble_cuda* _d_k4;
	double *_d_Ex;
	double *_d_Ey;
	double *_d_Ez;
	double *_d_Ax;
	double *_d_Ay;
	double *_d_Az;
        cdouble* _peierls_phase;
        void _clear_kn();

        void _update_k1_conv_fftw(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_k2_conv_fftw(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_k3_conv_fftw(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_k4_conv_fftw(const double ex, const double ey, const double ez, 
                                  const double ax, const double ay, const double az);
        void _update_heff(const double ex, const double ey, const double ez,
                          const double ax, const double ay, const double az);
        
        void _allocate();
        void _deallocate();
        void _calc_commutator(cdouble *A, cdouble *B, cdouble *C, const int n);
	void _init_device_h0_rbc_rho();
    public:
        Solver_cuda();
        Solver_cuda(Settings* settings,
               Grid* grid, 
               Hamiltonian* hamiltonian,
               RDM* rdm,
               BerryConnection** r_bc,
               Efield* efield,
               WannierTB* wannier);
        
        void step_rk4(const int ti, cdouble* peierls_phase);
        void init();
        ~Solver_cuda();
};
#endif
