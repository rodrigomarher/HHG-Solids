#pragma once

#include "settings.h"
#include "grid.h"
#include "hamiltonian.h"
#include "rdm.h"
#include "efield.h"
#include "berry_connection.h"
#include "wannier_tb.h"

class Solver{
    private:
        Settings* _settings;
        Grid* _grid;
        Hamiltonian* _hamiltonian;
        RDM* _rho;
        RDM* _rho_aux;
        Efield* _efield;
        BerryConnection** _r_bc;
        WannierTB* _wannier;

        int _num_points;
        int _num_orbitals;
        double _dt;

        cdouble **_heff;
        cdouble **_k1;
        cdouble **_k2;
        cdouble **_k3;
        cdouble **_k4;        
        cdouble* _peierls_phase;
        void _clear_kn();

        void _update_k1(const double ex, const double ey, const double ez, 
                        const double ax, const double ay, const double az);
        void _update_k2(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az);
        void _update_k3(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az);
        void _update_k4(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az);
        void _update_k1_no_blas(const double ex, const double ey, const double ez, 
                                const double ax, const double ay, const double az);
        void _update_k2_no_blas(const double ex, const double ey, const double ez, 
                                const double ax, const double ay, const double az);
        void _update_k3_no_blas(const double ex, const double ey, const double ez, 
                                const double ax, const double ay, const double az);
        void _update_k4_no_blas(const double ex, const double ey, const double ez, 
                                const double ax, const double ay, const double az);
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
    public:
        Solver();
        Solver(Settings* settings,
               Grid* grid, 
               Hamiltonian* hamiltonian,
               RDM* rdm,
               BerryConnection** r_bc,
               Efield* efield,
               WannierTB* wannier);
        
        void step_rk4(const int ti, cdouble* peierls_phase);
        void init();
        ~Solver();
};
