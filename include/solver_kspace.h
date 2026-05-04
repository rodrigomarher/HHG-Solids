#pragma once

#include "settings.h"
#include "grid.h"
#include "hamiltonian.h"
#include "rdm.h"
#include "efield.h"
#include "berry_connection.h"
#include "wannier_tb.h"

class Solver_kspace{
    private:
        Settings* _settings;
        Grid* _grid;
        Hamiltonian* _hamiltonian;
        RDM* _rho;
        Efield* _efield;
        BerryConnection** _r_bc;
        WannierTB* _wannier;

        int _num_points;
        int _num_orbitals;
        double _dt;

        cdouble **_heff_k;
        cdouble **_rho_k;
        cdouble **_h0_k;
        cdouble **_xbc_k;
        cdouble **_ybc_k;
        cdouble **_zbc_k;
        cdouble **_k1;
        cdouble **_k2;
        cdouble **_k3;
        cdouble **_k4;        
        void _clear_kn();

        void _update_heff(const double ex, const double ey, const double ez, 
                        const double ax, const double ay, const double az);
        void _update_k1(const double ex, const double ey, const double ez, 
                        const double ax, const double ay, const double az);
        void _update_k2(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az);
        void _update_k3(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az);
        void _update_k4(const double ex, const double ey, const double ez,
                        const double ax, const double ay, const double az);
        void _allocate();
        void _deallocate();
        void _calc_commutator(cdouble *A, cdouble *B, cdouble *C, const int n);
        void _convert_to_kspace(cdouble** data);
        void _convert_to_rspace(cdouble** data);
        void _fill_aux_data();
    public:
        Solver_kspace();
        Solver_kspace(Settings* settings,
               Grid* grid, 
               Hamiltonian* hamiltonian,
               RDM* rdm,
               BerryConnection** r_bc,
               Efield* efield,
               WannierTB* wannier);
        void init(); 
        void step_rk4(const int ti);
        ~Solver_kspace();
};
