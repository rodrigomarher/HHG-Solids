#pragma once
#include "constants.h"
#include "settings.h"
#include "grid.h"

class Efield{
    private:
              
        int _unit_system;
        double _lambda; 
        double _omega;
        double _tmax_field;
        double *_t;
        double _svec[3];
       
        Grid* _grid; 
        Settings* _settings;

        void _allocate();
        double _env_sin2(double ti);
    public:
        double *E_x;
        double *E_y;
        double *E_z;
        double *A_x;
        double *A_y;
        double *A_z;

        Efield(Settings* settings, Grid* grid, int unit_system);
        int unit_system();
        void convert_to_si();
        void convert_to_au();
        void write();
        void init_fields();
        ~Efield(); 
}; 
