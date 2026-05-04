#pragma once

#include <string>
#define nr 400
class Settings{
    public:
        Settings();
        int nr1 = nr;
        int nr2 = nr;
        int nr3 = 1;
        
        double tmax = 90.0;
        double dt = 10.98e-3/0.5;
        int nt = tmax/dt;

        double lambda_nm = 3000;
        double tmax_field = 80.0;
        double intensity_wcm2 = 5e10;
        double sx = 0.96;
        double sy = 0.25;
        double sz = 0.0;
        double phi_x = 0.0;
        double phi_y = 0.0;
        double phi_z = 0.0;

        int num_orb;
        int num_sites;
        std::string path_tb;
        std::string path_results = "results/";
        void print_settings();
};
