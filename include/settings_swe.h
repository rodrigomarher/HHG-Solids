#pragma once

#include <string>
#define nr 10
class Settings_swe{
    public:
        Settings_swe();
        int nr1 = nr;
        int nr2 = nr;
        int nr3 = 1;
        
        double tmax = 90.0;
        int nt = 4096;
        double dt = tmax/(double)dt;

        double lambda_nm = 3000;
        double tmax_field = 80.0;
        double intensity_wcm2 = 1e12;
        double sx = 0.0;
        double sy = 1.0;
        double sz = 0.0;
        double phi_x = 0.0;
        double phi_y = 0.0;
        double phi_z = 0.0;

        int num_orb;
        int num_sites;
        std::string path_tb;
        std::string path_results = "results/";
        void print_settings_swe();
};
