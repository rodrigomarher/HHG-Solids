#include <string>
#include "swe.h"
#include "settings.h"

extern "C"{
//    SWESim* SWESim_new(char* path_tb, Settings* settings){std::string str(path_tb); return new SWESim(str, settings);}
    SWESim* SWESim_new(){return new SWESim();}
    void SWESim_run_simulation(SWESim* swe){swe->run_simulation();}
    void SWESim_test_files(SWESim* swe){swe->test_files();}
    void SWESim_set_path_tb(SWESim* swe, char* path_tb){std::string str(path_tb); swe->set_path_tb(str);}
    void SWESim_set_settings(SWESim* swe, Settings* settings){swe->set_settings(settings);}
    void SWESim_init(SWESim* swe){swe->init();}
    void SWESim_restart(SWESim* swe){swe->restart();}
    void SWESim_get_current(SWESim* swe, double* t , cdouble* jx, cdouble* jy, cdouble* jz){swe->get_current(t, jx, jy, jz);}
    void SWESim_delete(SWESim* swe){delete swe;}

    Settings* Settings_new(){return new Settings();}
    void Settings_set_nr1(Settings* s, int val){s->nr1 = val;}
    void Settings_set_nr2(Settings* s, int val){s->nr2 = val;}
    void Settings_set_nr3(Settings* s, int val){s->nr3 = val;}
    void Settings_set_tmax(Settings* s, double val){s->tmax = val;} 
    void Settings_set_dt(Settings* s, double val){s->dt = val;}
    void Settings_set_intensity(Settings* s, double val){s->intensity_wcm2 = val;}
    void Settings_set_lambda(Settings* s, double val){s->lambda_nm = val;}
    void Settings_set_tmax_field(Settings* s, double val){s->tmax_field = val;}
    void Settings_set_pol_vec(Settings* s, double* val){s->sx = val[0]; s->sy = val[1]; s->sz = val[2];}
    void Settings_set_phi_vec(Settings* s, double* val){s->phi_x = val[0]; s->phi_y = val[1]; s->phi_z = val[2];}
}
