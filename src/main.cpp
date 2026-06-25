#include <iostream>
#include "vec3_util.h"
#include "constants.h"
#include "swe.h"
#include "settings_swe.h"

int main(){
    SWESim *swe;
    Settings_swe* settings_swe;
    
    settings_swe = new Settings_swe();
    std::string path = "../hmcase0_tb.dat";

    swe = new SWESim();
    swe->set_path_tb(path);
    swe->set_settings_swe(settings_swe);
    swe->init();
    //(*swe).test_files();
    (*swe).run_simulation();
    delete swe;
    delete settings_swe;
    return 0;
}
