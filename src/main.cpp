#include <iostream>
#include "vec3_util.h"
#include "constants.h"
#include "swe.h"
#include "settings.h"

int main(){
    SWESim *swe;
    Settings* settings;
    
    settings = new Settings();
    std::string path = "../hmcase0_tb.dat";

    swe = new SWESim();
    swe->set_path_tb(path);
    swe->set_settings(settings);
    swe->init();
    //(*swe).test_files();
    (*swe).run_simulation();
    delete swe;
    delete settings;
    return 0;
}
