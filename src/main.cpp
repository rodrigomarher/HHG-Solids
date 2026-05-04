#include <iostream>
#include "vec3_util.h"
#include "constants.h"
#include "swe.h"

int main(){
    SWESim *swe;
    swe = new SWESim("../hmcase0_tb.dat");
    //(*swe).test_files();
    (*swe).run_simulation();
    delete swe;
    return 0;
}
