#include <iostream>
#include "settings_swe.h"

Settings_swe::Settings_swe(){}

void Settings_swe::print_settings_swe(){
    std::cout<<"Grid parameters: "<<std::endl;
    std::cout<<"   nr1: "<<nr1<< ", nr2: "<< nr2<<", nr3: "<<nr3<< std::endl;
    std::cout<<std::endl;
    std::cout<<"Number of orbitals: "<<num_orb<<std::endl;

}
