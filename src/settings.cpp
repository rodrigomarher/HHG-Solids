#include <iostream>
#include "settings.h"

Settings::Settings(){}

void Settings::print_settings(){
    std::cout<<"Grid parameters: "<<std::endl;
    std::cout<<"   nr1: "<<nr1<< ", nr2: "<< nr2<<", nr3: "<<nr3<< std::endl;
    std::cout<<std::endl;
    std::cout<<"Number of orbitals: "<<num_orb<<std::endl;

}
