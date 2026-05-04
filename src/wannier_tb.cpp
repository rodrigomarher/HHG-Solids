#include <iostream>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include "constants.h"
#include "wannier_tb.h"


WannierTB::WannierTB(const std::string& path_tb){
    _path_tb = path_tb;
    parse_tb();
}

void WannierTB::parse_tb(){
    std::ifstream file(_path_tb);
    if(!file){
        throw std::runtime_error("Error opening wannier tb file.");
    }

    std::string line;
    std::getline(file, line);
    comment = line;

    for(int i=0; i<3; i++){
        std::getline(file, line);
        std::istringstream iss(line);
        iss>> lattice_vector[i][0] >> lattice_vector[i][1] >> lattice_vector[i][2];
    }

    std::getline(file, line);
    num_orb = std::stoi(line);
    std::getline(file, line);
    num_sites = std::stoi(line);

    degeneracy = new int[num_sites];
    int num_degeneracy_lines = (num_sites-1)/15+1;
    int cnt = 0;
    for (int i =0; i<num_degeneracy_lines; i++){
        std::getline(file, line);
        std::istringstream iss;
        iss.str(line);
        for(std::string tmp; std::getline(iss, tmp, ' ');){
            if(tmp != ""){
                degeneracy[cnt] = std::stoi(tmp);
                cnt++;
            } 
        }
    }
    
    std::getline(file, line);

    R_idx = new int*[num_sites];
    for (int i=0; i<num_sites; i++){
        R_idx[i] = new int[3];
    }

    hopping = new cdouble[num_sites*num_orb*num_orb];
    int size_block = num_orb * num_orb;
    for (int i=0; i<num_sites; i++){
        std::getline(file, line);
        std::istringstream iss(line);
        iss >> R_idx[i][0] >> R_idx[i][1] >> R_idx[i][2]; 

        for (int j=0; j<size_block; j++){
            int idx_orb1, idx_orb2, idx_hopping;
            double matrix_element_real, matrix_element_imag;
            std::getline(file, line);
            std::istringstream iss(line);
            iss >> idx_orb1 >> idx_orb2 >> matrix_element_real >> matrix_element_imag;
            idx_orb1 -= 1;
            idx_orb2 -= 1;
            // In wannier _tb.dat file, the first orbital index is the fast index
            idx_hopping = i*num_orb*num_orb + idx_orb2*num_orb + idx_orb1;
            hopping[idx_hopping] = cdouble(matrix_element_real, matrix_element_imag);
                
        }
        std::getline(file, line);

    }
    //std::getline(file, line);

    x_bc = new cdouble[num_sites*num_orb*num_orb];
    y_bc = new cdouble[num_sites*num_orb*num_orb];
    z_bc = new cdouble[num_sites*num_orb*num_orb];
    for (int i=0; i<num_sites; i++){
        std::getline(file, line);
        for (int j=0; j<size_block; j++){
            int idx_orb1, idx_orb2, idx_bc;
            double x_bc_real, x_bc_imag;
            double y_bc_real, y_bc_imag;
            double z_bc_real, z_bc_imag;

            std::getline(file, line);
            std::istringstream iss(line);
            iss >> idx_orb1 >> idx_orb2 >> x_bc_real >> x_bc_imag 
                                        >> y_bc_real >> y_bc_imag
                                        >> z_bc_real >> z_bc_imag;   
            idx_orb1 -= 1;
            idx_orb2 -= 1;
            // In wannier _tb.dat file, the first orbital index is the fast index
            idx_bc = i*num_orb*num_orb + idx_orb2*num_orb + idx_orb1;
            x_bc[idx_bc] = cdouble(x_bc_real, x_bc_imag);
            y_bc[idx_bc] = cdouble(y_bc_real, y_bc_imag);
            z_bc[idx_bc] = cdouble(z_bc_real, z_bc_imag);
                
        }
        std::getline(file, line);

    }
    std::getline(file, line);

}

void WannierTB::print_lattice(){
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"Lattice vectors:"<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
    for (int i=0; i<3;i++){
        std::cout<<"{"<< lattice_vector[i][0]<<", "<< lattice_vector[i][1] << ", "<< lattice_vector[i][2]<<"}"<< std::endl; 
    }


}


void WannierTB::print_hamiltonian(){
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"Hamiltonian matrix elements:"<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
    for(int i=0; i<num_sites;i++){
        std::cout<<"{"<<R_idx[i][0]<<", "<<R_idx[i][1]<<", "<<R_idx[i][2]<<"}"<<std::endl;
        for(int j=0; j<num_orb*num_orb; j++){
            std::cout<<"\t"<<hopping[i*num_orb*num_orb + j]<<std::endl;
        }
        std::cout<<std::endl;
    }
}

void WannierTB::print_berryconnection(){
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"Berry connection matrix elements:"<<std::endl;
    std::cout<<"----------------------------------------"<<std::endl;
    for(int i=0; i<num_sites;i++){
        std::cout<<"{"<<R_idx[i][0]<<", "<<R_idx[i][1]<<", "<<R_idx[i][2]<<"}"<<std::endl;
        for(int j=0; j<num_orb*num_orb; j++){
            std::cout << "\t" << x_bc[i*num_orb*num_orb + j] << " "
                              << y_bc[i*num_orb*num_orb + j] << " " 
                              << z_bc[i*num_orb*num_orb + j] << std::endl;
        }
        std::cout<<std::endl;
    }
}

void WannierTB::print_tb(){
    std::string unitsystem = _unit_system == SI ? "International system" : "Atomic units";
    std::cout<<"File path: "<<_path_tb<<std::endl;
    std::cout<<"Unit system: " << unitsystem <<std::endl;
    std::cout<<"Comment: "<<comment<<std::endl;
    std::cout<<std::endl;
    std::cout<<"Number of orbitals: " << num_orb << std::endl;
    std::cout<<"Number of sites: " << num_sites << std::endl;
    std::cout<<std::endl;
    std::cout<<"Site degeneracy: "<<std::endl; 
    for (int i=0; i<num_sites; i++){
        std::cout<<degeneracy[i]<<" ";
    }

    std::cout<<std::endl;
    print_lattice();
    std::cout<<std::endl;
    print_hamiltonian();
    std::cout<<std::endl;
    print_berryconnection(); 
    std::cout<<std::endl;

}

void WannierTB::convert_to_au(){
    if(_unit_system == AU){
        return;
    }
    for(int i=0; i<num_sites*num_orb*num_orb; i++){
        hopping[i] *= energy_ev2au;
        x_bc[i] *= length_A2au;
        y_bc[i] *= length_A2au;
        z_bc[i] *= length_A2au;
    }
    for(int i=0; i<3; i++){
        lattice_vector[i][0] *= length_A2au;
        lattice_vector[i][1] *= length_A2au;
        lattice_vector[i][2] *= length_A2au;
    }
    _unit_system = AU;
}

void WannierTB::convert_to_si(){
    if(_unit_system == SI){
        return;
    }
    for(int i=0; i<num_sites*num_orb*num_orb;i++){
        hopping[i] *= energy_au2ev;
        x_bc[i] *= length_au2A;
        y_bc[i] *= length_au2A;
        z_bc[i] *= length_au2A;
    }
    for(int i=0; i<3; i++){
        lattice_vector[i][0] *= length_au2A;
        lattice_vector[i][1] *= length_au2A;
        lattice_vector[i][2] *= length_au2A;
    }
    _unit_system = SI;
}

int WannierTB::unit_system(){
    return _unit_system;
}

int WannierTB::check_vec_in_wannier(int *vec){
    bool flag = false;
    for(int i=0; i<num_sites; i++){
        flag = (vec[0] == R_idx[i][0]) && (vec[1] == R_idx[i][1]) && (vec[2] == R_idx[i][2]);
        if(flag == true){
            return i;
        }
    }
    return -1;
}

WannierTB::~WannierTB(){
    delete[] degeneracy;
    for(int i=0; i<num_sites; i++){
        delete[] R_idx[i];
    }
    delete[] R_idx;
}

