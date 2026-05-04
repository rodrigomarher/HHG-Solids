#include <iostream>
#include "lapacke.h"
#include "velocity.h"
#include "vec3_util.h"
#include "constants.h"
#include "fftw_helper.h"


void Velocity::_setup(){
    }


void Velocity::convert_to_au(){
    if (_unit_system == AU){
        return;
    }
    for (int idx_r = 0; idx_r < _numpoints; idx_r++){    
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            cdouble value = _matrix->get(idx_r, iorb);
            _matrix->set( value*velocity_Afs2au, idx_r, iorb);
        }
    }    
    _unit_system = AU;
    _matrix->set_unit_system(_unit_system);
}

void Velocity::convert_to_si(){
    if (_unit_system == SI){
        return;
    }
    for (int idx_r = 0; idx_r < _numpoints; idx_r++){    
        for (int iorb = 0; iorb < _num_orbitals*_num_orbitals; iorb++){
            cdouble value = _matrix->get(idx_r, iorb);
            _matrix->set(value*velocity_au2Afs, idx_r, iorb);
        }
    }    
    _unit_system = SI;
    _matrix->set_unit_system(_unit_system);
}

void Velocity::setup(Hamiltonian* ham, BerryConnection* rbc, const int axis){
    _axis = axis;

    //Commutator calculation using convolution
    cdouble *h0_k = new cdouble[_num_orbitals*_num_orbitals*_numpoints];
    cdouble *rbc_k  = new cdouble[_num_orbitals*_num_orbitals*_numpoints];
    cdouble *comm_k = new cdouble[_num_orbitals*_num_orbitals*_numpoints];
    cdouble *tmp_1 = new cdouble[_numpoints];
    cdouble *tmp_2 = new cdouble[_numpoints];

    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_numpoints);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_numpoints);
    fftw_plan forward = fftw_plan_dft_3d(_settings->nr1, _settings->nr2, _settings->nr3, in, out, FFTW_FORWARD, FFTW_MEASURE);
    fftw_plan backward = fftw_plan_dft_3d(_settings->nr1, _settings->nr2, _settings->nr3, in, out, FFTW_BACKWARD, FFTW_MEASURE);

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for(int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_r = 0; idx_r < _numpoints; idx_r++){
                tmp_1[idx_r] = ham->data_ptr()[idx_r][iorb*_num_orbitals + jorb];
                tmp_2[idx_r] = rbc->data_ptr()[idx_r][iorb*_num_orbitals + jorb];
            }

            fft3(tmp_1, in, out, _numpoints, forward);
            fft3(tmp_2, in, out, _numpoints, forward);

            for(int idx_k = 0; idx_k < _numpoints; idx_k++){
                h0_k[iorb*_num_orbitals*_numpoints + jorb*_numpoints + idx_k] = tmp_1[idx_k];
                rbc_k[iorb*_num_orbitals*_numpoints + jorb*_numpoints + idx_k] = tmp_2[idx_k];
            }
        }
    }
 
    for(int idx_k=0; idx_k <_numpoints; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                cdouble comm_value = {0.0, 0.0};
                for (int korb = 0; korb<_num_orbitals; korb++){
                    cdouble ham_0 = h0_k[korb*_num_orbitals*_numpoints + jorb*_numpoints + idx_k];
                    cdouble ham_1 = h0_k[iorb*_num_orbitals*_numpoints + korb*_numpoints + idx_k];
                    
                    cdouble rbc_0 = rbc_k[iorb*_num_orbitals*_numpoints + korb*_numpoints + idx_k];
                    cdouble rbc_1 = rbc_k[korb*_num_orbitals*_numpoints + jorb*_numpoints + idx_k];
                    comm_value += rbc_0*ham_0 - rbc_1*ham_1;
                }
                comm_k[iorb*_num_orbitals*_numpoints + jorb*_numpoints + idx_k] = comm_value;
            }
        }
    } 

    for(int iorb = 0; iorb < _num_orbitals; iorb++){
        for (int jorb = 0; jorb < _num_orbitals; jorb++){
            for(int idx_k = 0; idx_k <_numpoints; idx_k++){
                tmp_1[idx_k] = comm_k[iorb*_num_orbitals*_numpoints + jorb*_numpoints + idx_k];
            }

            ifft3(tmp_1, in, out, _numpoints, backward);
            ifftshift(tmp_1, _settings->nr1, _settings->nr2); 
            for(int idx_r = 0; idx_r < _numpoints; idx_r++){
                double r0 = _grid->Rvecs(idx_r)[axis];
                _matrix->set(cdouble(0.0,1.0)*(1.0*r0*ham->data_ptr()[idx_r][iorb*_num_orbitals + jorb] - 1.0*tmp_1[idx_r]), idx_r, iorb*_num_orbitals + jorb);
            }
        }
    }

    fftw_free(in);
    fftw_free(out);
    delete[] h0_k;
    delete[] rbc_k;    
    delete[] comm_k;
    delete[] tmp_1;
    delete[] tmp_2;

}
