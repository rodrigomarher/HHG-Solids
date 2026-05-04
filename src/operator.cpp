#include <algorithm>
#include <iostream>
#include "cblas.h"
#include "operator.h"
#include "matrix_field.h"
#include "vec3_util.h"
#include "fftw_helper.h"

Operator::Operator(Settings* settings,Grid* grid,  WannierTB *wannier, int gauge, int space_type){
    _gauge = gauge; _unit_system = wannier->unit_system(); _space_type = space_type;
    _settings = settings; _wannier = wannier; _grid = grid; 
    _numpoints = settings->nr1*settings->nr2*_settings->nr3;
    _num_orbitals = settings->num_orb;
    
    _matrix = new MatrixField(settings, wannier);
}

int Operator::unit_system(){
    return _unit_system;
}

int Operator::gauge(){
    return _gauge;
}

int Operator::space_type(){
    return _space_type;
}

cdouble Operator::get(const int idx, const int iorb){
    return _matrix->get(idx, iorb);
}
 
void Operator::set(cdouble value, const int idx, const int iorb){
    _matrix->set(value, idx, iorb);
}

cdouble** Operator::data_ptr(){
    return  _matrix->data_ptr();
}
void Operator::convert_to_k(){
    if (_space_type == KSPACE){
        return;
    }
    cdouble* tmp1 = new cdouble[_numpoints];
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_numpoints);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_numpoints);
    fftw_plan forward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
        for(int idx_r = 0; idx_r <_numpoints; idx_r++){
            tmp1[idx_r] = _matrix->get(idx_r,iorb);           
        }
        ifftshift(tmp1, _settings->nr1, _settings->nr2);
        fft3(tmp1, in, out, _numpoints, forward);
        for(int idx_k = 0; idx_k< _numpoints; idx_k++){
            _matrix->set(tmp1[idx_k], idx_k, iorb);
        }
    }
    _space_type = KSPACE;
    delete[] tmp1;
    fftw_free(in);
    fftw_free(out);
}

void Operator::convert_to_r(){
    if (_space_type == RSPACE){
        return;
    }
    cdouble* tmp1 = new cdouble[_numpoints];
    fftw_complex* in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_numpoints);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*_numpoints);
    fftw_plan backward = fftw_plan_dft_2d(_settings->nr1, _settings->nr2, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    for(int iorb = 0; iorb<_num_orbitals*_num_orbitals; iorb++){
        for(int idx_k = 0; idx_k <_numpoints; idx_k++){
            tmp1[idx_k] = _matrix->get(idx_k,iorb);           
        }
        ifft3(tmp1, in, out, _numpoints, backward);
        ifftshift(tmp1, _settings->nr1, _settings->nr2);
        for(int idx_r = 0; idx_r< _numpoints; idx_r++){
            _matrix->set(tmp1[idx_r], idx_r, iorb);
        }
    }
    _space_type = RSPACE;
    delete[] tmp1;
    fftw_free(in);
    fftw_free(out);
}

//void Operator::convert_to_k(){
//    if(_space_type == KSPACE){
//        return;
//    }
//    cdouble **tmp_data;
//    tmp_data = new cdouble*[_numpoints];
//    for(int i=0; i<_numpoints; i++){
//        tmp_data[i] = new cdouble[_num_orbitals*_num_orbitals];
//    }
//    cdouble *phases;
//    phases = new cdouble[_numpoints];
//    for (int idx_k = 0; idx_k<_numpoints; idx_k++){
//        for (int idx_r = 0; idx_r<_numpoints; idx_r++){
//            phases[idx_r] = std::exp(cdouble(0,dot(_grid->Kvecs(idx_k), _grid->Rvecs(idx_r))));
//        }
//        for (int idx_orb = 0; idx_orb<_num_orbitals*_num_orbitals; idx_orb++){
//            cdouble sum = cdouble(0.0, 0.0);
//            for (int idx_r=0; idx_r < _numpoints; idx_r++){
//                sum += _matrix->get(idx_r, idx_orb)*phases[idx_r];
//            }
//            tmp_data[idx_k][idx_orb] = sum;
//        }
//    }
//
//    for (int idx_k=0; idx_k<_numpoints; idx_k++){
//        for (int idx_orb=0; idx_orb<_num_orbitals*_num_orbitals; idx_orb++){
//            _matrix->set( tmp_data[idx_k][idx_orb], idx_k, idx_orb);
//        }
//        delete[] tmp_data[idx_k];
//    }
//    delete[] tmp_data;
//    delete[] phases;
//    _space_type = KSPACE;
//    
//}

//void Operator::convert_to_r(){
//    if(_space_type == RSPACE){
//        return;
//    }
//    cdouble **tmp_data;
//    tmp_data = new cdouble*[_numpoints];
//    for(int i=0; i<_numpoints; i++){
//        tmp_data[i] = new cdouble[_num_orbitals*_num_orbitals];
//    }
//
//    cdouble *phases;
//    phases = new cdouble[_numpoints];
//    for (int idx_r = 0; idx_r<_numpoints; idx_r++){
//        for (int idx_k = 0; idx_k<_numpoints; idx_k++){
//            phases[idx_k] = std::exp(cdouble(0,-1.0*dot(_grid->Kvecs(idx_k), _grid->Rvecs(idx_r))));
//        }
//        for (int idx_orb = 0; idx_orb<_num_orbitals*_num_orbitals; idx_orb++){
//            cdouble sum = cdouble(0.0, 0.0);
//            for (int idx_k=0; idx_k < _numpoints; idx_k++){
//                sum += _matrix->get(idx_k, idx_orb)*phases[idx_k];
//            }
//            tmp_data[idx_r][idx_orb] = sum/(double)_numpoints;
//        }
//    }
//
//    for (int idx_r=0; idx_r<_numpoints; idx_r++){
//        for (int idx_orb=0; idx_orb<_num_orbitals*_num_orbitals; idx_orb++){
//            _matrix->set(tmp_data[idx_r][idx_orb], idx_r, idx_orb);
//        }
//        delete[] tmp_data[idx_r];
//    }
//    delete[] tmp_data;
//    delete[] phases;
//    _space_type = RSPACE;
//}

void Operator::write(const std::string& filename){
    _matrix->write_to_file(filename);
}

void Operator::convert_w_to_h(MatrixField *U){
    if (_space_type != KSPACE || _gauge != WGAUGE){
        std::cout<<"Error calculating diagonal transformation. Space type not KSPACE or gauge not W. "<<std::endl;
        std::exit(-1);
    }
    cdouble *matrix_temp_h = new cdouble [_num_orbitals*_num_orbitals]; 
    cdouble *matrix_temp_u = new cdouble [_num_orbitals*_num_orbitals]; 
    cdouble *matrix_temp_aux = new cdouble[_num_orbitals*_num_orbitals]; 
    cdouble alpha = {1.0,0.0};
    cdouble beta = {0.0,0.0};
    for(int idx_k = 0; idx_k<_numpoints; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb=0; jorb<_num_orbitals; jorb++){
                matrix_temp_h[iorb*_num_orbitals + jorb] = _matrix->get(idx_k, iorb*_num_orbitals +jorb);
                matrix_temp_u[iorb*_num_orbitals + jorb] = U->get(idx_k, iorb*_num_orbitals +jorb);
                matrix_temp_aux[iorb*_num_orbitals + jorb] = cdouble(0.0,0.0);
            }
        }
        cblas_zgemm(CblasRowMajor,
                                 CblasConjTrans,
                                 CblasNoTrans,
                                  _num_orbitals,
                                  _num_orbitals,
                                  _num_orbitals,
                                  &alpha,
                                  matrix_temp_u,
                                  _num_orbitals,
                                  matrix_temp_h,
                                  _num_orbitals,
                                  &beta,
                                  matrix_temp_aux,
                                  _num_orbitals);
        cblas_zgemm(CblasRowMajor,
                                 CblasNoTrans,
                                 CblasNoTrans,
                                  _num_orbitals,
                                  _num_orbitals,
                                  _num_orbitals,
                                  &alpha,
                                  matrix_temp_aux,
                                  _num_orbitals,
                                  matrix_temp_u,
                                  _num_orbitals,
                                  &beta,
                                  matrix_temp_h,
                                  _num_orbitals);
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                _matrix->set(matrix_temp_h[iorb*_num_orbitals + jorb], idx_k, iorb*_num_orbitals + jorb);
            }
        }
   } 
   delete[] matrix_temp_h;
   delete[] matrix_temp_u;
   delete[] matrix_temp_aux;
   _gauge = HGAUGE;
}

void Operator::convert_h_to_w(MatrixField *U){
    if (_space_type != KSPACE || _gauge != HGAUGE){
        std::cout<<"Error calculating diagonal transformation. Space type not KSPACE or gauge not H. "<<std::endl;
        std::exit(-1);
    }
    cdouble *matrix_temp_h = new cdouble [_num_orbitals*_num_orbitals]; 
    cdouble *matrix_temp_u = new cdouble [_num_orbitals*_num_orbitals]; 
    cdouble *matrix_temp_aux = new cdouble[_num_orbitals*_num_orbitals]; 

    cdouble alpha = {1.0,0.0};
    cdouble beta = {0.0,0.0};
    for(int idx_k = 0; idx_k<_numpoints; idx_k++){
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb=0; jorb<_num_orbitals; jorb++){
                matrix_temp_h[iorb*_num_orbitals + jorb] = _matrix->get(idx_k, iorb*_num_orbitals +jorb);
                matrix_temp_u[iorb*_num_orbitals + jorb] = U->get(idx_k, iorb*_num_orbitals +jorb);
                matrix_temp_aux[iorb*_num_orbitals + jorb] = cdouble(0.0,0.0);
            }
        }
        cblas_zgemm(CblasRowMajor,
                                 CblasNoTrans,
                                 CblasNoTrans,
                                  _num_orbitals,
                                  _num_orbitals,
                                  _num_orbitals,
                                  &alpha,
                                  matrix_temp_u,
                                  _num_orbitals,
                                  matrix_temp_h,
                                  _num_orbitals,
                                  &beta,
                                  matrix_temp_aux,
                                  _num_orbitals);
        cblas_zgemm(CblasRowMajor,
                                 CblasNoTrans,
                                 CblasConjTrans,
                                  _num_orbitals,
                                  _num_orbitals,
                                  _num_orbitals,
                                  &alpha,
                                  matrix_temp_aux,
                                  _num_orbitals,
                                  matrix_temp_u,
                                  _num_orbitals,
                                  &beta,
                                  matrix_temp_h,
                                  _num_orbitals);
        for(int iorb = 0; iorb<_num_orbitals; iorb++){
            for(int jorb = 0; jorb <_num_orbitals; jorb++){
                _matrix->set(matrix_temp_h[iorb*_num_orbitals + jorb], idx_k, iorb*_num_orbitals + jorb);
            }
        }
   } 
   delete[] matrix_temp_h;
   delete[] matrix_temp_u;
   delete[] matrix_temp_aux;
   _gauge = WGAUGE;
}

void Operator::convert_l_to_p(double vecpot_x, double vecpot_y, double vecpot_z){
    if (_space_type != RSPACE || _gauge != LGAUGE){
        std::cout<<"Error calculating Peierls transformation. Space type not RSPACE or gauge not L."<<std::endl;
        std::exit(-1);
    }
    for(int idx_r=0; idx_r<_numpoints; idx_r++){
        double* r_vec;
        r_vec = _grid->Rvecs(idx_r);
        cdouble peierls_phase = exp(cdouble(0.0,vecpot_x*r_vec[0] + vecpot_y*r_vec[1] + vecpot_z*r_vec[2]));
        for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            _matrix->set(peierls_phase*_matrix->get(idx_r, iorb), idx_r, iorb);
        }
    }
    _gauge = PGAUGE;
}

void Operator::convert_p_to_l(double vecpot_x, double vecpot_y, double vecpot_z){
    if (_space_type != RSPACE || _gauge != PGAUGE){
        std::cout<<"Error calculating inverse Peierls transformation. Space type not RSPACE or gauge not P."<<std::endl;
        std::exit(-1);
    }
    for(int idx_r=0; idx_r<_numpoints; idx_r++){
        double* r_vec;
        r_vec = _grid->Rvecs(idx_r);
        cdouble peierls_phase = exp(cdouble(0.0,-(vecpot_x*r_vec[0] + vecpot_y*r_vec[1] + vecpot_z*r_vec[2])));
        for(int iorb=0; iorb<_num_orbitals*_num_orbitals; iorb++){
            _matrix->set(peierls_phase*_matrix->get(idx_r, iorb), idx_r, iorb);
        }
    }
    _gauge = LGAUGE;
}

void Operator::copy_data(Operator *op){
    cdouble** data_dest = op->data_ptr();
    cdouble** data_src = data_ptr();
    for(int idx_r=0; idx_r<_numpoints; idx_r++){
        std::copy(&data_src[idx_r][0], &data_src[idx_r][_num_orbitals*_num_orbitals-1], &data_dest[idx_r][0]);
    }
}

Operator::~Operator(){
    delete _matrix;
}
