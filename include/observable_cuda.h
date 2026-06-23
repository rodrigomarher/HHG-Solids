#pragma once

#ifdef HAVE_CUDA

#include <cuda_runtime.h>
#include <cuComplex.h>
#include <vector>
#include "rdm_cuda.h"
#include "operator.h"
#include "grid.h"

#define cdouble_cuda cuDoubleComplex
#define cdouble std::complex<double>
class Observable_cuda{
    private:
        int _unit_system;
        int _num_points;
        int _num_orbitals;
        cdouble_cuda *_d_data;
        cdouble_cuda *_d_operator;
        double* _d_r_vec_x;
        double* _d_r_vec_y;
        double* _d_r_vec_z;
        cdouble* _data;
        cudaStream_t* _stream;
        
        Settings* _settings;
        Grid* _grid;
        RDM_cuda* _rho_cuda;
        Operator* _operator;
    
        void _allocate();
        void _deallocate();
    public:
        Observable_cuda(Settings* settings, Grid* grid, RDM_cuda* rho, Operator* op, cudaStream_t* stream);
        void init_device();
        void calculate(int ti, double ax, double ay, double az);
        void write(std::string filename);
        cdouble* get_ptr();
        cdouble_cuda* get_d_ptr();
        void copy_device_to_host();
        ~Observable_cuda();
};

#endif
