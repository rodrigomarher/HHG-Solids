#pragma once
#include <complex>

#define cdouble std::complex<double>
cdouble cubicInterpolate(cdouble p0, cdouble p1, cdouble p2, cdouble p3, double t); 

// Bicubic interpolation on periodic grid
cdouble bicubicInterpolate(
    cdouble** data, int Nx, int Ny,
    double u, double v, const int idx_orb); 
