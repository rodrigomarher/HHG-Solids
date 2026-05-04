#pragma once
#include <fftw3.h>
#include <complex>

#define cdouble std::complex<double>
void fft3(cdouble *data, 
          fftw_complex *in, fftw_complex *out,
          const int N, fftw_plan& plan);

void ifft3(cdouble* data,
           fftw_complex* in, fftw_complex* out,
           const int N, fftw_plan& plan);

void fftshift(cdouble *data, const int N1, const int N2);
void ifftshift(cdouble *data, const int N1, const int N2);
