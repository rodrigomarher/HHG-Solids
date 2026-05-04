#include <fftw3.h>
#include <complex>
#include "fftw_helper.h"

void fft3(cdouble *data, 
          fftw_complex *in, fftw_complex *out,
          const int N, fftw_plan& plan){
    for (int i=0; i<N; i++){
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }
    fftw_execute(plan);
    for (int i = 0; i < N; i++) {
        data[i] = cdouble(out[i][0], out[i][1]);
    } 
}

void ifft3(cdouble* data,
           fftw_complex* in, fftw_complex* out,
           const int N, fftw_plan& plan)
{

    for (int i = 0; i < N; i++) {
        in[i][0] = data[i].real();
        in[i][1] = data[i].imag();
    }

    fftw_execute(plan);

    // Normalize (FFTW does NOT normalize)
    for (int i = 0; i < N; i++) {
        data[i] = cdouble(out[i][0], out[i][1]) /(double)N;
    }
}

void fftshift(cdouble *data, const int N1, const int N2){
    cdouble *data_tmp = new cdouble[N1*N2];
    int h1 = N1/2;
    int h2 = N2/2;

    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            int ii = (i + h1) % N1;
            int jj = (j + h2) % N2;

            data_tmp[ii * N2 + jj] = data[i * N2 + j];
        }
    }

    // copy back
    for (int i = 0; i < N1 * N2; ++i) {
        data[i] = data_tmp[i];
    }

    delete[] data_tmp;
}

void ifftshift(cdouble *data, const int N1, const int N2){
    cdouble *data_tmp = new cdouble[N1*N2];
    int h1 = (N1+1)/2;
    int h2 = (N2+1)/2;

    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            int ii = (i + h1) % N1;
            int jj = (j + h2) % N2;

            data_tmp[ii * N2 + jj] = data[i * N2 + j];
        }
    }

    // copy back
    for (int i = 0; i < N1 * N2; ++i) {
        data[i] = data_tmp[i];
    }

    delete[] data_tmp;

}
