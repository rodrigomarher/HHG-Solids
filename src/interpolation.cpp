#include "interpolation.h"

inline int pmod(int x, int n) {
    int r = x % n;
    return (r < 0) ? r + n : r;
}

cdouble cubicInterpolate(cdouble p0, cdouble p1, cdouble p2, cdouble p3, double t) {
    cdouble a0 = -0.5*p0 + 1.5*p1 - 1.5*p2 + 0.5*p3;
    cdouble a1 = p0 - 2.5*p1 + 2.0*p2 - 0.5*p3;
    cdouble a2 = -0.5*p0 + 0.5*p2;
    cdouble a3 = p1;
    return ((a0*t + a1)*t + a2)*t + a3;
}

// Bicubic interpolation on periodic grid
cdouble bicubicInterpolate(
    cdouble** data, int Nx, int Ny,
    double u, double v, const int idx_orb) {
    // Convert fractional coords to grid index space
    double x = u * Nx;
    double y = v * Ny;

    int ix = (int)floor(x);
    int iy = (int)floor(y);

    double tx = x - (double)ix;
    double ty = y - (double)iy;

    cdouble arr[4];

    for (int m = -1; m <= 2; m++) {
        cdouble col[4];
        for (int n = -1; n <= 2; n++) {
            int xi = pmod(ix + n, Nx);
            int yi = pmod(iy + m, Ny);
            col[n + 1] = data[yi * Nx + xi][idx_orb];
        }
        arr[m + 1] = cubicInterpolate(col[0], col[1], col[2], col[3], tx);
    }

    return cubicInterpolate(arr[0], arr[1], arr[2], arr[3], ty);
}

