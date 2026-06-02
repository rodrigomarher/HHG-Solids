/**
 * prop_q.c
 *
 * C implementation of the far-field propagation integral and the underlying
 * cubic B-spline evaluator (equivalent to scipy.interpolate.splev with k=3).
 *
 * Build as a shared library via CMake (see CMakeLists.txt).
 * Call from Python via the ctypes wrapper (prop_q_ctypes.py).
 */

#include "prop_q.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

/* -------------------------------------------------------------------------
 * Constants that mirror the Python code
 * ---------------------------------------------------------------------- */
#define N_PHI        256
#define RHO_START    22.0
#define RHO_STEP      0.5
/* np.arange(22, 35, 0.5) has ceil((35-22)/0.5) = 26 elements */
#define N_RHO        26
#define PI           3.14159265358979323846
#define SPLINE_DEG   3      /* cubic */

/* -------------------------------------------------------------------------
 * bspline_eval – de Boor's algorithm for a single cubic B-spline
 *
 * Reproduces scipy.interpolate.splev(x, (t, c, 3), ext=3)
 * (ext=3 = clamp to boundary value rather than extrapolate).
 *
 * Algorithm reference:
 *   C. de Boor, "A Practical Guide to Splines", Springer 2001, Chapter X.
 * ---------------------------------------------------------------------- */
double bspline_eval(const double *t, int n_t, const double *c, double x)
{
    const int k = SPLINE_DEG;
    /* Number of usable control points = n_t - k - 1 */
    const int n = n_t - k - 1;

    /* Clamp x to the valid domain [t[k], t[n]] */
    const double x_min = t[k];
    const double x_max = t[n];
    if (x < x_min) x = x_min;
    if (x > x_max) x = x_max;

    /* ------------------------------------------------------------------ *
     * Locate the knot span:  find the largest index 'span' in [k, n-1]   *
     * such that t[span] <= x.                                              *
     * Special case: x == x_max maps to the last interior span (n-1).      *
     * ------------------------------------------------------------------ */
    int span = k;
    while (span < n - 1 && t[span + 1] <= x)
        ++span;

    /* ------------------------------------------------------------------ *
     * De Boor's triangular algorithm.                                      *
     * d[j] starts as c[span-k+j] for j = 0..k, then is updated in-place. *
     * ------------------------------------------------------------------ */
    double d[SPLINE_DEG + 1];
    for (int j = 0; j <= k; ++j)
        d[j] = c[span - k + j];

    for (int r = 1; r <= k; ++r) {
        for (int j = k; j >= r; --j) {
            const int l     = span - k + j;
            const double lo = t[l];
            const double hi = t[l + k - r + 1];
            const double denom = hi - lo;
            const double alpha = (denom > 1e-15) ? (x - lo) / denom : 0.0;
            d[j] = (1.0 - alpha) * d[j - 1] + alpha * d[j];
        }
    }

    return d[k];
}

/* -------------------------------------------------------------------------
 * prop_q_c – main far-field integral
 * ---------------------------------------------------------------------- */
void prop_q_c(
    const double *t_real, int n_t_real, const double *c_real,
    const double *t_imag, int n_t_imag, const double *c_imag,
    double q,
    const double *theta, int n_theta,
    const double *omega, int n_omega,
    double *out_real,
    double *out_imag)
{
    /* ------------------------------------------------------------------ *
     * Pre-compute phi grid: linspace(0, 2*pi, 128)                        *
     * ------------------------------------------------------------------ */
    double phi[N_PHI];
    for (int p = 0; p < N_PHI; ++p)
        phi[p] = 2.0 * PI * p / (N_PHI - 1);

    /* ------------------------------------------------------------------ *
     * Pre-compute rho grid: arange(22, 35, 0.5) → 26 values              *
     * ------------------------------------------------------------------ */
    double rho_array[N_RHO];
    for (int r = 0; r < N_RHO; ++r)
        rho_array[r] = RHO_START + r * RHO_STEP;

    /* ------------------------------------------------------------------ *
     * Pre-compute spline values for every phi point so we evaluate the    *
     * B-spline only N_PHI times instead of N_PHI * n_theta * n_omega.    *
     * ------------------------------------------------------------------ */
    double spline_real[N_PHI];
    double spline_imag[N_PHI];
    for (int p = 0; p < N_PHI; ++p) {
        spline_real[p] = bspline_eval(t_real, n_t_real, c_real, phi[p]);
        spline_imag[p] = bspline_eval(t_imag, n_t_imag, c_imag, phi[p]);
    }

    const double k_val = q * 2.0 * PI / 3.0;

    /* Zero-initialise output arrays */
    memset(out_real, 0, (size_t)(n_theta * n_omega) * sizeof(double));
    memset(out_imag, 0, (size_t)(n_theta * n_omega) * sizeof(double));

    /* ------------------------------------------------------------------ *
     * Main quadrature loops                                                *
     *                                                                      *
     * far_field_q[i,j] += rho * f(phi) * exp(-i*k*theta[i]*rho*cos(omega[j]-phi))
     *                                                                      *
     * where f(phi) = spline_real(phi) + i*spline_imag(phi).               *
     *                                                                      *
     * The i/j loops are parallelised with OpenMP when available.           *
     * ------------------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(dynamic)
#endif
    for (int i = 0; i < n_theta; ++i) {
        for (int j = 0; j < n_omega; ++j) {
            double acc_re = 0.0;
            double acc_im = 0.0;

            for (int p = 0; p < N_PHI; ++p) {
                const double phii    = phi[p];
                const double f_re    = spline_real[p];
                const double f_im    = spline_imag[p];
                const double om_phi  = omega[j] - phii;   /* omega[j] - phi */

                for (int ri = 0; ri < N_RHO; ++ri) {
                    const double rho   = rho_array[ri];
                    /*  phase = -k * theta[i] * rho * cos(omega[j] - phi)  */
                    const double phase = -k_val * theta[i] * rho * cos(om_phi);
                    const double cp    = cos(phase);
                    const double sp    = sin(phase);

                    /*  rho * (f_re + i*f_im) * (cp + i*sp)                */
                    acc_re += rho * (f_re * cp - f_im * sp);
                    acc_im += rho * (f_re * sp + f_im * cp);
                }
            }

            out_real[i * n_omega + j] = acc_re;
            out_imag[i * n_omega + j] = acc_im;
        }
    }

    /* Mirror the per-row print from the original Python */
    //for (int i = 0; i < n_theta; ++i)
        //printf("%d  theta = %g\n", i, theta[i]);
}
