#ifndef PROP_Q_H
#define PROP_Q_H

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Evaluate a cubic (degree-3) B-spline at a single point x using de Boor's algorithm.
 *
 * Matches scipy.interpolate.splev(x, (t, c, 3)).
 *
 * @param t     Knot vector (length n_t).  Must include repeated boundary knots as
 *              produced by scipy.interpolate.splrep / UnivariateSpline.
 * @param n_t   Number of knots.
 * @param c     B-spline coefficients (length n_t; only the first n_t-4 are used).
 * @param x     Evaluation point.  Clamped to [t[3], t[n_t-4]] if outside range.
 * @return      Spline value at x.
 */
double bspline_eval(const double *t, int n_t, const double *c, double x);

/**
 * Propagate far-field integral for mode q.
 *
 * C equivalent of the Python prop_q() function.  Computes:
 *
 *   far_field_q[i,j] = sum_{phi} sum_{rho} rho * f(phi)
 *                      * exp(-i*k*theta[i]*rho*cos(omega[j]-phi))
 *
 * where f(phi) = spline_real(phi) + i*spline_imag(phi),
 * phi  in linspace(0, 2*pi, 128),
 * rho  in arange(22, 35, 0.5)  (26 values),
 * k    = q * 2*pi / 3.
 *
 * @param t_real   Knot vector for the real part spline.
 * @param n_t_real Number of knots (real part).
 * @param c_real   Coefficients for the real part spline.
 * @param t_imag   Knot vector for the imaginary part spline.
 * @param n_t_imag Number of knots (imaginary part).
 * @param c_imag   Coefficients for the imaginary part spline.
 * @param q        Mode index (determines wave-number k = q*2*pi/3).
 * @param theta    Array of theta values (length n_theta).
 * @param n_theta  Number of theta values.
 * @param omega    Array of omega values (length n_omega).
 * @param n_omega  Number of omega values.
 * @param out_real Output real part,      flat row-major array [n_theta * n_omega].
 * @param out_imag Output imaginary part, flat row-major array [n_theta * n_omega].
 */
void prop_q_c(
    const double *t_real, int n_t_real, const double *c_real,
    const double *t_imag, int n_t_imag, const double *c_imag,
    double q,
    const double *theta, int n_theta,
    const double *omega, int n_omega,
    double *out_real,
    double *out_imag
);

#ifdef __cplusplus
}
#endif

#endif /* PROP_Q_H */
