"""
prop_q_ctypes.py
================
ctypes interface to the ``prop_q`` shared library compiled from ``prop_q.c``.

Scipy ``splrep`` produces a (t, c, k) tuple whose knots *t* and coefficients
*c* are passed directly to the C ``bspline_eval`` / ``prop_q_c`` routines.

Expected workflow
-----------------
    from scipy.interpolate import splrep
    import numpy as np
    from prop_q_ctypes import prop_q

    phi_data = np.linspace(0, 2*np.pi, 256)

    # Fit separate splines for real and imaginary parts
    tck_real = splrep(phi_data, your_complex_field.real, k=3, s=0)
    tck_imag = splrep(phi_data, your_complex_field.imag, k=3, s=0)
    coefs_interp = (tck_real, tck_imag)

    theta = np.linspace(0, np.pi, 50)
    omega = np.linspace(0, 2*np.pi, 64)

    result = prop_q(coefs_interp, angles=None, freq=None, q=2,
                    theta=theta, omega=omega)
    # result.shape == (50, 64), dtype complex128
"""

import ctypes
import os
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Library loading
# ---------------------------------------------------------------------------

def _find_library() -> ctypes.CDLL:
    """Search for libprop_q in common build locations."""
    here = Path(__file__).resolve().parent

    if sys.platform == "win32":
        candidates = ["prop_q.dll"]
    elif sys.platform == "darwin":
        candidates = ["libprop_q.dylib"]
    else:
        candidates = ["libprop_q.so"]

    search_roots = [
        here,                        # same directory as this script (CMake puts it here)
        here / "build",
        here / "build" / "lib",
        here / "build" / "Release",
        here / "build" / "Debug",
        Path("/usr/local/lib"),
        Path("/usr/lib"),
    ]

    for root in search_roots:
        for name in candidates:
            path = root / name
            if path.exists():
                return ctypes.CDLL(str(path))

    tried = "\n  ".join(str(r / c) for r in search_roots for c in candidates)
    raise FileNotFoundError(
        f"Could not find the prop_q shared library.\n"
        f"Tried:\n  {tried}\n\n"
        f"Build it first:\n"
        f"  mkdir -p build && cd build && cmake .. && cmake --build ."
    )


_lib = _find_library()

# ---------------------------------------------------------------------------
# Function signatures
# ---------------------------------------------------------------------------

_dbl_p = ctypes.POINTER(ctypes.c_double)

# double bspline_eval(const double *t, int n_t, const double *c, double x)
_lib.bspline_eval.restype  = ctypes.c_double
_lib.bspline_eval.argtypes = [_dbl_p, ctypes.c_int, _dbl_p, ctypes.c_double]

# void prop_q_c(t_real, n_t_real, c_real,
#               t_imag, n_t_imag, c_imag,
#               q,
#               theta, n_theta,
#               omega, n_omega,
#               out_real, out_imag)
_lib.prop_q_c.restype  = None
_lib.prop_q_c.argtypes = [
    _dbl_p, ctypes.c_int, _dbl_p,   # real spline  (t, n_t, c)
    _dbl_p, ctypes.c_int, _dbl_p,   # imag spline  (t, n_t, c)
    ctypes.c_double,                 # q
    _dbl_p, ctypes.c_int,           # theta, n_theta
    _dbl_p, ctypes.c_int,           # omega, n_omega
    _dbl_p,                          # out_real
    _dbl_p,                          # out_imag
]

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------

def _ptr(arr: np.ndarray) -> ctypes.POINTER(ctypes.c_double):
    """Return a ctypes double pointer to a contiguous float64 array."""
    assert arr.dtype == np.float64 and arr.flags["C_CONTIGUOUS"]
    return arr.ctypes.data_as(_dbl_p)


def _c64(arr: np.ndarray) -> np.ndarray:
    """Ensure contiguous float64 array."""
    return np.ascontiguousarray(arr, dtype=np.float64)

# ---------------------------------------------------------------------------
# Public Python API
# ---------------------------------------------------------------------------

def bspline_eval_py(tck: tuple, x: float) -> float:
    """
    Evaluate a single cubic B-spline at scalar *x*.

    Parameters
    ----------
    tck : tuple
        ``(t, c, k)`` as returned by ``scipy.interpolate.splrep``.
        *k* must be 3 (cubic).
    x : float
        Evaluation point.

    Returns
    -------
    float
    """
    t_arr, c_arr, deg = tck
    if deg != 3:
        raise ValueError(f"Only cubic (k=3) splines are supported, got k={deg}")
    t = _c64(np.asarray(t_arr))
    c = _c64(np.asarray(c_arr))
    return _lib.bspline_eval(_ptr(t), len(t), _ptr(c), float(x))


def prop_q(
    coefs_interp: tuple,
    angles,        # kept for API compatibility with Python original; unused
    freq,          # kept for API compatibility with Python original; unused
    q: float,
    theta: np.ndarray = None,
    omega: np.ndarray = None,
) -> np.ndarray:
    """
    Far-field propagation integral for mode *q*.

    Parameters
    ----------
    coefs_interp : tuple of two (t, c, k) tuples
        ``coefs_interp[0]`` – spline for the **real** part of the field.
        ``coefs_interp[1]`` – spline for the **imaginary** part of the field.
        Both must be cubic (``k=3``) and produced by
        ``scipy.interpolate.splrep`` over phi in [0, 2*pi].
    angles : ignored
        Kept for drop-in compatibility with the original Python signature.
    freq : ignored
        Kept for drop-in compatibility with the original Python signature.
    q : float
        Mode index.  Determines the wave-number ``k = q * 2*pi / 3``.
    theta : ndarray, shape (N,)
        Angular positions (radians).
    omega : ndarray, shape (M,)
        Azimuthal angles (radians).

    Returns
    -------
    far_field_q : ndarray, shape (N, M), dtype complex128
    """
    if theta is None or omega is None:
        raise ValueError("theta and omega must be provided.")

    # Unpack spline tuples
    t_re, c_re, deg_re = coefs_interp[0]
    t_im, c_im, deg_im = coefs_interp[1]
    if deg_re != 3 or deg_im != 3:
        raise ValueError(
            f"Only cubic (k=3) splines are supported; "
            f"got k_real={deg_re}, k_imag={deg_im}"
        )

    t_re = _c64(np.asarray(t_re))
    c_re = _c64(np.asarray(c_re))
    t_im = _c64(np.asarray(t_im))
    c_im = _c64(np.asarray(c_im))
    theta = _c64(np.asarray(theta))
    omega = _c64(np.asarray(omega))

    n_theta = int(theta.shape[0])
    n_omega = int(omega.shape[0])

    out_real = np.zeros(n_theta * n_omega, dtype=np.float64)
    out_imag = np.zeros(n_theta * n_omega, dtype=np.float64)

    _lib.prop_q_c(
        _ptr(t_re), len(t_re), _ptr(c_re),
        _ptr(t_im), len(t_im), _ptr(c_im),
        float(q),
        _ptr(theta), n_theta,
        _ptr(omega),  n_omega,
        _ptr(out_real),
        _ptr(out_imag),
    )

    return (
        out_real.reshape(n_theta, n_omega)
        + 1j * out_imag.reshape(n_theta, n_omega)
    )
