# HHG-Solids

**High-Order Harmonic Generation in 2D Materials from Wannier Tight-Binding Hamiltonians**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)](https://en.cppreference.com/w/cpp/17)
[![CMake](https://img.shields.io/badge/CMake-%3E%3D3.15-green.svg)](https://cmake.org/)
---

## Overview

`HHG-Solids` is a high-performance C++ simulation framework for computing **high-order harmonic generation (HHG)** spectra in two-dimensional solid-state materials. Starting from *ab initio*-derived Wannier tight-binding (TB) Hamiltonians (e.g., produced by DFT calculations followed by a wannierization process using [Wannier90](http://wannier.org/)), the code propagates the density matrix in both real space and k-space to obtain the time-dependent current density **J(t)**, as the fundamental observable for HHG.

The solver is written in C++17 and exposes a clean Python interface via `ctypes`, making it straightforward to set up and analyse simulations from Python scripts or Jupyter notebooks.

### Key capabilities

- Reads Wannier90-formatted tight-binding Hamiltonian files (`*_tb.dat`)
- Constructs Bloch Hamiltonians, Berry connections, and velocity matrices on an arbitrary **r-grid** and **k-grid**
- Propagates the **density matrix** in time under a strong laser field (SBE solver)
- Supports **arbitrary laser polarisation** and pulse shapes (in the future the laser field will be passed as an argument)
- Python wrapper (`pyswe.py`) for scripted workflows
- Shared-library build (`libwannier`) for embedding in external codes

---

## Repository Structure

```
HHG-Solids/
├── include/                  # C++ header files
│   ├── settings.hpp          # Simulation parameter container
│   ├── wannier_tb.hpp        # Wannier TB Hamiltonian reader/builder
│   ├── hamiltonian.hpp       # Bloch Hamiltonian & diagonalisation
│   ├── berry_connection.hpp  # Berry connection matrix elements
│   ├── velocity.hpp          # Velocity matrix elements
│   ├── efield.hpp            # Laser electric field
│   ├── swe.hpp               # Time propagation
│   ├── solver.hpp            # Real-space solver (solves one time step)
│   ├── solver_kspace.hpp     # k-space solver (solves one time step)
│   ├── rdm.hpp               # Reduced density matrix
│   ├── observable.hpp        # Observables (current, etc.) base class
│   ├── operator.hpp          # Operator base class
│   ├── grid.hpp              # r-point, k-point and temporal grid
│   ├── matrix_field.hpp      # Field-of-matrices utilities
│   ├── fftw_helper.hpp       # FFTW wrapper
│   ├── interpolation.hpp     # Interpolation utilities
│   ├── vec3_util.hpp         # 3-vector utilities
│   └── cwrapper.hpp          # C-linkage wrapper for Python ctypes
│
├── src/                      # C++ implementation files
│   ├── main.cpp             
│   ├── settings.cpp
│   ├── wannier_tb.cpp
│   ├── hamiltonian.cpp
│   ├── berry_connection.cpp
│   ├── velocity.cpp
│   ├── efield.cpp
│   ├── swe.cpp
│   ├── solver.cpp
│   ├── solver_kspace.cpp
│   ├── rdm.cpp
│   ├── observable.cpp
│   ├── operator.cpp
│   ├── grid.cpp
│   ├── matrix_field.cpp
│   ├── fftw_helper.cpp
│   ├── interpolation.cpp
│   └── cwrapper.cpp          # C-linkage exported symbols
│
├── scripts/                  # Utility / post-processing / example scripts
├── pyswe.py                  # Python ctypes interface
├── test_swe.py               # Minimal usage example
├── hmcase0_tb.dat            # Example Wannier TB Hamiltonian (hBN- / haldane-like)
├── CMakeLists.txt            # CMake build configuration
├── .gitignore
└── LICENSE                   # MIT licence
```

---

## Dependencies

| Library | Minimum version | Notes |
|---------|----------------|-------|
| CMake   | 3.15           | Build system |
| GCC / Clang | C++17 | Compiler |
| OpenBLAS | any recent | Dense linear algebra |
| LAPACK  | 3.x            | Eigenvalue solvers (`zheev`) |
| FFTW3   | 3.3            | Fast Fourier transforms |
| Python  | ≥ 3.8          | Optional, for Python interface |
| NumPy   | ≥ 1.20         | Optional, for Python interface |
| Mpi4py  | ≥ 4.0.0        | Optional, for `scripts/pol_scan.py` example |

### Installing dependencies (macOS / Homebrew)

```bash
brew install openblas lapack fftw
```

### Installing dependencies (Linux / apt)

```bash
sudo apt install libopenblas-dev liblapacke-dev libfftw3-dev libomp-dev
```

---

## Building

```bash
# 1. Clone the repository
git clone https://github.com/rodrigomarher/HHG-Solids.git
cd HHG-Solids

# 2. Edit library paths in CMakeLists.txt if they differ from the defaults
#    (see LAPACK_DIR, OPENBLAS_DIR, FFTW_DIR variables)

# 3. Configure and build
mkdir build && cd build
cmake ..
cmake --build . -j$(nproc)
```

This produces:
- `build/libwannier.so` (or `.dylib` on macOS) — shared library
- `build/main` — standalone executable

> **Note (macOS):** The `CMakeLists.txt` currently hard-codes Homebrew Cellar paths. Either update them to match your installation or export the appropriate `CMAKE_PREFIX_PATH`.

---

## Running a simulation

### Python (recommended)

```python
import numpy as np
from pyswe import Settings, SWE

param = {
    "path_lib":    "build/libwannier.dylib",  # or .so on Linux
    "path_tb":     "hmcase0_tb.dat",           # Wannier90 _tb.dat file
    "nr1": 200, "nr2": 200, "nr3": 1,          # k-grid dimensions
    "tmax":        90.0,                        # total propagation time (a.u.)
    "dt":          21.97e-3,                    # time step (a.u.)
    "intensity":   1e12,                        # peak intensity (W/cm²)
    "lambda":      3000.0,                      # central wavelength (nm)
    "tmax_field":  80.0,                        # field duration (a.u.)
    "pol_vec":     np.array([0.0, 1.0, 0.0]),  # polarisation direction
    "phi_vec":     np.array([0.0, 0.0, 0.0]),  # carrier-envelope phase vector
}

settings = Settings(param)
swe      = SWE(settings)
swe.run_simulation()

t, jx, jy, jz = swe.get_current()  # retrieve time-resolved current
swe.delete()                        # free C++ heap memory
```

See `test_swe.py` for a complete minimal example.

### Standalone executable

```bash
./build/main
```

The executable reads its configuration from the same `Settings` object initialised inside `main.cpp`. Edit `src/main.cpp` to adjust parameters for standalone runs.

---

## Input file format

`hmcase0_tb.dat` follows the **Wannier90 `_tb.dat`** convention:

```
<comment line>
<num_wann>
<nrpts>
<degeneracy weights …>
<R1 R2 R3  m  n  Re(H_mn)  Im(H_mn)>
…
```

Any Wannier90-compatible tight-binding Hamiltonian can be used as input.

---

## Output

The primary output is the time-resolved **current density** returned by `get_current()`:

| Array | Type | Description |
|-------|------|-------------|
| `t`   | `float64[nt]`    | Time grid (a.u.) |
| `jx`  | `complex128[nt]` | x-component of J(t) |
| `jy`  | `complex128[nt]` | y-component of J(t) |
| `jz`  | `complex128[nt]` | z-component of J(t) |

The **HHG spectrum** is obtained by Fourier-transforming the current:

```python
import numpy as np

dt = t[1] - t[0]
freq = np.fft.rfftfreq(len(t), d=dt)
spectrum = np.abs(np.fft.rfft(jy.real))**2
```

---

## Citation

If you use this code in published work, please cite the repository:

```bibtex
@software{HHG-Solids,
  author  = {Martín Hernández, Rodrigo},
  title   = {{HHG-Solids}: High-Order Harmonic Generation in 2D Materials
             from Wannier Tight-Binding Hamiltonians},
  year    = {2026},
  url     = {https://github.com/rodrigomarher/HHG-Solids},
  license = {MIT}
}
```

---

## Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or pull request on GitHub.

---

## License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.
