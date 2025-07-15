
# Solver – Numerical ODE Solver Library in C

## Overview

**Solver** is a modular, header-configurable numerical library written in C that provides implementations of several Runge-Kutta methods for solving Ordinary Differential Equations (ODEs). It is designed for flexibility and numerical precision, with support for multiple real number representations.

This solver library was developed as part of a master's thesis focused on modeling plant growth through a digital twin approach. It serves as the numerical core for simulating physiological processes in plants by solving ordinary differential equations (ODEs) that describe biomass accumulation, resource allocation, and other dynamic behaviors over time.
---

## Implemented Features

- **RKF78 (Runge-Kutta-Fehlberg 7(8))** – adaptive step-size method with high accuracy.
- **ODE45** – a commonly used adaptive Runge-Kutta method (similar to MATLAB's ode45).
- **Classic Runge-Kutta** – standard fixed-step implementation.
- **Type-flexible real number arithmetic** – supports `float`, `double`, and `long double`.
- The library defines a macro **REAL(x)** to cast numeric literals to the correct real_t type (float, double, or long double) depending on the build configuration.


---

## Project Structure

```
solver/
├── include/         # Public headers (solver.h)
├── src/             # Core solver implementations
├── tests/           # Unit tests using Criterion framework
├── CMakeLists.txt   # CMake build configuration
├── solver.pc.in     # pkg-config template
└── readme.MD        # Project documentation
```

---

## Build Instructions

### Prerequisites

- CMake ≥ 3.15
- C compiler (GCC or Clang recommended)

### Build and Install

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=/usr/local
cmake --build build
sudo cmake --install build
```

---

## Usage

To use this library, you must define the real number precision before including the headers:

```c
#define SOLVER_REAL_AS double
#include <solver.h>
```

You can then call the solver like:

```c
real_t result = solve(SOLVER_ODE45, REAL(1e-8),
                      y0, t0, t_end, max_dt,
                      rhs_function, user_params);
```

### Linking

Depending on your precision, link against the appropriate variant:

- `-lsolver_float`
- `-lsolver_double`
- `-lsolver_longdouble`

---

## Current Status

- [x] Header-only interface for users
- [x] Adaptive solver infrastructure
- [x] CMake and pkg-config integration
- [x] Precision variants (float/double/long double)
- [x] Criterion-based test suite

---

## Future Work

- [ ] Error estimation and step rejection diagnostics
- [ ] Parallel execution strategies for large systems



