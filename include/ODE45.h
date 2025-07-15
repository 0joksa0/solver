#ifndef ODE45_H
#define ODE45_H

#include "solver.h"

real_t ode45_solver(
    real_t y0, real_t t0, real_t t_end,
    real_t h_init, real_t tol,
    ODEFunction f, void* params);


#endif // !ODE45_H
