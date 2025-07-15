#ifndef RKF78_H
#define RKF78_H

#include "solver.h"

real_t rkf78_solver(
    real_t y0, real_t t0, real_t t_end,
    real_t h_init, real_t tol,
    ODEFunction f, void* params);


#endif // !RKF78_H


