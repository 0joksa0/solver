#ifndef RUNGE_CUTTA_H
#define RUNGE_CUTTA_H
#include "solver.h"

real_t runge_kutta_4(real_t y_n, real_t t_n, real_t h, ODEFunction f, void* params);

void vector_runge_kutta_4_step(
    real_t t,
    real_t h,
    void* state,
    size_t state_size,
    RHS rhs,
    void* user_ctx);

void rk4_vector_solve(
    real_t* x,
    size_t n,
    real_t t0,
    real_t t_end,
    real_t h,
    RHS rhs,
    void* rhs_ctx,
    VectorObserver obs,
    void* obs_ctx
);

#endif
