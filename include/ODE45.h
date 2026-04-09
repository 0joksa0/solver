#ifndef ODE45_H
#define ODE45_H

#include "solver.h"

real_t ode45_solver(
    real_t y0, real_t t0, real_t t_end,
    real_t h_init, real_t tol,
    ODEFunction f, void* params);

void ode45_vector_step(
    real_t* state,
    size_t dim,
    real_t* t,
    real_t* h,
    real_t tol,
    RHS rhs,
    void* ctx);
void ode45_vector_solve(
    real_t* x, // in/out state
    size_t n, // dimenzija
    real_t t0,
    real_t t_end,
    real_t h_init,
    real_t rtol,
    real_t atol,
    RHS rhs,
    void* rhs_ctx,
    VectorObserver obs, // može NULL
    void* obs_ctx,
    SolverStats* stats,
    SolverLogger* logger);

#endif // !ODE45_H
