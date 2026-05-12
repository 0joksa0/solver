#ifndef BACKENDS_RK4_BACKEND_H
#define BACKENDS_RK4_BACKEND_H

#include "solver.h"

real_t rk4_backend_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx);

SolverStepResult rk4_backend_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx);

void rk4_backend_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

SolverStepResult rk4_backend_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

#endif
