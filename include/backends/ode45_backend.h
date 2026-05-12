#ifndef BACKENDS_ODE45_BACKEND_H
#define BACKENDS_ODE45_BACKEND_H

#include "solver.h"

real_t ode45_backend_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx);

SolverStepResult ode45_backend_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx);

void ode45_backend_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

SolverStepResult ode45_backend_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

#endif
