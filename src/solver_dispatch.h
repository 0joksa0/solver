#ifndef SOLVER_DISPATCH_H
#define SOLVER_DISPATCH_H

#include "solver.h"

real_t solver_dispatch_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx);

SolverStepResult solver_dispatch_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx);

void solver_dispatch_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

SolverStepResult solver_dispatch_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

#endif
