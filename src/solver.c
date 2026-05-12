#include "solver.h"

#include "solver_dispatch.h"

real_t solver_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx)
{
    return solver_dispatch_scalar_run(config, y0, rhs, rhs_ctx);
}

SolverStepResult solver_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx)
{
    return solver_dispatch_scalar_step(config, y, rhs, rhs_ctx);
}

void solver_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    solver_dispatch_vector_run(config, state, n, rhs, rhs_ctx);
}

SolverStepResult solver_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    return solver_dispatch_vector_step(config, state, n, rhs, rhs_ctx);
}
