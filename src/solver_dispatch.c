#include "solver_dispatch.h"

#include "backends/ode45_backend.h"
#include "backends/rk4_backend.h"
#include "backends/rkf78_backend.h"

#include <stdio.h>
#include <stdlib.h>

static void solver_dispatch_panic(SolverType type, const char* op)
{
    fprintf(stderr, "[solver_dispatch] unsupported solver type %d for %s\n", (int)type, op);
    exit(EXIT_FAILURE);
}

real_t solver_dispatch_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx)
{
    switch (config->type) {
    case SOLVER_RK4:
        return rk4_backend_scalar_run(config, y0, rhs, rhs_ctx);
    case SOLVER_ODE45:
        return ode45_backend_scalar_run(config, y0, rhs, rhs_ctx);
    case SOLVER_RKF78:
        return rkf78_backend_scalar_run(config, y0, rhs, rhs_ctx);
    default:
        solver_dispatch_panic(config->type, "scalar_run");
        return y0;
    }
}

SolverStepResult solver_dispatch_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx)
{
    switch (config->type) {
    case SOLVER_RK4:
        return rk4_backend_scalar_step(config, y, rhs, rhs_ctx);
    case SOLVER_ODE45:
        return ode45_backend_scalar_step(config, y, rhs, rhs_ctx);
    case SOLVER_RKF78:
        return rkf78_backend_scalar_step(config, y, rhs, rhs_ctx);
    default:
        solver_dispatch_panic(config->type, "scalar_step");
        return (SolverStepResult){0};
    }
}

void solver_dispatch_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    switch (config->type) {
    case SOLVER_RK4:
        rk4_backend_vector_run(config, state, n, rhs, rhs_ctx);
        return;
    case SOLVER_ODE45:
        ode45_backend_vector_run(config, state, n, rhs, rhs_ctx);
        return;
    case SOLVER_RKF78:
        rkf78_backend_vector_run(config, state, n, rhs, rhs_ctx);
        return;
    default:
        solver_dispatch_panic(config->type, "vector_run");
        return;
    }
}

SolverStepResult solver_dispatch_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    switch (config->type) {
    case SOLVER_RK4:
        return rk4_backend_vector_step(config, state, n, rhs, rhs_ctx);
    case SOLVER_ODE45:
        return ode45_backend_vector_step(config, state, n, rhs, rhs_ctx);
    case SOLVER_RKF78:
        return rkf78_backend_vector_step(config, state, n, rhs, rhs_ctx);
    default:
        solver_dispatch_panic(config->type, "vector_step");
        return (SolverStepResult){0};
    }
}
