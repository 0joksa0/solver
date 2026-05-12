#include "test_support.h"

Test(SolverScalarAPI, RunComputesODE45Solution)
{
    SolverRunConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .t_end = REAL(2.0),
        .h_init = REAL(0.2),
        .tol = REAL(1e-8),
        .plugins = NULL
    };
    real_t expected = REAL((2.0 * 2.0 * 2.0) / 3.0);
    real_t y_end = solver_scalar_run(&config, REAL(0.0), rhs_poly, NULL);

    cr_assert(RABS(y_end - expected) < REAL(1e-6));
}

Test(SolverScalarAPI, StepReturnsAdaptiveMetadata)
{
    SolverStepConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .h = REAL(0.1),
        .tol = REAL(1e-8),
        .plugins = NULL
    };
    real_t y = REAL(1.0);
    real_t k = REAL(-0.5);
    SolverStepResult result = solver_scalar_step(&config, &y, rhs_decay_scalar, &k);

    cr_assert(result.accepted);
    cr_assert_gt((long double)result.t_out, 0.0L);
    cr_assert_gt((long double)result.h_used, 0.0L);
    cr_assert_gt((long double)result.h_next, 0.0L);
    cr_assert_eq(result.stopped_by_plugin, 0);
}
