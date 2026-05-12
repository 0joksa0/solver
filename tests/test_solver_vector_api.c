#include "test_support.h"

Test(SolverVectorAPI, RunIntegratesHarmonicOscillatorWithoutPlugins)
{
    HarmonicOscillatorParams params = { .omega = REAL(1.0) };
    real_t state[HARMONIC_OSCILLATOR_DIM] = { REAL(1.0), REAL(0.0) };
    SolverRunConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .t_end = REAL(2.0) * REAL(3.14159265358979323846L),
        .h_init = REAL(0.05),
        .tol = REAL(1e-8),
        .plugins = NULL
    };

    solver_vector_run(&config, state, HARMONIC_OSCILLATOR_DIM, harmonic_oscillator_rhs, &params);

    cr_assert(fabsl((long double)state[0] - 1.0L) < 1e-4L);
    cr_assert(fabsl((long double)state[1]) < 1e-4L);
}

Test(SolverVectorAPI, StepAdvancesStateAndReturnsMetadata)
{
    SolverStepConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .h = REAL(0.1),
        .tol = REAL(1e-8),
        .plugins = NULL
    };
    real_t state[1] = { REAL(1.0) };
    real_t k = REAL(-0.5);
    SolverStepResult result = solver_vector_step(&config, state, 1, rhs_decay_vector, &k);

    cr_assert(result.accepted);
    cr_assert_lt((long double)state[0], 1.0L);
    cr_assert_gt((long double)result.t_out, 0.0L);
    cr_assert_gt((long double)result.h_used, 0.0L);
    cr_assert_gt((long double)result.h_next, 0.0L);
    cr_assert_eq(result.stopped_by_plugin, 0);
}
