#include "test_support.h"

Test(DefaultRuntime, VectorRunTracksObserverStatsAndLogger)
{
    SolverPluginManager manager;
    SolverDefaultRuntimeConfig runtime_config;
    ObserverProbe probe = {0};
    SolverStats stats;
    SolverLogger logger;
    HarmonicOscillatorParams params = { .omega = REAL(1.0) };
    real_t state[HARMONIC_OSCILLATOR_DIM] = { REAL(1.0), REAL(0.0) };
    SolverRunConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .t_end = REAL(2.0) * REAL(3.14159265358979323846L),
        .h_init = REAL(0.05),
        .tol = REAL(1e-8),
        .plugins = &manager
    };

    cr_assert_eq(solver_plugin_manager_init(&manager), 0);
    cr_assert_eq(solver_logger_init(&logger, "vector_steps.csv", "vector_summary.csv"), 0);

    runtime_config.observer = probe_observer;
    runtime_config.observer_ctx = &probe;
    runtime_config.stats = &stats;
    runtime_config.logger = &logger;

    cr_assert_eq(solver_default_runtime_init(&manager, &runtime_config), 0);

    solver_vector_run(&config, state, HARMONIC_OSCILLATOR_DIM, harmonic_oscillator_rhs, &params);
    solver_logger_close(&logger);

    cr_assert(probe.calls > 1);
    cr_assert(stats.n_steps > 0);
    cr_assert(fabsl((long double)probe.last_t - (long double)config.t_end) < 1e-6L);

    remove("vector_steps.csv");
    remove("vector_summary.csv");
}
