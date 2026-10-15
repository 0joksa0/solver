#include "test_support.h"

#include <framework/async_logger.h>
#include <stdio.h>
#include <string.h>

static int file_contains(const char* path, const char* needle)
{
    FILE* f = fopen(path, "r");
    char buf[512];
    int found = 0;

    if (!f) {
        return 0;
    }

    while (fgets(buf, sizeof(buf), f)) {
        if (strstr(buf, needle)) {
            found = 1;
            break;
        }
    }

    fclose(f);
    return found;
}

Test(AsyncLogger, FlushesQueuedStepAndSummaryToWrappedSolverLogger)
{
    SolverLogger sync_logger;
    AsyncSolverLogger async_logger;
    SolverStats stats;
    StepInfo info = {
        .dt_old = REAL(0.25),
        .dt_new = REAL(0.50),
        .err_norm = REAL(0.001),
        .accepted = 1
    };

    remove("async_steps.csv");
    remove("async_summary.csv");

    cr_assert_eq(solver_logger_init(&sync_logger, "async_steps.csv", "async_summary.csv"), 0);
    solver_stats_init(&stats);
    stats.n_steps = 1;
    stats.n_accept = 1;
    stats.dt_sum = REAL(0.25);
    stats.dt_min = REAL(0.25);
    stats.dt_max = REAL(0.25);
    stats.err_max = REAL(0.001);

    cr_assert_eq(async_solver_logger_init(&async_logger, &sync_logger, 4), 0);
    cr_assert_eq(async_solver_logger_start(&async_logger), 0);
    cr_assert_eq(async_solver_logger_log_step(&async_logger, REAL(1.0), &info), 0);
    cr_assert_eq(async_solver_logger_log_summary(&async_logger, &stats), 0);
    cr_assert_eq(async_solver_logger_flush(&async_logger), 0);
    async_solver_logger_stop(&async_logger);
    solver_logger_close(&sync_logger);

    cr_assert(file_contains("async_steps.csv", "1.000000"));
    cr_assert(file_contains("async_summary.csv", "1,1,0"));

    remove("async_steps.csv");
    remove("async_summary.csv");
}

Test(AsyncLogger, DefaultRuntimeCanWriteThroughAsyncLogger)
{
    SolverPluginManager manager;
    SolverDefaultRuntimeConfig runtime_config;
    SolverStats stats;
    SolverLogger sync_logger;
    AsyncSolverLogger async_logger;
    HarmonicOscillatorParams params = { .omega = REAL(1.0) };
    real_t state[HARMONIC_OSCILLATOR_DIM] = { REAL(1.0), REAL(0.0) };
    SolverRunConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .t_end = REAL(0.5),
        .h_init = REAL(0.05),
        .tol = REAL(1e-8),
        .plugins = &manager
    };

    remove("async_runtime_steps.csv");
    remove("async_runtime_summary.csv");

    cr_assert_eq(solver_plugin_manager_init(&manager), 0);
    cr_assert_eq(solver_logger_init(&sync_logger, "async_runtime_steps.csv", "async_runtime_summary.csv"), 0);
    cr_assert_eq(async_solver_logger_init(&async_logger, &sync_logger, 8), 0);
    cr_assert_eq(async_solver_logger_start(&async_logger), 0);

    runtime_config.observer = NULL;
    runtime_config.observer_ctx = NULL;
    runtime_config.stats = &stats;
    runtime_config.logger = NULL;
    runtime_config.async_logger = &async_logger;

    cr_assert_eq(solver_default_runtime_init(&manager, &runtime_config), 0);

    solver_vector_run(&config, state, HARMONIC_OSCILLATOR_DIM, harmonic_oscillator_rhs, &params);
    cr_assert_eq(async_solver_logger_flush(&async_logger), 0);
    async_solver_logger_stop(&async_logger);
    solver_logger_close(&sync_logger);

    cr_assert(file_contains("async_runtime_steps.csv", "accepted"));
    cr_assert(file_contains("async_runtime_summary.csv", "steps"));
    cr_assert(stats.n_steps > 0);

    remove("async_runtime_steps.csv");
    remove("async_runtime_summary.csv");
}
