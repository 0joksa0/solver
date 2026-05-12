#include "test_support.h"

Test(PluginStop, VectorStepCanBeStoppedByPlugin)
{
    SolverPluginManager manager;
    PluginTrace trace = {0};
    TracePluginCtx stop_ctx = {
        .trace = &trace,
        .marker = 9,
        .result = SOLVER_PLUGIN_STOP
    };
    SolverPlugin stop_plugin = {
        .name = "stopper",
        .plugin_ctx = &stop_ctx,
        .on_step_accepted = trace_on_step_accepted,
        .on_run_finish = trace_on_run_finish
    };
    SolverStepConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .h = REAL(0.1),
        .tol = REAL(1e-8),
        .plugins = &manager
    };
    real_t state[1] = { REAL(1.0) };
    real_t k = REAL(-0.5);
    SolverStepResult result;

    cr_assert_eq(solver_plugin_manager_init(&manager), 0);
    cr_assert_eq(solver_plugin_manager_add(&manager, &stop_plugin), 0);

    result = solver_vector_step(&config, state, 1, rhs_decay_vector, &k);

    cr_assert(result.accepted);
    cr_assert_eq(result.stopped_by_plugin, 1);
    cr_assert_eq(trace.count, 2);
}
