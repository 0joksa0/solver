#include "framework/solver_runtime.h"

static SolverEvent solver_runtime_make_event(
    SolverRuntime* runtime,
    SolverEventType type,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n,
    const StepInfo* info,
    size_t nfev_delta,
    size_t step_index)
{
    SolverEvent event;

    event.type = type;
    event.t_start = t_start;
    event.t_end = t_end;
    event.state = state;
    event.state_dim = n;
    event.step_info = info;
    event.step_index = step_index;
    event.accepted_steps = runtime->accepted_steps;
    event.rejected_steps = runtime->rejected_steps;
    event.nfev_delta = nfev_delta;

    return event;
}

void solver_runtime_init(
    SolverRuntime* runtime,
    SolverPluginManager* plugins)
{
    runtime->plugins = plugins;
    runtime->step_index = 0;
    runtime->accepted_steps = 0;
    runtime->rejected_steps = 0;
    runtime->stopped_by_plugin = 0;
}

SolverPluginResult solver_runtime_emit_run_start(
    SolverRuntime* runtime,
    real_t t,
    const real_t* state,
    size_t n)
{
    SolverEvent event = solver_runtime_make_event(
        runtime, SOLVER_EVENT_RUN_START, t, t, state, n, NULL, 0, 0);

    if (solver_plugin_manager_emit(runtime->plugins, &event) == SOLVER_PLUGIN_STOP) {
        runtime->stopped_by_plugin = 1;
        return SOLVER_PLUGIN_STOP;
    }

    return SOLVER_PLUGIN_CONTINUE;
}

SolverPluginResult solver_runtime_emit_step_attempt(
    SolverRuntime* runtime,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n)
{
    SolverEvent event = solver_runtime_make_event(
        runtime, SOLVER_EVENT_STEP_ATTEMPT, t_start, t_end, state, n, NULL, 0, runtime->step_index + 1);

    if (solver_plugin_manager_emit(runtime->plugins, &event) == SOLVER_PLUGIN_STOP) {
        runtime->stopped_by_plugin = 1;
        return SOLVER_PLUGIN_STOP;
    }

    return SOLVER_PLUGIN_CONTINUE;
}

SolverPluginResult solver_runtime_emit_step_accepted(
    SolverRuntime* runtime,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n,
    const StepInfo* info,
    size_t nfev_delta)
{
    SolverEvent event;

    runtime->step_index++;
    runtime->accepted_steps++;
    event = solver_runtime_make_event(
        runtime, SOLVER_EVENT_STEP_ACCEPTED, t_start, t_end, state, n, info, nfev_delta, runtime->step_index);

    if (solver_plugin_manager_emit(runtime->plugins, &event) == SOLVER_PLUGIN_STOP) {
        runtime->stopped_by_plugin = 1;
        return SOLVER_PLUGIN_STOP;
    }

    return SOLVER_PLUGIN_CONTINUE;
}

SolverPluginResult solver_runtime_emit_step_rejected(
    SolverRuntime* runtime,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n,
    const StepInfo* info,
    size_t nfev_delta)
{
    SolverEvent event;

    runtime->step_index++;
    runtime->rejected_steps++;
    event = solver_runtime_make_event(
        runtime, SOLVER_EVENT_STEP_REJECTED, t_start, t_end, state, n, info, nfev_delta, runtime->step_index);

    if (solver_plugin_manager_emit(runtime->plugins, &event) == SOLVER_PLUGIN_STOP) {
        runtime->stopped_by_plugin = 1;
        return SOLVER_PLUGIN_STOP;
    }

    return SOLVER_PLUGIN_CONTINUE;
}

void solver_runtime_emit_run_finish(
    SolverRuntime* runtime,
    real_t t,
    const real_t* state,
    size_t n)
{
    SolverEvent event = solver_runtime_make_event(
        runtime, SOLVER_EVENT_RUN_FINISH, t, t, state, n, NULL, 0, runtime->step_index);

    (void)solver_plugin_manager_emit(runtime->plugins, &event);
}
