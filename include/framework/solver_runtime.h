#ifndef FRAMEWORK_SOLVER_RUNTIME_H
#define FRAMEWORK_SOLVER_RUNTIME_H

#include "framework/solver_plugin.h"

typedef struct SolverRuntime {
    SolverPluginManager* plugins;
    size_t step_index;
    size_t accepted_steps;
    size_t rejected_steps;
    int stopped_by_plugin;
} SolverRuntime;

void solver_runtime_init(
    SolverRuntime* runtime,
    SolverPluginManager* plugins);

SolverPluginResult solver_runtime_emit_run_start(
    SolverRuntime* runtime,
    real_t t,
    const real_t* state,
    size_t n);

SolverPluginResult solver_runtime_emit_step_attempt(
    SolverRuntime* runtime,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n);

SolverPluginResult solver_runtime_emit_step_accepted(
    SolverRuntime* runtime,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n,
    const StepInfo* info,
    size_t nfev_delta);

SolverPluginResult solver_runtime_emit_step_rejected(
    SolverRuntime* runtime,
    real_t t_start,
    real_t t_end,
    const real_t* state,
    size_t n,
    const StepInfo* info,
    size_t nfev_delta);

void solver_runtime_emit_run_finish(
    SolverRuntime* runtime,
    real_t t,
    const real_t* state,
    size_t n);

#endif
