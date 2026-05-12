#include "backends/rk4_backend.h"

#include "framework/solver_runtime.h"

#include <stdio.h>
#include <stdlib.h>

typedef struct RK4ScalarAdapterCtx {
    ODEFunction rhs;
    void* rhs_ctx;
} RK4ScalarAdapterCtx;

typedef enum RK4Mode {
    RK4_MODE_RUN = 0,
    RK4_MODE_STEP = 1
} RK4Mode;

static void rk4_scalar_rhs_adapter(
    real_t t,
    const real_t* state,
    real_t* dstate,
    size_t n,
    void* ctx)
{
    RK4ScalarAdapterCtx* adapter = (RK4ScalarAdapterCtx*)ctx;
    (void)n;

    dstate[0] = adapter->rhs(t, state[0], adapter->rhs_ctx);
}

static void* rk4_xcalloc(size_t count, size_t size)
{
    void* ptr = calloc(count, size);

    if (!ptr) {
        fprintf(stderr, "[rk4_backend] allocation failed\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

static SolverStepResult rk4_core_run(
    RK4Mode mode,
    const SolverPluginManager* plugins,
    real_t* state,
    size_t n,
    real_t t0,
    real_t t_end,
    real_t h,
    VectorRHS rhs,
    void* rhs_ctx)
{
    SolverRuntime runtime;
    SolverStepResult result;
    real_t t = t0;
    real_t* k1;
    real_t* k2;
    real_t* k3;
    real_t* k4;
    real_t* yt;

    result.t_out = t0;
    result.h_used = REAL(0.0);
    result.h_next = h;
    result.accepted = 0;
    result.stopped_by_plugin = 0;

    if (n == 0) {
        return result;
    }

    solver_runtime_init(&runtime, (SolverPluginManager*)plugins);

    k1 = rk4_xcalloc(n, sizeof(real_t));
    k2 = rk4_xcalloc(n, sizeof(real_t));
    k3 = rk4_xcalloc(n, sizeof(real_t));
    k4 = rk4_xcalloc(n, sizeof(real_t));
    yt = rk4_xcalloc(n, sizeof(real_t));

    if (solver_runtime_emit_run_start(&runtime, t, state, n) == SOLVER_PLUGIN_STOP) {
        goto finish;
    }

    while (mode == RK4_MODE_STEP || t < t_end) {
        StepInfo info;
        real_t step_h = h;
        real_t step_t0 = t;
        size_t i;

        if (mode == RK4_MODE_RUN && t + step_h > t_end) {
            step_h = t_end - t;
        }

        if (step_h <= REAL(0.0)) {
            break;
        }

        if (solver_runtime_emit_step_attempt(&runtime, step_t0, step_t0 + step_h, state, n) == SOLVER_PLUGIN_STOP) {
            break;
        }

        rhs(t, state, k1, n, rhs_ctx);

        for (i = 0; i < n; ++i) {
            yt[i] = state[i] + step_h * REAL(0.5) * k1[i];
        }
        rhs(t + step_h * REAL(0.5), yt, k2, n, rhs_ctx);

        for (i = 0; i < n; ++i) {
            yt[i] = state[i] + step_h * REAL(0.5) * k2[i];
        }
        rhs(t + step_h * REAL(0.5), yt, k3, n, rhs_ctx);

        for (i = 0; i < n; ++i) {
            yt[i] = state[i] + step_h * k3[i];
        }
        rhs(t + step_h, yt, k4, n, rhs_ctx);

        for (i = 0; i < n; ++i) {
            state[i] += (step_h / REAL(6.0)) * (k1[i] + REAL(2.0) * k2[i] + REAL(2.0) * k3[i] + k4[i]);
        }

        t += step_h;
        info.err_norm = REAL(0.0);
        info.dt_old = step_h;
        info.dt_new = step_h;
        info.accepted = 1;

        result.t_out = t;
        result.h_used = step_h;
        result.h_next = step_h;
        result.accepted = 1;

        if (solver_runtime_emit_step_accepted(&runtime, step_t0, t, state, n, &info, 4) == SOLVER_PLUGIN_STOP) {
            break;
        }

        if (mode == RK4_MODE_STEP) {
            break;
        }
    }

finish:
    solver_runtime_emit_run_finish(&runtime, t, state, n);
    result.stopped_by_plugin = runtime.stopped_by_plugin;

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(yt);

    return result;
}

real_t rk4_backend_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx)
{
    RK4ScalarAdapterCtx adapter;
    real_t state[1];

    state[0] = y0;
    adapter.rhs = rhs;
    adapter.rhs_ctx = rhs_ctx;

    rk4_core_run(
        RK4_MODE_RUN,
        config->plugins,
        state,
        1,
        config->t0,
        config->t_end,
        config->h_init,
        rk4_scalar_rhs_adapter,
        &adapter);

    return state[0];
}

SolverStepResult rk4_backend_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx)
{
    RK4ScalarAdapterCtx adapter;
    real_t state[1];
    SolverStepResult result;

    state[0] = *y;
    adapter.rhs = rhs;
    adapter.rhs_ctx = rhs_ctx;

    result = rk4_core_run(
        RK4_MODE_STEP,
        config->plugins,
        state,
        1,
        config->t0,
        config->t0 + config->h,
        config->h,
        rk4_scalar_rhs_adapter,
        &adapter);

    *y = state[0];
    return result;
}

void rk4_backend_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    (void)rk4_core_run(
        RK4_MODE_RUN,
        config->plugins,
        state,
        n,
        config->t0,
        config->t_end,
        config->h_init,
        rhs,
        rhs_ctx);
}

SolverStepResult rk4_backend_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    return rk4_core_run(
        RK4_MODE_STEP,
        config->plugins,
        state,
        n,
        config->t0,
        config->t0 + config->h,
        config->h,
        rhs,
        rhs_ctx);
}
