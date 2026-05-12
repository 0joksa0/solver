#include "backends/rkf78_backend.h"

#include "framework/solver_runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define RKF78_STAGES 13

typedef struct RKF78ScalarAdapterCtx {
    ODEFunction rhs;
    void* rhs_ctx;
} RKF78ScalarAdapterCtx;

typedef enum RKF78Mode {
    RKF78_MODE_RUN = 0,
    RKF78_MODE_STEP = 1
} RKF78Mode;

static const real_t rkf78_c[RKF78_STAGES] = {
    REAL(0.0),
    REAL(2.0) / REAL(27.0),
    REAL(1.0) / REAL(9.0),
    REAL(1.0) / REAL(6.0),
    REAL(5.0) / REAL(12.0),
    REAL(0.5),
    REAL(5.0) / REAL(6.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0)
};

static const real_t rkf78_a[RKF78_STAGES][RKF78_STAGES] = {
    { 0 },
    { REAL(2.0) / REAL(27.0) },
    { REAL(1.0) / REAL(36.0), REAL(1.0) / REAL(12.0) },
    { REAL(1.0) / REAL(24.0), REAL(0.0), REAL(1.0) / REAL(8.0) },
    { REAL(5.0) / REAL(12.0), REAL(0.0), -REAL(25.0) / REAL(16.0), REAL(25.0) / REAL(16.0) },
    { REAL(1.0) / REAL(20.0), REAL(0.0), REAL(0.0), REAL(1.0) / REAL(4.0), REAL(1.0) / REAL(5.0) },
    { -REAL(25.0) / REAL(108.0), REAL(0.0), REAL(0.0), REAL(125.0) / REAL(108.0), -REAL(65.0) / REAL(27.0), REAL(125.0) / REAL(54.0) },
    { REAL(31.0) / REAL(300.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(61.0) / REAL(225.0), -REAL(2.0) / REAL(9.0), REAL(13.0) / REAL(900.0) },
    { REAL(2.0), REAL(0.0), REAL(0.0), -REAL(53.0) / REAL(6.0), REAL(704.0) / REAL(45.0), -REAL(107.0) / REAL(9.0), REAL(67.0) / REAL(90.0), REAL(3.0) },
    { -REAL(91.0) / REAL(108.0), REAL(0.0), REAL(0.0), REAL(23.0) / REAL(108.0), -REAL(976.0) / REAL(135.0), REAL(311.0) / REAL(54.0), -REAL(19.0) / REAL(60.0), REAL(17.0) / REAL(6.0), -REAL(1.0) / REAL(12.0) },
    { REAL(2383.0) / REAL(4100.0), REAL(0.0), REAL(0.0), -REAL(341.0) / REAL(164.0), REAL(4496.0) / REAL(1025.0), -REAL(301.0) / REAL(82.0), REAL(2133.0) / REAL(4100.0), REAL(45.0) / REAL(82.0), REAL(45.0) / REAL(164.0), REAL(18.0) / REAL(41.0) },
    { REAL(3.0) / REAL(205.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), -REAL(6.0) / REAL(41.0), -REAL(3.0) / REAL(205.0), -REAL(3.0) / REAL(41.0), REAL(3.0) / REAL(41.0), REAL(6.0) / REAL(41.0), REAL(0.0) },
    { -REAL(1777.0) / REAL(4100.0), REAL(0.0), REAL(0.0), -REAL(341.0) / REAL(164.0), REAL(4496.0) / REAL(1025.0), -REAL(289.0) / REAL(82.0), REAL(2193.0) / REAL(4100.0), REAL(51.0) / REAL(82.0), REAL(33.0) / REAL(164.0), REAL(12.0) / REAL(41.0), REAL(0.0), REAL(1.0) }
};

static const real_t rkf78_b7[RKF78_STAGES] = {
    REAL(41.0) / REAL(840.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(34.0) / REAL(105.0), REAL(9.0) / REAL(35.0),
    REAL(9.0) / REAL(35.0), REAL(9.0) / REAL(280.0), REAL(9.0) / REAL(280.0), REAL(41.0) / REAL(840.0), REAL(0.0), REAL(0.0)
};

static const real_t rkf78_b8[RKF78_STAGES] = {
    REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(34.0) / REAL(105.0), REAL(9.0) / REAL(35.0),
    REAL(9.0) / REAL(35.0), REAL(9.0) / REAL(280.0), REAL(9.0) / REAL(280.0), REAL(0.0), REAL(0.0), REAL(0.0)
};

static void rkf78_scalar_rhs_adapter(
    real_t t,
    const real_t* state,
    real_t* dstate,
    size_t n,
    void* ctx)
{
    RKF78ScalarAdapterCtx* adapter = (RKF78ScalarAdapterCtx*)ctx;
    (void)n;

    dstate[0] = adapter->rhs(t, state[0], adapter->rhs_ctx);
}

static void* rkf78_xcalloc(size_t count, size_t size)
{
    void* ptr = calloc(count, size);

    if (!ptr) {
        fprintf(stderr, "[rkf78_backend] allocation failed\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

static void rkf78_map_tol(real_t tol, real_t* rtol, real_t* atol)
{
    real_t safe_tol = RMAX(tol, REAL_EPSILON);

    *rtol = safe_tol;
    *atol = safe_tol;
}

static SolverStepResult rkf78_core_run(
    RKF78Mode mode,
    const SolverPluginManager* plugins,
    real_t* state,
    size_t n,
    real_t t0,
    real_t t_end,
    real_t h_init,
    real_t tol,
    VectorRHS rhs,
    void* rhs_ctx)
{
    SolverRuntime runtime;
    SolverStepResult result;
    real_t t = t0;
    real_t h = h_init;
    real_t rtol;
    real_t atol;
    const real_t safety = REAL(0.9);
    const real_t min_scale = REAL(0.1);
    const real_t max_scale = REAL(4.0);
    real_t* k_storage;
    real_t* yt;
    real_t* y7;
    real_t* y8;
    size_t stage;

    result.t_out = t0;
    result.h_used = REAL(0.0);
    result.h_next = h_init;
    result.accepted = 0;
    result.stopped_by_plugin = 0;

    if (n == 0) {
        return result;
    }

    rkf78_map_tol(tol, &rtol, &atol);
    solver_runtime_init(&runtime, (SolverPluginManager*)plugins);

    k_storage = rkf78_xcalloc((size_t)RKF78_STAGES * n, sizeof(real_t));
    yt = rkf78_xcalloc(n, sizeof(real_t));
    y7 = rkf78_xcalloc(n, sizeof(real_t));
    y8 = rkf78_xcalloc(n, sizeof(real_t));

    if (solver_runtime_emit_run_start(&runtime, t, state, n) == SOLVER_PLUGIN_STOP) {
        goto finish;
    }

    while (mode == RKF78_MODE_STEP || t < t_end) {
        real_t err_sum = REAL(0.0);
        real_t step_t0 = t;
        real_t step_h = h;
        real_t err;
        real_t scale;
        real_t next_h;
        StepInfo info;
        size_t i;

        if (mode == RKF78_MODE_RUN && t + step_h > t_end) {
            step_h = t_end - t;
        }

        if (step_h <= REAL(0.0)) {
            break;
        }

        if (solver_runtime_emit_step_attempt(&runtime, step_t0, step_t0 + step_h, state, n) == SOLVER_PLUGIN_STOP) {
            break;
        }

        for (stage = 0; stage < (size_t)RKF78_STAGES; ++stage) {
            real_t* k_stage = k_storage + stage * n;

            if (stage == 0) {
                rhs(t, state, k_stage, n, rhs_ctx);
            } else {
                for (i = 0; i < n; ++i) {
                    real_t sum = REAL(0.0);
                    size_t j;

                    for (j = 0; j < stage; ++j) {
                        sum += rkf78_a[stage][j] * k_storage[j * n + i];
                    }
                    yt[i] = state[i] + sum;
                }
                rhs(t + rkf78_c[stage] * step_h, yt, k_stage, n, rhs_ctx);
            }

            for (i = 0; i < n; ++i) {
                k_stage[i] *= step_h;
            }
        }

        for (i = 0; i < n; ++i) {
            real_t scale_i;
            real_t e;

            y7[i] = state[i];
            y8[i] = state[i];

            for (stage = 0; stage < (size_t)RKF78_STAGES; ++stage) {
                y7[i] += rkf78_b7[stage] * k_storage[stage * n + i];
                y8[i] += rkf78_b8[stage] * k_storage[stage * n + i];
            }

            scale_i = atol + rtol * RMAX(RABS(state[i]), RABS(y8[i]));
            e = (y8[i] - y7[i]) / scale_i;
            err_sum += e * e;
        }

        err = RSQRT(err_sum / REAL(n));
        scale = safety * RPOW(REAL(1.0) / RMAX(err, REAL_EPSILON), REAL(1.0) / REAL(8.0));
        scale = RMAX(min_scale, RMIN(scale, max_scale));
        next_h = step_h * scale;

        info.err_norm = err;
        info.dt_old = step_h;
        info.dt_new = next_h;
        info.accepted = (err <= REAL(1.0)) || (RABS(step_h) <= REAL_EPSILON);

        h = next_h;
        result.h_next = h;

        if (info.accepted) {
            memcpy(state, y8, n * sizeof(real_t));
            t += step_h;

            result.t_out = t;
            result.h_used = step_h;
            result.accepted = 1;

            if (solver_runtime_emit_step_accepted(&runtime, step_t0, t, state, n, &info, RKF78_STAGES) == SOLVER_PLUGIN_STOP) {
                break;
            }

            if (mode == RKF78_MODE_STEP) {
                break;
            }

            continue;
        }

        result.t_out = t;
        result.h_used = REAL(0.0);
        if (solver_runtime_emit_step_rejected(&runtime, step_t0, step_t0 + step_h, state, n, &info, RKF78_STAGES) == SOLVER_PLUGIN_STOP) {
            break;
        }
    }

finish:
    solver_runtime_emit_run_finish(&runtime, t, state, n);
    result.stopped_by_plugin = runtime.stopped_by_plugin;

    free(k_storage);
    free(yt);
    free(y7);
    free(y8);

    return result;
}

real_t rkf78_backend_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx)
{
    RKF78ScalarAdapterCtx adapter;
    real_t state[1];

    state[0] = y0;
    adapter.rhs = rhs;
    adapter.rhs_ctx = rhs_ctx;

    rkf78_core_run(
        RKF78_MODE_RUN,
        config->plugins,
        state,
        1,
        config->t0,
        config->t_end,
        config->h_init,
        config->tol,
        rkf78_scalar_rhs_adapter,
        &adapter);

    return state[0];
}

SolverStepResult rkf78_backend_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx)
{
    RKF78ScalarAdapterCtx adapter;
    real_t state[1];
    SolverStepResult result;

    state[0] = *y;
    adapter.rhs = rhs;
    adapter.rhs_ctx = rhs_ctx;

    result = rkf78_core_run(
        RKF78_MODE_STEP,
        config->plugins,
        state,
        1,
        config->t0,
        config->t0 + config->h,
        config->h,
        config->tol,
        rkf78_scalar_rhs_adapter,
        &adapter);

    *y = state[0];
    return result;
}

void rkf78_backend_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    (void)rkf78_core_run(
        RKF78_MODE_RUN,
        config->plugins,
        state,
        n,
        config->t0,
        config->t_end,
        config->h_init,
        config->tol,
        rhs,
        rhs_ctx);
}

SolverStepResult rkf78_backend_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    return rkf78_core_run(
        RKF78_MODE_STEP,
        config->plugins,
        state,
        n,
        config->t0,
        config->t0 + config->h,
        config->h,
        config->tol,
        rhs,
        rhs_ctx);
}
