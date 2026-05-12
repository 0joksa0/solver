#include "backends/ode45_backend.h"

#include "framework/solver_runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct ODE45ScalarAdapterCtx {
    ODEFunction rhs;
    void* rhs_ctx;
} ODE45ScalarAdapterCtx;

typedef enum ODE45Mode {
    ODE45_MODE_RUN = 0,
    ODE45_MODE_STEP = 1
} ODE45Mode;

static void ode45_scalar_rhs_adapter(
    real_t t,
    const real_t* state,
    real_t* dstate,
    size_t n,
    void* ctx)
{
    ODE45ScalarAdapterCtx* adapter = (ODE45ScalarAdapterCtx*)ctx;
    (void)n;

    dstate[0] = adapter->rhs(t, state[0], adapter->rhs_ctx);
}

static void* ode45_xcalloc(size_t count, size_t size)
{
    void* ptr = calloc(count, size);

    if (!ptr) {
        fprintf(stderr, "[ode45_backend] allocation failed\n");
        exit(EXIT_FAILURE);
    }

    return ptr;
}

static void ode45_map_tol(real_t tol, real_t* rtol, real_t* atol)
{
    real_t safe_tol = RMAX(tol, REAL_EPSILON);

    *rtol = safe_tol;
    *atol = safe_tol;
}

static SolverStepResult ode45_core_run(
    ODE45Mode mode,
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
    const real_t kP = REAL(0.2);
    const real_t kI = REAL(0.08);
    real_t err_prev = REAL(1.0);
    real_t* k1;
    real_t* k2;
    real_t* k3;
    real_t* k4;
    real_t* k5;
    real_t* k6;
    real_t* k7;
    real_t* yt;
    real_t* y5;
    real_t* y4;

    result.t_out = t0;
    result.h_used = REAL(0.0);
    result.h_next = h_init;
    result.accepted = 0;
    result.stopped_by_plugin = 0;

    if (n == 0) {
        return result;
    }

    ode45_map_tol(tol, &rtol, &atol);
    solver_runtime_init(&runtime, (SolverPluginManager*)plugins);

    k1 = ode45_xcalloc(n, sizeof(real_t));
    k2 = ode45_xcalloc(n, sizeof(real_t));
    k3 = ode45_xcalloc(n, sizeof(real_t));
    k4 = ode45_xcalloc(n, sizeof(real_t));
    k5 = ode45_xcalloc(n, sizeof(real_t));
    k6 = ode45_xcalloc(n, sizeof(real_t));
    k7 = ode45_xcalloc(n, sizeof(real_t));
    yt = ode45_xcalloc(n, sizeof(real_t));
    y5 = ode45_xcalloc(n, sizeof(real_t));
    y4 = ode45_xcalloc(n, sizeof(real_t));

    if (solver_runtime_emit_run_start(&runtime, t, state, n) == SOLVER_PLUGIN_STOP) {
        goto finish;
    }

    while (mode == ODE45_MODE_STEP || t < t_end) {
        real_t err_sum = REAL(0.0);
        real_t step_t0 = t;
        real_t step_h = h;
        real_t err;
        real_t err_safe;
        real_t err_prev_safe;
        real_t scale;
        real_t next_h;
        StepInfo info;
        size_t i;

        if (mode == ODE45_MODE_RUN && t + step_h > t_end) {
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
            k1[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            yt[i] = state[i] + k1[i] * (REAL(1.0) / REAL(5.0));
        }
        rhs(t + step_h * (REAL(1.0) / REAL(5.0)), yt, k2, n, rhs_ctx);
        for (i = 0; i < n; ++i) {
            k2[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            yt[i] = state[i]
                + k1[i] * (REAL(3.0) / REAL(40.0))
                + k2[i] * (REAL(9.0) / REAL(40.0));
        }
        rhs(t + step_h * (REAL(3.0) / REAL(10.0)), yt, k3, n, rhs_ctx);
        for (i = 0; i < n; ++i) {
            k3[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            yt[i] = state[i]
                + k1[i] * (REAL(44.0) / REAL(45.0))
                - k2[i] * (REAL(56.0) / REAL(15.0))
                + k3[i] * (REAL(32.0) / REAL(9.0));
        }
        rhs(t + step_h * (REAL(4.0) / REAL(5.0)), yt, k4, n, rhs_ctx);
        for (i = 0; i < n; ++i) {
            k4[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            yt[i] = state[i]
                + k1[i] * (REAL(19372.0) / REAL(6561.0))
                - k2[i] * (REAL(25360.0) / REAL(2187.0))
                + k3[i] * (REAL(64448.0) / REAL(6561.0))
                - k4[i] * (REAL(212.0) / REAL(729.0));
        }
        rhs(t + step_h * (REAL(8.0) / REAL(9.0)), yt, k5, n, rhs_ctx);
        for (i = 0; i < n; ++i) {
            k5[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            yt[i] = state[i]
                + k1[i] * (REAL(9017.0) / REAL(3168.0))
                - k2[i] * (REAL(355.0) / REAL(33.0))
                + k3[i] * (REAL(46732.0) / REAL(5247.0))
                + k4[i] * (REAL(49.0) / REAL(176.0))
                - k5[i] * (REAL(5103.0) / REAL(18656.0));
        }
        rhs(t + step_h, yt, k6, n, rhs_ctx);
        for (i = 0; i < n; ++i) {
            k6[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            yt[i] = state[i]
                + k1[i] * (REAL(35.0) / REAL(384.0))
                + k3[i] * (REAL(500.0) / REAL(1113.0))
                + k4[i] * (REAL(125.0) / REAL(192.0))
                - k5[i] * (REAL(2187.0) / REAL(6784.0))
                + k6[i] * (REAL(11.0) / REAL(84.0));
        }
        rhs(t + step_h, yt, k7, n, rhs_ctx);
        for (i = 0; i < n; ++i) {
            k7[i] *= step_h;
        }

        for (i = 0; i < n; ++i) {
            real_t scale_i;
            real_t e;

            y5[i] = state[i]
                + k1[i] * (REAL(35.0) / REAL(384.0))
                + k3[i] * (REAL(500.0) / REAL(1113.0))
                + k4[i] * (REAL(125.0) / REAL(192.0))
                - k5[i] * (REAL(2187.0) / REAL(6784.0))
                + k6[i] * (REAL(11.0) / REAL(84.0));

            y4[i] = state[i]
                + k1[i] * (REAL(5179.0) / REAL(57600.0))
                + k3[i] * (REAL(7571.0) / REAL(16695.0))
                + k4[i] * (REAL(393.0) / REAL(640.0))
                - k5[i] * (REAL(92097.0) / REAL(339200.0))
                + k6[i] * (REAL(187.0) / REAL(2100.0))
                + k7[i] * (REAL(1.0) / REAL(40.0));

            scale_i = atol + rtol * RMAX(RABS(state[i]), RABS(y5[i]));
            e = (y5[i] - y4[i]) / scale_i;
            err_sum += e * e;
        }

        err = RSQRT(err_sum / REAL(n));
        err_safe = RMAX(err, REAL_EPSILON);
        err_prev_safe = RMAX(err_prev, REAL_EPSILON);
        scale = safety * RPOW(err_safe, -kP) * RPOW(err_prev_safe, kI);
        scale = RMAX(min_scale, RMIN(scale, max_scale));
        next_h = step_h * scale;

        info.err_norm = err;
        info.dt_old = step_h;
        info.dt_new = next_h;
        info.accepted = (err <= REAL(1.0));

        h = next_h;
        result.h_next = h;

        if (info.accepted) {
            memcpy(state, y5, n * sizeof(real_t));
            t += step_h;
            err_prev = err_safe;

            result.t_out = t;
            result.h_used = step_h;
            result.accepted = 1;

            if (solver_runtime_emit_step_accepted(&runtime, step_t0, t, state, n, &info, 7) == SOLVER_PLUGIN_STOP) {
                break;
            }

            if (mode == ODE45_MODE_STEP) {
                break;
            }

            continue;
        }

        result.t_out = t;
        result.h_used = REAL(0.0);
        if (solver_runtime_emit_step_rejected(&runtime, step_t0, step_t0 + step_h, state, n, &info, 7) == SOLVER_PLUGIN_STOP) {
            break;
        }
    }

finish:
    solver_runtime_emit_run_finish(&runtime, t, state, n);
    result.stopped_by_plugin = runtime.stopped_by_plugin;
    if (result.accepted && result.t_out == t0) {
        result.t_out = t;
    }
    if (!result.accepted && result.h_next == REAL(0.0)) {
        result.h_next = h;
    }

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(k5);
    free(k6);
    free(k7);
    free(yt);
    free(y5);
    free(y4);

    return result;
}

real_t ode45_backend_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx)
{
    ODE45ScalarAdapterCtx adapter;
    real_t state[1];

    state[0] = y0;
    adapter.rhs = rhs;
    adapter.rhs_ctx = rhs_ctx;

    ode45_core_run(
        ODE45_MODE_RUN,
        config->plugins,
        state,
        1,
        config->t0,
        config->t_end,
        config->h_init,
        config->tol,
        ode45_scalar_rhs_adapter,
        &adapter);

    return state[0];
}

SolverStepResult ode45_backend_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx)
{
    ODE45ScalarAdapterCtx adapter;
    real_t state[1];
    SolverStepResult result;

    state[0] = *y;
    adapter.rhs = rhs;
    adapter.rhs_ctx = rhs_ctx;

    result = ode45_core_run(
        ODE45_MODE_STEP,
        config->plugins,
        state,
        1,
        config->t0,
        config->t0 + config->h,
        config->h,
        config->tol,
        ode45_scalar_rhs_adapter,
        &adapter);

    *y = state[0];
    return result;
}

void ode45_backend_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    (void)ode45_core_run(
        ODE45_MODE_RUN,
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

SolverStepResult ode45_backend_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx)
{
    return ode45_core_run(
        ODE45_MODE_STEP,
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
