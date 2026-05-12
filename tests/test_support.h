#ifndef TEST_SUPPORT_H
#define TEST_SUPPORT_H

#include <criterion/criterion.h>
#include <framework/default_runtime.h>
#include <logger.h>
#include <solver.h>

#include "models/harmonic_oscillator.h"

typedef struct ObserverProbe {
    size_t calls;
    real_t last_t;
} ObserverProbe;

typedef struct PluginTrace {
    int markers[16];
    size_t count;
} PluginTrace;

typedef struct TracePluginCtx {
    PluginTrace* trace;
    int marker;
    SolverPluginResult result;
} TracePluginCtx;

/* @md
## ODE
$y'(t) = t^2$.
## Explicit solution
$y(t) = \frac{t^3}{3} + C$.
For $y(0) = 0$:
$y(t) = \frac{t^3}{3}$
*/
static real_t rhs_poly(real_t t, real_t y, void* unused)
{
    (void)y;
    (void)unused;
    return t * t;
}

/* @md
## ODE
$y'(t) = k y(t)$.
## Explicit solution
$y(t) = y_0 \cdot e^{k (t - t_0)}$.
*/
static real_t rhs_decay_scalar(real_t t, real_t y, void* k_ptr)
{
    real_t k = *(real_t*)k_ptr;
    (void)t;
    return k * y;
}

/* @md
## Vector form
$x'(t) = k x(t)$.
## Explicit solution
$x(t) = x_0 \cdot e^{k (t - t_0)}$.
*/
static void rhs_decay_vector(real_t t, const real_t* state, real_t* dstate, size_t n, void* k_ptr)
{
    real_t k = *(real_t*)k_ptr;
    (void)t;
    (void)n;
    dstate[0] = k * state[0];
}

static void probe_observer(real_t t, const real_t* x, size_t n, void* ctx)
{
    ObserverProbe* probe = (ObserverProbe*)ctx;
    (void)x;
    (void)n;
    probe->calls++;
    probe->last_t = t;
}

static SolverPluginResult trace_on_step_accepted(void* plugin_ctx, const SolverEvent* event)
{
    TracePluginCtx* ctx = (TracePluginCtx*)plugin_ctx;
    (void)event;
    ctx->trace->markers[ctx->trace->count++] = ctx->marker;
    return ctx->result;
}

static SolverPluginResult trace_on_run_finish(void* plugin_ctx, const SolverEvent* event)
{
    TracePluginCtx* ctx = (TracePluginCtx*)plugin_ctx;
    (void)event;
    ctx->trace->markers[ctx->trace->count++] = ctx->marker;
    return ctx->result;
}

#endif
