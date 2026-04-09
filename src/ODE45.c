#include "ODE45.h"
#include "logger.h"
#include "observers.h"
#include <stdlib.h>
#include <string.h>

real_t ode45_solver(
    real_t y0, real_t t0, real_t t_end,
    real_t h_init, real_t tol,
    ODEFunction f, void* params)
{
    real_t t = t0;
    real_t y = y0;
    real_t h = h_init;

    const real_t safety = REAL(0.9);
    const real_t min_scale = REAL(0.1);
    const real_t max_scale = REAL(5.0);

    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }

        real_t k1 = h * f(t, y, params);
        real_t k2 = h * f(t + h * (REAL(1.0) / REAL(5.0)), y + k1 * (REAL(1.0) / REAL(5.0)), params);
        real_t k3 = h * f(t + h * (REAL(3.0) / REAL(10.0)), y + k1 * (REAL(3.0) / REAL(40.0)) + k2 * (REAL(9.0) / REAL(40.0)), params);
        real_t k4 = h * f(t + h * (REAL(4.0) / REAL(5.0)), y + k1 * (REAL(44.0) / REAL(45.0)) - k2 * (REAL(56.0) / REAL(15.0)) + k3 * (REAL(32.0) / REAL(9.0)), params);
        real_t k5 = h * f(t + h * (REAL(8.0) / REAL(9.0)), y + k1 * (REAL(19372.0) / REAL(6561.0)) - k2 * (REAL(25360.0) / REAL(2187.0)) + k3 * (REAL(64448.0) / REAL(6561.0)) - k4 * (REAL(212.0) / REAL(729.0)), params);
        real_t k6 = h * f(t + h, y + k1 * (REAL(9017.0) / REAL(3168.0)) - k2 * (REAL(355.0) / REAL(33.0)) + k3 * (REAL(46732.0) / REAL(5247.0)) + k4 * (REAL(49.0) / REAL(176.0)) - k5 * (REAL(5103.0) / REAL(18656.0)), params);
        real_t k7 = h * f(t + h, y + k1 * (REAL(35.0) / REAL(384.0)) + REAL(0.0) + k3 * (REAL(500.0) / REAL(1113.0)) + k4 * (REAL(125.0) / REAL(192.0)) - k5 * (REAL(2187.0) / REAL(6784.0)) + k6 * (REAL(11.0) / REAL(84.0)), params);

        real_t y5 = y + k1 * (REAL(35.0) / REAL(384.0)) + k3 * (REAL(500.0) / REAL(1113.0))
            + k4 * (REAL(125.0) / REAL(192.0)) - k5 * (REAL(2187.0) / REAL(6784.0))
            + k6 * (REAL(11.0) / REAL(84.0));

        real_t y4 = y + k1 * (REAL(5179.0) / REAL(57600.0)) + k3 * (REAL(7571.0) / REAL(16695.0))
            + k4 * (REAL(393.0) / REAL(640.0)) - k5 * (REAL(92097.0) / REAL(339200.0))
            + k6 * (REAL(187.0) / REAL(2100.0)) + k7 * (REAL(1.0) / REAL(40.0));

        real_t err = RABS(y5 - y4);

        if (err <= tol || h < REAL(1e-6)) {
            t += h;
            y = y5;
        } else {
        }

        real_t scale = safety * RPOW(tol / (err + REAL_EPSILON), REAL(0.2));
        scale = RMAX(min_scale, RMIN(scale, max_scale));
        h *= scale;
    }

    return y;
}

void ode45_vector_step(
    real_t* y,
    size_t n,
    real_t* t,
    real_t* h,
    real_t tol,
    RHS f,
    void* ctx)
{
    const real_t safety = REAL(0.9);
    const real_t min_scale = REAL(0.1);
    const real_t max_scale = REAL(5.0);

    real_t* k1 = calloc(n, sizeof(real_t));
    real_t* k2 = calloc(n, sizeof(real_t));
    real_t* k3 = calloc(n, sizeof(real_t));
    real_t* k4 = calloc(n, sizeof(real_t));
    real_t* k5 = calloc(n, sizeof(real_t));
    real_t* k6 = calloc(n, sizeof(real_t));
    real_t* k7 = calloc(n, sizeof(real_t));

    real_t* yt = calloc(n, sizeof(real_t));
    real_t* y5 = calloc(n, sizeof(real_t));
    real_t* y4 = calloc(n, sizeof(real_t));

    while (1) {
        /* k1 */
        f(*t, y, k1, ctx);
        for (size_t i = 0; i < n; ++i)
            k1[i] *= *h;

        /* k2 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = y[i] + k1[i] * (REAL(1.0) / REAL(5.0));
        f(*t + *h * (REAL(1.0) / REAL(5.0)), yt, k2, ctx);
        for (size_t i = 0; i < n; ++i)
            k2[i] *= *h;

        /* k3 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = y[i]
                + k1[i] * (REAL(3.0) / REAL(40.0))
                + k2[i] * (REAL(9.0) / REAL(40.0));
        f(*t + *h * (REAL(3.0) / REAL(10.0)), yt, k3, ctx);
        for (size_t i = 0; i < n; ++i)
            k3[i] *= *h;

        /* k4 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = y[i]
                + k1[i] * (REAL(44.0) / REAL(45.0))
                - k2[i] * (REAL(56.0) / REAL(15.0))
                + k3[i] * (REAL(32.0) / REAL(9.0));
        f(*t + *h * (REAL(4.0) / REAL(5.0)), yt, k4, ctx);
        for (size_t i = 0; i < n; ++i)
            k4[i] *= *h;

        /* k5 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = y[i]
                + k1[i] * (REAL(19372.0) / REAL(6561.0))
                - k2[i] * (REAL(25360.0) / REAL(2187.0))
                + k3[i] * (REAL(64448.0) / REAL(6561.0))
                - k4[i] * (REAL(212.0) / REAL(729.0));
        f(*t + *h * (REAL(8.0) / REAL(9.0)), yt, k5, ctx);
        for (size_t i = 0; i < n; ++i)
            k5[i] *= *h;

        /* k6 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = y[i]
                + k1[i] * (REAL(9017.0) / REAL(3168.0))
                - k2[i] * (REAL(355.0) / REAL(33.0))
                + k3[i] * (REAL(46732.0) / REAL(5247.0))
                + k4[i] * (REAL(49.0) / REAL(176.0))
                - k5[i] * (REAL(5103.0) / REAL(18656.0));
        f(*t + *h, yt, k6, ctx);
        for (size_t i = 0; i < n; ++i)
            k6[i] *= *h;

        /* k7 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = y[i]
                + k1[i] * (REAL(35.0) / REAL(384.0))
                + k3[i] * (REAL(500.0) / REAL(1113.0))
                + k4[i] * (REAL(125.0) / REAL(192.0))
                - k5[i] * (REAL(2187.0) / REAL(6784.0))
                + k6[i] * (REAL(11.0) / REAL(84.0));
        f(*t + *h, yt, k7, ctx);
        for (size_t i = 0; i < n; ++i)
            k7[i] *= *h;

        /* y5 i y4 */
        real_t err = REAL(0.0);
        for (size_t i = 0; i < n; ++i) {
            y5[i] = y[i]
                + k1[i] * (REAL(35.0) / REAL(384.0))
                + k3[i] * (REAL(500.0) / REAL(1113.0))
                + k4[i] * (REAL(125.0) / REAL(192.0))
                - k5[i] * (REAL(2187.0) / REAL(6784.0))
                + k6[i] * (REAL(11.0) / REAL(84.0));

            y4[i] = y[i]
                + k1[i] * (REAL(5179.0) / REAL(57600.0))
                + k3[i] * (REAL(7571.0) / REAL(16695.0))
                + k4[i] * (REAL(393.0) / REAL(640.0))
                - k5[i] * (REAL(92097.0) / REAL(339200.0))
                + k6[i] * (REAL(187.0) / REAL(2100.0))
                + k7[i] * (REAL(1.0) / REAL(40.0));

            err = RMAX(err, RABS(y5[i] - y4[i]));
        }

        if (err <= tol || *h < REAL(1e-6)) {
            memcpy(y, y5, n * sizeof(real_t));
            *t += *h;
            break;
        }

        real_t scale = safety * RPOW(tol / (err + REAL_EPSILON), REAL(0.2));
        scale = RMAX(min_scale, RMIN(scale, max_scale));
        *h *= scale;
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
}

void ode45_vector_solve(
    real_t* x,
    size_t n,
    real_t t0,
    real_t t_end,
    real_t h_init,
    real_t rtol,
    real_t atol,
    RHS rhs,
    void* rhs_ctx,
    VectorObserver obs,
    void* obs_ctx,
    SolverStats* stats,
    SolverLogger* logger)
{
    real_t t = t0;
    real_t h = h_init;

    if (stats) {
        solver_stats_init(stats);
        solver_stats_start_timer(stats);
    }

    const real_t safety = REAL(0.9);
    const real_t min_scale = REAL(0.001);
    const real_t max_scale = REAL(1.5);

    /* PI controller gains for RK45 (order 5) */
    const real_t kP = REAL(0.2); /* 1/(p+1) with p=4 error order -> 1/5 */
    const real_t kI = REAL(0.08); /* typical: 0.4/(p+1) = 0.4/5 */

    /* keep previous error (start with 1 so it doesn't blow up) */
    real_t err_prev = REAL(1.0);

    /* radni baferi */
    real_t* k1 = calloc(n, sizeof(real_t));
    real_t* k2 = calloc(n, sizeof(real_t));
    real_t* k3 = calloc(n, sizeof(real_t));
    real_t* k4 = calloc(n, sizeof(real_t));
    real_t* k5 = calloc(n, sizeof(real_t));
    real_t* k6 = calloc(n, sizeof(real_t));
    real_t* k7 = calloc(n, sizeof(real_t));

    real_t* yt = calloc(n, sizeof(real_t));
    real_t* y5 = calloc(n, sizeof(real_t));
    real_t* y4 = calloc(n, sizeof(real_t));

    if (obs)
        obs(t, x, n, obs_ctx);

    while (t < t_end) {

        if (t + h > t_end)
            h = t_end - t;

        /* k1 */
        rhs(t, x, k1, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k1[i] *= h;

        /* k2 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = x[i] + k1[i] * (REAL(1.0) / REAL(5.0));
        rhs(t + h * (REAL(1.0) / REAL(5.0)), yt, k2, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k2[i] *= h;

        /* k3 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = x[i]
                + k1[i] * (REAL(3.0) / REAL(40.0))
                + k2[i] * (REAL(9.0) / REAL(40.0));
        rhs(t + h * (REAL(3.0) / REAL(10.0)), yt, k3, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k3[i] *= h;

        /* k4 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = x[i]
                + k1[i] * (REAL(44.0) / REAL(45.0))
                - k2[i] * (REAL(56.0) / REAL(15.0))
                + k3[i] * (REAL(32.0) / REAL(9.0));
        rhs(t + h * (REAL(4.0) / REAL(5.0)), yt, k4, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k4[i] *= h;

        /* k5 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = x[i]
                + k1[i] * (REAL(19372.0) / REAL(6561.0))
                - k2[i] * (REAL(25360.0) / REAL(2187.0))
                + k3[i] * (REAL(64448.0) / REAL(6561.0))
                - k4[i] * (REAL(212.0) / REAL(729.0));
        rhs(t + h * (REAL(8.0) / REAL(9.0)), yt, k5, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k5[i] *= h;

        /* k6 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = x[i]
                + k1[i] * (REAL(9017.0) / REAL(3168.0))
                - k2[i] * (REAL(355.0) / REAL(33.0))
                + k3[i] * (REAL(46732.0) / REAL(5247.0))
                + k4[i] * (REAL(49.0) / REAL(176.0))
                - k5[i] * (REAL(5103.0) / REAL(18656.0));
        rhs(t + h, yt, k6, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k6[i] *= h;

        /* k7 */
        for (size_t i = 0; i < n; ++i)
            yt[i] = x[i]
                + k1[i] * (REAL(35.0) / REAL(384.0))
                + k3[i] * (REAL(500.0) / REAL(1113.0))
                + k4[i] * (REAL(125.0) / REAL(192.0))
                - k5[i] * (REAL(2187.0) / REAL(6784.0))
                + k6[i] * (REAL(11.0) / REAL(84.0));
        rhs(t + h, yt, k7, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k7[i] *= h;

        /* error + accept/reject */
        real_t err_sum = REAL(0.0);
        for (size_t i = 0; i < n; ++i) {
            y5[i] = x[i]
                + k1[i] * (REAL(35.0) / REAL(384.0))
                + k3[i] * (REAL(500.0) / REAL(1113.0))
                + k4[i] * (REAL(125.0) / REAL(192.0))
                - k5[i] * (REAL(2187.0) / REAL(6784.0))
                + k6[i] * (REAL(11.0) / REAL(84.0));

            y4[i] = x[i]
                + k1[i] * (REAL(5179.0) / REAL(57600.0))
                + k3[i] * (REAL(7571.0) / REAL(16695.0))
                + k4[i] * (REAL(393.0) / REAL(640.0))
                - k5[i] * (REAL(92097.0) / REAL(339200.0))
                + k6[i] * (REAL(187.0) / REAL(2100.0))
                + k7[i] * (REAL(1.0) / REAL(40.0));

            real_t scale = atol + rtol * RMAX(RABS(x[i]), RABS(y5[i]));
            real_t e = (y5[i] - y4[i]) / scale;

            err_sum += e * e;
        }
        real_t err = RSQRT(err_sum / REAL(n));

        int accepted = (err <= REAL(1.0));

        real_t dt_old = h;

        /* PI step-size controller */
        real_t err_safe = RMAX(err, REAL_EPSILON);
        real_t err_prev_safe = RMAX(err_prev, REAL_EPSILON);

        /* scale = safety * err^{-kP} * err_prev^{kI} */
        real_t scale = safety
            * RPOW(err_safe, -kP)
            * RPOW(err_prev_safe, kI);

        /* clamp */
        scale = RMAX(min_scale, RMIN(scale, max_scale));

        real_t dt_new = h * scale;

        if (stats) {
            StepInfo info;
            info.err_norm = err;
            info.dt_old = dt_old;
            info.dt_new = dt_new;
            info.accepted = accepted;
            solver_stats_update(stats, &info);
            if (logger)
                solver_logger_log_step(logger, t, &info);
        }

        if (accepted) {
            memcpy(x, y5, n * sizeof(real_t));
            t += h;

            /* update PI memory on accepted step */
            err_prev = err_safe;

            if (obs)
                obs(t, x, n, obs_ctx);
        }
        if (!accepted) {
            /* extra shrink on reject */
            h = h * RMAX(min_scale, safety * REAL(0.5));
            continue;
        }

        h = dt_new;
    }

    if (stats) {
        solver_stats_stop_timer(stats);
        solver_stats_finalize(stats);
        if (logger)
            solver_logger_log_summary(logger, stats);
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
}
