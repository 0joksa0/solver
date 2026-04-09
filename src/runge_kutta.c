#include "runge_cutta.h"
#include <stdio.h>
#include <stdlib.h>

real_t runge_kutta_4(real_t y_n, real_t t_n, real_t h, ODEFunction f, void* params)
{
    real_t k1 = h * f(t_n, y_n, params);
    real_t k2 = h * f(t_n + (h * REAL(0.5)), y_n + (k1 * REAL(0.5)), params);
    real_t k3 = h * f(t_n + (h * REAL(0.5)), y_n + (k2 * REAL(0.5)), params);
    real_t k4 = h * f(t_n + h, y_n + k3, params);

    return y_n + ((REAL(1.0) / REAL(6.0)) * (k1 + (REAL(2.0) * k2) + (REAL(2.0) * k3) + k4));
}

void vector_runge_kutta_4_step(
    real_t t,
    real_t h,
    void* state,
    size_t state_size,
    RHS rhs,
    void* user_ctx)
{

    size_t n = state_size / sizeof(real_t);
    real_t* y = (real_t*)state;

    real_t* k1 = calloc(n, sizeof(real_t));
    real_t* k2 = calloc(n, sizeof(real_t));
    real_t* k3 = calloc(n, sizeof(real_t));
    real_t* k4 = calloc(n, sizeof(real_t));
    real_t* yt = calloc(n, sizeof(real_t));

    // k1
    rhs(t, y, k1, user_ctx);

    // k2
    for (size_t i = 0; i < n; i++)
        yt[i] = y[i] + h * 0.5 * k1[i];
    rhs(t + h * 0.5, yt, k2, user_ctx);

    // k3
    for (size_t i = 0; i < n; i++)
        yt[i] = y[i] + h * 0.5 * k2[i];
    rhs(t + h * 0.5, yt, k3, user_ctx);

    // k4
    for (size_t i = 0; i < n; i++)
        yt[i] = y[i] + h * k3[i];
    rhs(t + h, yt, k4, user_ctx);

    // update
    for (size_t i = 0; i < n; i++)
        y[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);

    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(yt);
}

void rk4_vector_solve(
    real_t* x,
    size_t n,
    real_t t0,
    real_t t_end,
    real_t h,
    RHS rhs,
    void* rhs_ctx,
    VectorObserver obs,
    void* obs_ctx)
{
    real_t t = t0;

    real_t* k1 = calloc(n, sizeof(real_t));
    real_t* k2 = calloc(n, sizeof(real_t));
    real_t* k3 = calloc(n, sizeof(real_t));
    real_t* k4 = calloc(n, sizeof(real_t));
    real_t* xt = calloc(n, sizeof(real_t));

    if (!k1 || !k2 || !k3 || !k4 || !xt) {
        fprintf(stderr, "[rk4_vector_solve] allocation failed\n");
        exit(EXIT_FAILURE);
    }

    if (obs) {
        obs(t, x, n, obs_ctx);
    }

    while (t < t_end) {
        if (t + h > t_end)
            h = t_end - t;

        /* k1 */
        rhs(t, x, k1, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k1[i] *= h;

        /* k2 */
        for (size_t i = 0; i < n; ++i)
            xt[i] = x[i] + REAL(0.5) * k1[i];
        rhs(t + REAL(0.5) * h, xt, k2, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k2[i] *= h;

        /* k3 */
        for (size_t i = 0; i < n; ++i)
            xt[i] = x[i] + REAL(0.5) * k2[i];
        rhs(t + REAL(0.5) * h, xt, k3, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k3[i] *= h;

        /* k4 */
        for (size_t i = 0; i < n; ++i)
            xt[i] = x[i] + k3[i];
        rhs(t + h, xt, k4, rhs_ctx);
        for (size_t i = 0; i < n; ++i)
            k4[i] *= h;

        /* update state */
        for (size_t i = 0; i < n; ++i) {
            x[i] += (k1[i]
                        + REAL(2.0) * k2[i]
                        + REAL(2.0) * k3[i]
                        + k4[i])
                / REAL(6.0);
        }

        t += h;
        if (obs) {
            obs(t, x, n, obs_ctx);
        }
        /* observer (ako postoji) */
    }
    free(k1);
    free(k2);
    free(k3);
    free(k4);
    free(xt);
}
