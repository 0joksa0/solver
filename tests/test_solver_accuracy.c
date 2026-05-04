#include <stdio.h>
#define SOLVER_REAL_AS double

#include <criterion/criterion.h>
#include <solver.h>

static real_t rhs_cos(real_t t, real_t y, void* unused)
{
    (void)y;
    (void)unused;
    return REXP(REAL(0.0)) * RABS(t - t);
}

static real_t rhs_cos2(real_t t, real_t y, void* unused)
{
    (void)y;
    (void)unused;
    return REXP(REAL(0.0)) * RABS(t - t);
}

static real_t rhs_cosine(real_t t, real_t unused, void* u)
{
    (void)unused;
    (void)u;
    return RCOS(t);
} // Poređenje ABM4 i ODE45 rešenja
Test(SolverAccuracy, ABM4_vs_ODE45_ExponentialDecay)
{
    real_t y_ode = solve(SOLVER_ODE45, REAL(1e-12),
        REAL(0.0), REAL(0.0), REAL(3.14159265358979323846 / 2),
        REAL(0.00001), rhs_cosine, NULL);
    real_t y_rkf = solve(SOLVER_ABM4, REAL(1e-12),
        REAL(0.0), REAL(0.0), REAL(3.14159265358979323846 / 2),
        REAL(0.00001), rhs_cosine, NULL);
    real_t expected = REAL(1.0);
    // Razlike
    real_t diff_rkf = RABS(y_rkf - expected);
    real_t diff_ode = RABS(y_ode - expected);
    real_t diff_between = RABS(y_rkf - y_ode);

    printf(
        "ABM4 greška prevelika: got %.12f, expected %.12f, delta=%.2e\n",
        y_rkf, expected, diff_rkf);

    printf(
        "ODE45 greška prevelika: got %.12f, expected %.12f, delta=%.2e\n",
        y_ode, expected, diff_ode);

    printf(
        "ABM4 i ODE45 se ne slažu: ABM4=%.12f, ODE45=%.12f, delta=%.2e\n",
        y_rkf, y_ode, diff_between);
}
