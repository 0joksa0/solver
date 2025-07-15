#define SOLVER_REAL_AS long double

#include <criterion/criterion.h>
#include <solver.h>

static real_t rhs_poly(real_t t, real_t y, void* unused)
{
    (void)y;
    (void)unused;
    return t * t;
}

Test(Solver, Polynomial_T3over3_ODE45)
{
    real_t y0 = REAL(0.0);
    real_t y_end = solve(SOLVER_ODE45, REAL(1e-8),
        y0, REAL(0.0), REAL(2.0),
        REAL(0.2), rhs_poly, NULL);

    real_t expected = REAL((2.0 * 2.0 * 2.0) / 3.0);
    cr_assert(RABS(y_end - expected) < REAL(1e-6),
        "Expected %.10f, got %.10f\n",
        (double)expected, (double)y_end);
}

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
}

Test(Solver, SineSolution_ODE45)
{
    real_t y_end = solve(SOLVER_ODE45, REAL(1e-9),
        REAL(0.0), REAL(0.0), REAL(3.14159265358979323846 / 2),
        REAL(0.1), rhs_cosine, NULL);

    real_t expected = REAL(1.0);
    cr_assert(RABS(y_end - expected) < REAL(1e-7),
        "Expected %.10f, got %.10f\n",
        (double)expected, (double)y_end);
}

static real_t rhs_logistic(real_t t, real_t y, void* p)
{
    (void)t;
    real_t r = ((real_t*)p)[0];
    real_t K = ((real_t*)p)[1];
    return r * y * (REAL(1.0) - y / K);
}

Test(Solver, Logistic_ODE45)
{
    real_t p[2] = { REAL(1.0), REAL(10.0) };
    real_t y0 = REAL(1.0);
    real_t t_end = REAL(5.0);

    real_t y_end = solve(SOLVER_ODE45, REAL(1e-8),
        y0, REAL(0.0), t_end, REAL(0.2),
        rhs_logistic, p);

    real_t K = p[1], r = p[0];
    real_t expected = K / (REAL(1.0) + ((K - y0) / y0) * REXP(-r * t_end));

    cr_assert(RABS(y_end - expected) < REAL(1e-6),
        "Expected %.10f, got %.10f\n",
        (double)expected, (double)y_end);
}

static real_t rhs_decay(real_t t, real_t y, void* k_ptr)
{
    (void)t;
    real_t k = *(real_t*)k_ptr;
    return k * y;
}

Test(Solver, ExponentialDecay_Short_ODE45)
{
    real_t k = REAL(-0.5);
    real_t y0 = REAL(1.0);
    real_t y_end = solve(SOLVER_ODE45, REAL(1e-8),
        y0, REAL(0.0), REAL(2.0),
        REAL(0.1), rhs_decay, &k);

    real_t expected = y0 * REXP(k * REAL(2.0));
    cr_assert(RABS(y_end - expected) < REAL(1e-6),
        "Expected %.10f, got %.10f\n",
        (double)expected, (double)y_end);
}
