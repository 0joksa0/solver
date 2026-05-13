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
} 
