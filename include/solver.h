#ifndef SOLVER_H
#define SOLVER_H

#ifndef SOLVER_REAL_AS
#define SOLVER_REAL_AS float
#endif

typedef SOLVER_REAL_AS real_t;

#include <float.h>
#include <math.h>

#if defined(SOLVER_REAL_IS_DOUBLE)
#define REAL_EPSILON DBL_EPSILON
#define RSQRT sqrt
#define RABS fabs
#define RMAX fmax
#define RMIN fmin
#define RPOW pow
#define REXP exp
#define RLOG log
#define RSIN sin
#define RCOS cos
#define RTAN tan
#define RASIN asin
#define RACOS acos
#define RATAN atan
#define RSINH sinh
#define RCOSH cosh
#define RTANH tanh
#define REXP2 exp2
#define RPOW10 exp10

#elif defined(SOLVER_REAL_IS_LONG_DOUBLE)
#define REAL_EPSILON LDBL_EPSILON
#define RSQRT sqrtl
#define RABS fabsl
#define RMAX fmaxl
#define RMIN fminl
#define RPOW powl
#define REXP expl
#define RLOG logl
#define RSIN sinl
#define RCOS cosl
#define RTAN tanl
#define RASIN asinl
#define RACOS acosl
#define RATAN atanl
#define RSINH sinhl
#define RCOSH coshl
#define RTANH tanhl
#define REXP2 exp2l
#define RPOW10 exp10l

#else
#define REAL_EPSILON FLT_EPSILON
#define RSQRT sqrtf
#define RABS fabsf
#define RMAX fmaxf
#define RMIN fminf
#define RPOW powf
#define REXP expf
#define RLOG logf
#define RSIN sinf
#define RCOS cosf
#define RTAN tanf
#define RASIN asinf
#define RACOS acosf
#define RATAN atanf
#define RSINH sinhf
#define RCOSH coshf
#define RTANH tanhf
#define REXP2 exp2f
#define RPOW10 exp10f
#endif

#define REAL(x) ((real_t)(x))

#include <stddef.h>

typedef enum {
    SOLVER_RK4,
    SOLVER_ODE45,
    SOLVER_RKF78
} SolverType;

typedef real_t (*ODEFunction)(real_t t, real_t y, void* params);

real_t solve(SolverType type, real_t tol, real_t y0, real_t t0, real_t t_end, real_t h_init, ODEFunction f, void* params);

real_t solve_step(SolverType type, real_t tol, real_t y_n, real_t t_n, real_t h, ODEFunction f, void* params);

#endif /* SOLVER_H */
