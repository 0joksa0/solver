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
#define REAL_MAX DBL_MAX
#define REAL_MIN DBL_MIN
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
#define REAL_MAX LDBL_MAX
#define REAL_MIN LDBL_MIN
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
#define REAL_MAX FLT_MAX
#define REAL_MIN FLT_MIN
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

typedef real_t (*ODEFunction)(
    real_t t,
    real_t y,
    void* params);

typedef void (*RHS)(
    real_t t,
    const void* state,
    void* dstate,
    void* ctx);

typedef void (*VectorObserver)(
    real_t t,
    const real_t* x,
    size_t n,
    void* user);

typedef struct SolverStats SolverStats;
typedef struct SolverLogger SolverLogger;

real_t solve(
    SolverType type,
    real_t tol,
    real_t y0,
    real_t t0,
    real_t t_end,
    real_t h_init,
    ODEFunction f,
    void* ctx);

real_t solve_step(
    SolverType type,
    real_t tol,
    real_t y_n,
    real_t t_n,
    real_t h,
    ODEFunction f,
    void* ctx);

void vector_solve_step(
    SolverType type,
    real_t tol,
    void* state,
    size_t state_size,
    real_t t_n,
    real_t h,
    RHS rhs,
    void* ctx);

void vector_solve(
    SolverType type,
    real_t* x, // in/out state
    size_t n, // broj real_t (dimenzija)
    real_t t0,
    real_t t_end,
    real_t h_init,
    real_t tol,
    RHS rhs,
    void* rhs_ctx,
    VectorObserver obs, // može NULL
    void* obs_ctx,
    SolverStats* stats,
    SolverLogger* logger

);

#endif /* SOLVER_H */
