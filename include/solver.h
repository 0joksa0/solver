#ifndef SOLVER_H
#define SOLVER_H

#ifndef SOLVER_REAL_AS
#define SOLVER_REAL_AS float
#endif

typedef SOLVER_REAL_AS real_t;

#include <float.h>
#include <math.h>
#include <stddef.h>

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

typedef enum SolverType {
    SOLVER_RK4 = 0,
    SOLVER_ODE45,
    SOLVER_RKF78
} SolverType;

typedef real_t (*ODEFunction)(
    real_t t,
    real_t y,
    void* params);

typedef void (*VectorRHS)(
    real_t t,
    const real_t* state,
    real_t* dstate,
    size_t n,
    void* ctx);

typedef void (*VectorObserver)(
    real_t t,
    const real_t* x,
    size_t n,
    void* user);

typedef struct SolverPluginManager SolverPluginManager;

typedef struct SolverRunConfig {
    SolverType type;
    real_t t0;
    real_t t_end;
    real_t h_init;
    real_t tol;
    SolverPluginManager* plugins;
} SolverRunConfig;

typedef struct SolverStepConfig {
    SolverType type;
    real_t t0;
    real_t h;
    real_t tol;
    SolverPluginManager* plugins;
} SolverStepConfig;

typedef struct SolverStepResult {
    real_t t_out;
    real_t h_used;
    real_t h_next;
    int accepted;
    int stopped_by_plugin;
} SolverStepResult;

real_t solver_scalar_run(
    const SolverRunConfig* config,
    real_t y0,
    ODEFunction rhs,
    void* rhs_ctx);

SolverStepResult solver_scalar_step(
    const SolverStepConfig* config,
    real_t* y,
    ODEFunction rhs,
    void* rhs_ctx);

void solver_vector_run(
    const SolverRunConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

SolverStepResult solver_vector_step(
    const SolverStepConfig* config,
    real_t* state,
    size_t n,
    VectorRHS rhs,
    void* rhs_ctx);

#endif
