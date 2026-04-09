#include "solver.h"
#include "ODE45.h"
#include "RKF78.h"
#include "runge_cutta.h"
#include <stdio.h>
#include <stdlib.h>

real_t solve(SolverType type, real_t tol, real_t y0, real_t t0, real_t t_end, real_t h_init, ODEFunction f, void* params)
{
    switch (type) {
    case SOLVER_RK4: {
        real_t t = t0;
        real_t y = y0;
        real_t h = h_init;
        while (t < t_end) {
            if (t + h > t_end)
                h = t_end - t;
            y = runge_kutta_4(y, t, h, f, params);
            t += h;
        }
        return y;
    }
    case SOLVER_ODE45:
        return ode45_solver(y0, t0, t_end, h_init, tol, f, params);
    case SOLVER_RKF78:
        return rkf78_solver(y0, t0, t_end, h_init, tol, f, params);
    default:
        fprintf(stderr, "[solve] Unknown solver type\n");
        exit(EXIT_FAILURE);
    }
}

real_t solve_step(SolverType type, real_t tol, real_t y_n, real_t t_n, real_t h, ODEFunction f, void* params)
{
    switch (type) {
    case SOLVER_RK4:
        return runge_kutta_4(y_n, t_n, h, f, params);
    case SOLVER_ODE45:
        return ode45_solver(y_n, t_n, t_n + h, h, tol, f, params);
    case SOLVER_RKF78:
        return rkf78_solver(y_n, t_n, t_n + h, h, tol, f, params);
    default:
        fprintf(stderr, "[solve_step] Unknown solver type\n");
        exit(EXIT_FAILURE);
    }
}

void vector_solve_step(
    SolverType type,
    real_t tol,
    void* state,
    size_t state_size,
    real_t t_n,
    real_t h,
    RHS rhs,
    void* ctx)
{
    switch (type) {
    case SOLVER_RK4:
        vector_runge_kutta_4_step(
            t_n,
            h,
            (real_t*)state,
            state_size,
            rhs,
            ctx);
        break;

    case SOLVER_ODE45: {
        real_t t = t_n;
        real_t hloc = h;
        ode45_vector_step(
            (real_t*)state,
            state_size / sizeof(real_t),
            &t,
            &hloc,
            tol,
            rhs,
            ctx);
        break;
    }
    case SOLVER_RKF78:
        fprintf(stderr, "[vector_solve_step] Adaptive vector solvers not implemented yet\n");
        exit(EXIT_FAILURE);
    default:
        fprintf(stderr, "[solve_step] Unknown solver type\n");
        exit(EXIT_FAILURE);
    }
}

void vector_solve(
    SolverType type,
    real_t* x,
    size_t n,
    real_t t0,
    real_t t_end,
    real_t h_init,
    real_t tol,
    RHS rhs,
    void* rhs_ctx,
    VectorObserver obs,
    void* obs_ctx,
    SolverStats* stats,
    SolverLogger* logger)

{
    switch (type) {

    case SOLVER_RK4:
        rk4_vector_solve(
            x, n,
            t0, t_end, h_init,
            rhs, rhs_ctx,
            obs, obs_ctx);

        break;

    case SOLVER_ODE45: {
        real_t rtol = REAL(1e-6);
        real_t atol = REAL(1e-6);

        ode45_vector_solve(
            x, n,
            t0, t_end,
            h_init, rtol, atol,
            rhs, rhs_ctx,
            obs, obs_ctx, stats, logger);
        break;
    }

    case SOLVER_RKF78:

        break;

    default:
        fprintf(stderr, "[vector_solve] Unknown solver type\n");
        abort();
    }
}
