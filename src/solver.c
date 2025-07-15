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
