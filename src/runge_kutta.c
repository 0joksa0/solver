#include "runge_cutta.h"
#include <stdio.h>




real_t runge_kutta_4(real_t y_n, real_t t_n, real_t h, ODEFunction f, void* params)
{
    real_t k1 = h * f(t_n, y_n, params);
    real_t k2 = h * f(t_n + (h * REAL(0.5)), y_n + (k1 * REAL(0.5)), params);
    real_t k3 = h * f(t_n + (h * REAL(0.5)), y_n + (k2 * REAL(0.5)), params);
    real_t k4 = h * f(t_n + h, y_n + k3, params);

    return y_n + ((REAL(1.0) / REAL(6.0)) * (k1 + (REAL(2.0) * k2) + (REAL(2.0) * k3) + k4));
}
