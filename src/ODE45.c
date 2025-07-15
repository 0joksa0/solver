#include "ODE45.h"

real_t ode45_solver(
    real_t y0, real_t t0, real_t t_end,
    real_t h_init, real_t tol,
    ODEFunction f, void* params)
{
    real_t t = t0;
    real_t y = y0;
    real_t h = h_init;

    const real_t safety = REAL(0.9);
    const real_t min_scale = REAL(0.1);
    const real_t max_scale = REAL(5.0);

    while (t < t_end) {
        if (t + h > t_end) {
            h = t_end - t;
        }

        real_t k1 = h * f(t, y, params);
        real_t k2 = h * f(t + h * (REAL(1.0) / REAL(5.0)), y + k1 * (REAL(1.0) / REAL(5.0)), params);
        real_t k3 = h * f(t + h * (REAL(3.0) / REAL(10.0)), y + k1 * (REAL(3.0) / REAL(40.0)) + k2 * (REAL(9.0) / REAL(40.0)), params);
        real_t k4 = h * f(t + h * (REAL(4.0) / REAL(5.0)), y + k1 * (REAL(44.0) / REAL(45.0)) - k2 * (REAL(56.0) / REAL(15.0)) + k3 * (REAL(32.0) / REAL(9.0)), params);
        real_t k5 = h * f(t + h * (REAL(8.0) / REAL(9.0)), y + k1 * (REAL(19372.0) / REAL(6561.0)) - k2 * (REAL(25360.0) / REAL(2187.0))
                                                  + k3 * (REAL(64448.0) / REAL(6561.0)) - k4 * (REAL(212.0) / REAL(729.0)), params);
        real_t k6 = h * f(t + h, y + k1 * (REAL(9017.0) / REAL(3168.0)) - k2 * (REAL(355.0) / REAL(33.0))
                                      + k3 * (REAL(46732.0) / REAL(5247.0)) + k4 * (REAL(49.0) / REAL(176.0))
                                      - k5 * (REAL(5103.0) / REAL(18656.0)), params);
        real_t k7 = h * f(t + h, y + k1 * (REAL(35.0) / REAL(384.0)) + REAL(0.0)
                                      + k3 * (REAL(500.0) / REAL(1113.0)) + k4 * (REAL(125.0) / REAL(192.0))
                                      - k5 * (REAL(2187.0) / REAL(6784.0)) + k6 * (REAL(11.0) / REAL(84.0)), params);

        real_t y5 = y + k1 * (REAL(35.0) / REAL(384.0)) + k3 * (REAL(500.0) / REAL(1113.0))
                     + k4 * (REAL(125.0) / REAL(192.0)) - k5 * (REAL(2187.0) / REAL(6784.0))
                     + k6 * (REAL(11.0) / REAL(84.0));

        real_t y4 = y + k1 * (REAL(5179.0) / REAL(57600.0)) + k3 * (REAL(7571.0) / REAL(16695.0))
                     + k4 * (REAL(393.0) / REAL(640.0)) - k5 * (REAL(92097.0) / REAL(339200.0))
                     + k6 * (REAL(187.0) / REAL(2100.0)) + k7 * (REAL(1.0) / REAL(40.0));

        real_t err = RABS(y5 - y4);

        if (err <= tol || h < REAL(1e-6)) {
            t += h;
            y = y5;
        } else {
        }

        real_t scale = safety * RPOW(tol / (err + REAL_EPSILON), REAL(0.2));
        scale = RMAX(min_scale, RMIN(scale, max_scale));
        h *= scale;
    }

    return y;
}

