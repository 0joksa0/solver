#include "RKF78.h"

#define STAGES 13

static const real_t c[STAGES] = {
    REAL(0.0),
    REAL(2.0) / REAL(27.0),
    REAL(1.0) / REAL(9.0),
    REAL(1.0) / REAL(6.0),
    REAL(5.0) / REAL(12.0),
    REAL(0.5),
    REAL(5.0) / REAL(6.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0),
    REAL(1.0)
};

static const real_t a[STAGES][STAGES] = {
    { 0 },
    { REAL(2.0) / REAL(27.0) },
    { REAL(1.0) / REAL(36.0), REAL(1.0) / REAL(12.0) },
    { REAL(1.0) / REAL(24.0), REAL(0.0), REAL(1.0) / REAL(8.0) },
    { REAL(5.0) / REAL(12.0), REAL(0.0), -REAL(25.0) / REAL(16.0), REAL(25.0) / REAL(16.0) },
    { REAL(1.0) / REAL(20.0), REAL(0.0), REAL(0.0), REAL(1.0) / REAL(4.0), REAL(1.0) / REAL(5.0) },
    { -REAL(25.0) / REAL(108.0), REAL(0.0), REAL(0.0), REAL(125.0) / REAL(108.0), -REAL(65.0) / REAL(27.0), REAL(125.0) / REAL(54.0) },
    { REAL(31.0) / REAL(300.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(61.0) / REAL(225.0), -REAL(2.0) / REAL(9.0), REAL(13.0) / REAL(900.0) },
    { REAL(2.0), REAL(0.0), REAL(0.0), -REAL(53.0) / REAL(6.0), REAL(704.0) / REAL(45.0), -REAL(107.0) / REAL(9.0), REAL(67.0) / REAL(90.0), REAL(3.0) },
    { -REAL(91.0) / REAL(108.0), REAL(0.0), REAL(0.0), REAL(23.0) / REAL(108.0), -REAL(976.0) / REAL(135.0), REAL(311.0) / REAL(54.0), -REAL(19.0) / REAL(60.0), REAL(17.0) / REAL(6.0), -REAL(1.0) / REAL(12.0) },
    { REAL(2383.0) / REAL(4100.0), REAL(0.0), REAL(0.0), -REAL(341.0) / REAL(164.0), REAL(4496.0) / REAL(1025.0), -REAL(301.0) / REAL(82.0), REAL(2133.0) / REAL(4100.0), REAL(45.0) / REAL(82.0), REAL(45.0) / REAL(164.0), REAL(18.0) / REAL(41.0) },
    { REAL(3.0) / REAL(205.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), -REAL(6.0) / REAL(41.0), -REAL(3.0) / REAL(205.0), -REAL(3.0) / REAL(41.0), REAL(3.0) / REAL(41.0), REAL(6.0) / REAL(41.0), REAL(0.0) },
    { -REAL(1777.0) / REAL(4100.0), REAL(0.0), REAL(0.0), -REAL(341.0) / REAL(164.0), REAL(4496.0) / REAL(1025.0), -REAL(289.0) / REAL(82.0), REAL(2193.0) / REAL(4100.0), REAL(51.0) / REAL(82.0), REAL(33.0) / REAL(164.0), REAL(12.0) / REAL(41.0), REAL(0.0), REAL(1.0) }
};

static const real_t b7[STAGES] = {
    REAL(41.0) / REAL(840.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(34.0) / REAL(105.0), REAL(9.0) / REAL(35.0),
    REAL(9.0) / REAL(35.0), REAL(9.0) / REAL(280.0), REAL(9.0) / REAL(280.0), REAL(41.0) / REAL(840.0), REAL(0.0), REAL(0.0)
};

static const real_t b8[STAGES] = {
    REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(0.0), REAL(34.0) / REAL(105.0), REAL(9.0) / REAL(35.0),
    REAL(9.0) / REAL(35.0), REAL(9.0) / REAL(280.0), REAL(9.0) / REAL(280.0), REAL(0.0), REAL(0.0), REAL(0.0)
};

real_t rkf78_solver(
    real_t y0, real_t t0, real_t t_end,
    real_t h_init, real_t tol,
    ODEFunction f, void* params)

{
    real_t t = t0;
    real_t y = y0;
    real_t h = h_init;

    const real_t safety = REAL(0.9);
    const real_t min_scale = REAL(0.1);
    const real_t max_scale = REAL(4.0);

    while (t < t_end) {
        if (t + h > t_end)
            h = t_end - t;

        real_t k[STAGES];

        k[0] = h * f(t, y, params);

        for (int i = 1; i < STAGES; ++i) {
            real_t sum = REAL(0.0);
            for (int j = 0; j < i; ++j) {
                sum += a[i][j] * k[j];
            }
            k[i] = h * f(t + c[i] * h, y + sum, params);
        }

        real_t y7 = y;
        real_t y8 = y;

        for (int i = 0; i < STAGES; ++i) {
            y7 += b7[i] * k[i];
            y8 += b8[i] * k[i];
        }

        real_t err = RABS(y8 - y7);
        if (err < tol || h < REAL(1e-6)) {
            t += h;
            y = y8;
        } else {
        }

        real_t scale = safety * RPOW(tol / (err + REAL_EPSILON), REAL(1.0) / REAL(8.0));
        scale = RMAX(min_scale, RMIN(scale, max_scale));

        h *= scale;
    }

    return y;
}
