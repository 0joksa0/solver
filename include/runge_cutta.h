#ifndef RUNGE_CUTTA_H
#define RUNGE_CUTTA_H
#include "solver.h"

real_t runge_kutta_4(real_t y_n, real_t t_n, real_t h, ODEFunction f, void* params);


#endif
