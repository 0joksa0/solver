#include "models/harmonic_oscillator.h"

#include <assert.h>

void harmonic_oscillator_pack_state(real_t* x, const HarmonicOscillatorState* state)
{
    x[HARMONIC_OSCILLATOR_POSITION] = state->position;
    x[HARMONIC_OSCILLATOR_VELOCITY] = state->velocity;
}

void harmonic_oscillator_unpack_state(HarmonicOscillatorState* state, const real_t* x)
{
    state->position = x[HARMONIC_OSCILLATOR_POSITION];
    state->velocity = x[HARMONIC_OSCILLATOR_VELOCITY];
}

/* @md
## Harmonic oscillator
First-order system:
$x'(t) = v(t)$
$v'(t) = -\omega^2 x(t)$
Equivalent second-order form:
$x''(t) + \omega^2 x(t) = 0$
## Explicit solution
$x(t) = x_0 \cos(\omega (t - t_0)) + \frac{v_0}{\omega} \sin(\omega (t - t_0))$
$v(t) = -\omega x_0 \sin(\omega (t - t_0)) + v_0 \cos(\omega (t - t_0))$
*/
void harmonic_oscillator_rhs(real_t t, const real_t* x, real_t* dxdt, size_t n, void* ctx)
{
    (void)t;
    assert(n == HARMONIC_OSCILLATOR_DIM);

    HarmonicOscillatorParams* params = (HarmonicOscillatorParams*)ctx;
    real_t omega_sq = params->omega * params->omega;

    dxdt[HARMONIC_OSCILLATOR_POSITION] = x[HARMONIC_OSCILLATOR_VELOCITY];
    dxdt[HARMONIC_OSCILLATOR_VELOCITY] = -omega_sq * x[HARMONIC_OSCILLATOR_POSITION];
}
