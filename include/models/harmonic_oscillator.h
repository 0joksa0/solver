#ifndef MODELS_HARMONIC_OSCILLATOR_H
#define MODELS_HARMONIC_OSCILLATOR_H

#include <solver.h>

typedef enum {
    HARMONIC_OSCILLATOR_POSITION = 0,
    HARMONIC_OSCILLATOR_VELOCITY,
    HARMONIC_OSCILLATOR_DIM
} HarmonicOscillatorStateIndex;

typedef struct {
    real_t position;
    real_t velocity;
} HarmonicOscillatorState;

typedef struct {
    real_t omega;
} HarmonicOscillatorParams;

void harmonic_oscillator_pack_state(real_t* x, const HarmonicOscillatorState* state);
void harmonic_oscillator_unpack_state(HarmonicOscillatorState* state, const real_t* x);
void harmonic_oscillator_rhs(real_t t, const real_t* x, real_t* dxdt, size_t n, void* ctx);

#endif
