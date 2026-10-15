#ifndef FRAMEWORK_SENSITIVITY_H
#define FRAMEWORK_SENSITIVITY_H

#include "solver.h"

#include <stdio.h>
#include <stddef.h>

typedef enum SensitivityParameterScale {
    SENSITIVITY_SCALE_LINEAR = 0,
    SENSITIVITY_SCALE_LOG,
    SENSITIVITY_SCALE_FIXED
} SensitivityParameterScale;

typedef enum SensitivityParameterRole {
    SENSITIVITY_ROLE_MODEL_PARAMETER = 0,
    SENSITIVITY_ROLE_ENVIRONMENT,
    SENSITIVITY_ROLE_SCENARIO,
    SENSITIVITY_ROLE_INITIAL_CONDITION,
    SENSITIVITY_ROLE_PHYSICAL_CONSTANT
} SensitivityParameterRole;

typedef struct SensitivityParameter {
    const char* name;
    char name_storage[64];
    real_t lower;
    real_t upper;
    SensitivityParameterScale scale;
    SensitivityParameterRole role;
} SensitivityParameter;

typedef int (*SensitivityEvaluateFunction)(
    const real_t* parameters,
    size_t parameter_count,
    real_t* out_value,
    void* ctx);

typedef struct SensitivityProblem {
    const SensitivityParameter* parameters;
    size_t parameter_count;
    SensitivityEvaluateFunction evaluate;
    void* ctx;
} SensitivityProblem;

typedef enum SensitivitySamplerKind {
    SENSITIVITY_SAMPLER_LOW_DISCREPANCY = 0,
    SENSITIVITY_SAMPLER_SOBOL_SEQUENCE,
    SENSITIVITY_SAMPLER_MONTE_CARLO
} SensitivitySamplerKind;

typedef struct MorrisConfig {
    size_t levels;
    size_t trajectories;
    unsigned int seed;
} MorrisConfig;

typedef struct MorrisResult {
    real_t mu;
    real_t mu_star;
    real_t sigma;
    size_t effects_count;
} MorrisResult;

typedef struct SobolConfig {
    size_t sample_count;
    unsigned int seed;
    SensitivitySamplerKind sampler;
} SobolConfig;

typedef struct SobolResult {
    real_t first_order;
    real_t total_order;
} SobolResult;

typedef struct SensitivityBatchConfig {
    size_t max_workers;
} SensitivityBatchConfig;

int sensitivity_scale_unit_sample(
    const SensitivityParameter* parameters,
    size_t parameter_count,
    const real_t* unit_sample,
    real_t* out_scaled);

size_t sensitivity_count_parameters_by_role(
    const SensitivityParameter* parameters,
    size_t parameter_count,
    SensitivityParameterRole role);

int sensitivity_copy_parameters_by_role(
    const SensitivityParameter* parameters,
    size_t parameter_count,
    SensitivityParameterRole role,
    SensitivityParameter* out_parameters,
    size_t out_capacity,
    size_t* out_count);

int sensitivity_generate_unit_samples(
    SensitivitySamplerKind sampler,
    size_t sample_count,
    size_t dimension,
    unsigned int seed,
    real_t* out_samples);

int sensitivity_morris_analyze(
    const SensitivityProblem* problem,
    const MorrisConfig* config,
    MorrisResult* out_results);

int sensitivity_sobol_analyze(
    const SensitivityProblem* problem,
    const SobolConfig* config,
    SobolResult* out_results);

int sensitivity_evaluate_unit_samples_batch(
    const SensitivityProblem* problem,
    const SensitivityBatchConfig* config,
    const real_t* unit_samples,
    size_t sample_count,
    real_t* out_values);

const char* sensitivity_scale_string(SensitivityParameterScale scale);
const char* sensitivity_role_string(SensitivityParameterRole role);

int sensitivity_write_morris_csv(
    FILE* file,
    const SensitivityParameter* parameters,
    const MorrisResult* results,
    size_t count);

int sensitivity_write_sobol_csv(
    FILE* file,
    const SensitivityParameter* parameters,
    const SobolResult* results,
    size_t count,
    size_t sample_count);

int sensitivity_load_parameter_csv(
    const char* path,
    SensitivityParameter* out_parameters,
    size_t out_capacity,
    size_t* out_count);

#endif
