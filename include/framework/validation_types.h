#ifndef FRAMEWORK_VALIDATION_TYPES_H
#define FRAMEWORK_VALIDATION_TYPES_H

#include "framework/model_types.h"

typedef enum ValidationStatus {
    VALIDATION_STATUS_OK = 0,
    VALIDATION_STATUS_INVALID_ARGUMENT,
    VALIDATION_STATUS_INVALID_SERIES,
    VALIDATION_STATUS_INVALID_MAPPING,
    VALIDATION_STATUS_UNSUPPORTED_INTERPOLATION,
    VALIDATION_STATUS_OUT_OF_RANGE,
    VALIDATION_STATUS_NO_SAMPLES,
    VALIDATION_STATUS_ALLOCATION_FAILED
} ValidationStatus;

typedef enum ValidationInterpolationMethod {
    VALIDATION_INTERP_LINEAR = 0
} ValidationInterpolationMethod;

typedef enum ValidationOutOfRangePolicy {
    VALIDATION_RANGE_ERROR = 0,
    VALIDATION_RANGE_SKIP
} ValidationOutOfRangePolicy;

typedef enum ValidationMetricKind {
    VALIDATION_METRIC_MSE = 0,
    VALIDATION_METRIC_RMSE,
    VALIDATION_METRIC_MAE
} ValidationMetricKind;

typedef struct ValidationAlignConfig {
    ValidationInterpolationMethod interpolation;
    ValidationOutOfRangePolicy out_of_range_policy;
} ValidationAlignConfig;

typedef struct ValidationSeries {
    const real_t* time;
    const real_t* values;
    size_t sample_count;
    size_t channel_count;
    const char* const* channel_names;
} ValidationSeries;

typedef struct ValidationTrajectory {
    const real_t* time;
    const real_t* outputs;
    size_t sample_count;
    size_t output_count;
} ValidationTrajectory;

typedef struct ValidationChannelMapEntry {
    const char* reference_name;
    size_t model_output_index;
} ValidationChannelMapEntry;

typedef struct ValidationMetricResult {
    ValidationStatus status;
    ValidationMetricKind metric;
    size_t reference_channel_index;
    size_t model_output_index;
    size_t sample_count;
    size_t skipped_count;
    ValidationInterpolationMethod interpolation;
    ValidationOutOfRangePolicy out_of_range_policy;
    real_t value;
} ValidationMetricResult;

const char* validation_status_string(ValidationStatus status);

#endif
