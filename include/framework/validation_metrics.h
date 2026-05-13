#ifndef FRAMEWORK_VALIDATION_METRICS_H
#define FRAMEWORK_VALIDATION_METRICS_H

#include "framework/validation_types.h"

ValidationStatus validation_compute_metric(
    const ValidationSeries* reference,
    const ValidationSeries* predicted,
    size_t reference_channel_index,
    size_t predicted_channel_index,
    const ValidationAlignConfig* align_config,
    ValidationMetricKind metric,
    ValidationMetricResult* out_result);

#endif
