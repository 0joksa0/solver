#ifndef FRAMEWORK_VALIDATION_MODEL_H
#define FRAMEWORK_VALIDATION_MODEL_H

#include "framework/validation_mapping.h"
#include "framework/validation_metrics.h"

ValidationStatus validation_series_from_model_trajectory(
    const ModelDescriptor* model,
    const ValidationTrajectory* trajectory,
    ValidationSeries* out_series);

ValidationStatus validation_compute_model_metric(
    const ModelDescriptor* model,
    const ValidationSeries* reference,
    const ValidationTrajectory* predicted,
    const ValidationChannelMapEntry* overrides,
    size_t override_count,
    size_t reference_channel_index,
    const ValidationAlignConfig* align_config,
    ValidationMetricKind metric,
    ValidationMetricResult* out_result);

ValidationStatus validation_compute_model_metric_by_name(
    const ModelDescriptor* model,
    const ValidationSeries* reference,
    const ValidationTrajectory* predicted,
    const ValidationChannelMapEntry* overrides,
    size_t override_count,
    const char* reference_channel_name,
    const ValidationAlignConfig* align_config,
    ValidationMetricKind metric,
    ValidationMetricResult* out_result);

#endif
