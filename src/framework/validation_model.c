#include "framework/validation_model.h"

ValidationStatus validation_series_from_model_trajectory(
    const ModelDescriptor* model,
    const ValidationTrajectory* trajectory,
    ValidationSeries* out_series)
{
    if (!model || !trajectory || !out_series) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    if (!trajectory->time || !trajectory->outputs || model->output_count == 0 || trajectory->output_count != model->output_count) {
        return VALIDATION_STATUS_INVALID_SERIES;
    }

    *out_series = (ValidationSeries){
        .time = trajectory->time,
        .values = trajectory->outputs,
        .sample_count = trajectory->sample_count,
        .channel_count = trajectory->output_count,
        .channel_names = NULL
    };

    return VALIDATION_STATUS_OK;
}

ValidationStatus validation_compute_model_metric(
    const ModelDescriptor* model,
    const ValidationSeries* reference,
    const ValidationTrajectory* predicted,
    const ValidationChannelMapEntry* overrides,
    size_t override_count,
    size_t reference_channel_index,
    const ValidationAlignConfig* align_config,
    ValidationMetricKind metric,
    ValidationMetricResult* out_result)
{
    ValidationSeries predicted_series;
    ValidationStatus status;
    size_t model_output_index;

    if (!model || !reference || !predicted || !align_config || !out_result) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    status = validation_series_from_model_trajectory(model, predicted, &predicted_series);
    if (status != VALIDATION_STATUS_OK) {
        return status;
    }

    status = validation_map_channel(
        model,
        reference,
        overrides,
        override_count,
        reference_channel_index,
        &model_output_index);
    if (status != VALIDATION_STATUS_OK) {
        return status;
    }

    return validation_compute_metric(
        reference,
        &predicted_series,
        reference_channel_index,
        model_output_index,
        align_config,
        metric,
        out_result);
}

ValidationStatus validation_compute_model_metric_by_name(
    const ModelDescriptor* model,
    const ValidationSeries* reference,
    const ValidationTrajectory* predicted,
    const ValidationChannelMapEntry* overrides,
    size_t override_count,
    const char* reference_channel_name,
    const ValidationAlignConfig* align_config,
    ValidationMetricKind metric,
    ValidationMetricResult* out_result)
{
    size_t reference_channel_index;
    ValidationStatus status;

    if (!reference_channel_name) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    status = validation_find_series_channel(
        reference,
        reference_channel_name,
        &reference_channel_index);
    if (status != VALIDATION_STATUS_OK) {
        return status;
    }

    return validation_compute_model_metric(
        model,
        reference,
        predicted,
        overrides,
        override_count,
        reference_channel_index,
        align_config,
        metric,
        out_result);
}
