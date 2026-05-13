#include "framework/validation_metrics.h"

static ValidationStatus validation_series_check(
    const ValidationSeries* series,
    size_t channel_index)
{
    if (!series || !series->time || !series->values) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    if (series->sample_count == 0 || series->channel_count == 0 || channel_index >= series->channel_count) {
        return VALIDATION_STATUS_INVALID_SERIES;
    }

    return VALIDATION_STATUS_OK;
}

static real_t validation_series_value(
    const ValidationSeries* series,
    size_t sample_index,
    size_t channel_index)
{
    return series->values[sample_index * series->channel_count + channel_index];
}

static ValidationStatus validation_interp_linear(
    const ValidationSeries* predicted,
    size_t predicted_channel_index,
    real_t t,
    real_t* out_value)
{
    size_t i;

    if (!out_value) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    if (t < predicted->time[0] || t > predicted->time[predicted->sample_count - 1]) {
        return VALIDATION_STATUS_OUT_OF_RANGE;
    }

    for (i = 0; i + 1 < predicted->sample_count; ++i) {
        real_t t0 = predicted->time[i];
        real_t t1 = predicted->time[i + 1];

        if (t >= t0 && t <= t1) {
            real_t y0 = validation_series_value(predicted, i, predicted_channel_index);
            real_t y1 = validation_series_value(predicted, i + 1, predicted_channel_index);

            if (t1 == t0) {
                *out_value = y0;
            } else {
                real_t alpha = (t - t0) / (t1 - t0);
                *out_value = y0 + alpha * (y1 - y0);
            }
            return VALIDATION_STATUS_OK;
        }
    }

    if (t == predicted->time[predicted->sample_count - 1]) {
        *out_value = validation_series_value(predicted, predicted->sample_count - 1, predicted_channel_index);
        return VALIDATION_STATUS_OK;
    }

    return VALIDATION_STATUS_OUT_OF_RANGE;
}

ValidationStatus validation_compute_metric(
    const ValidationSeries* reference,
    const ValidationSeries* predicted,
    size_t reference_channel_index,
    size_t predicted_channel_index,
    const ValidationAlignConfig* align_config,
    ValidationMetricKind metric,
    ValidationMetricResult* out_result)
{
    ValidationStatus status;
    real_t accum = REAL(0.0);
    size_t used = 0;
    size_t skipped = 0;
    size_t i;

    if (!align_config || !out_result) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    status = validation_series_check(reference, reference_channel_index);
    if (status != VALIDATION_STATUS_OK) {
        return status;
    }
    status = validation_series_check(predicted, predicted_channel_index);
    if (status != VALIDATION_STATUS_OK) {
        return status;
    }

    if (align_config->interpolation != VALIDATION_INTERP_LINEAR) {
        return VALIDATION_STATUS_UNSUPPORTED_INTERPOLATION;
    }

    *out_result = (ValidationMetricResult){
        .status = VALIDATION_STATUS_OK,
        .metric = metric,
        .reference_channel_index = reference_channel_index,
        .model_output_index = predicted_channel_index,
        .sample_count = 0,
        .skipped_count = 0,
        .interpolation = align_config->interpolation,
        .out_of_range_policy = align_config->out_of_range_policy,
        .value = REAL(0.0)
    };

    for (i = 0; i < reference->sample_count; ++i) {
        real_t predicted_value;
        real_t reference_value = validation_series_value(reference, i, reference_channel_index);
        real_t diff;

        status = validation_interp_linear(predicted, predicted_channel_index, reference->time[i], &predicted_value);
        if (status == VALIDATION_STATUS_OUT_OF_RANGE) {
            if (align_config->out_of_range_policy == VALIDATION_RANGE_SKIP) {
                skipped++;
                continue;
            }
            out_result->status = status;
            out_result->skipped_count = skipped;
            return status;
        }
        if (status != VALIDATION_STATUS_OK) {
            out_result->status = status;
            out_result->skipped_count = skipped;
            return status;
        }

        diff = predicted_value - reference_value;
        switch (metric) {
        case VALIDATION_METRIC_MSE:
        case VALIDATION_METRIC_RMSE:
            accum += diff * diff;
            break;
        case VALIDATION_METRIC_MAE:
            accum += RABS(diff);
            break;
        default:
            return VALIDATION_STATUS_INVALID_ARGUMENT;
        }
        used++;
    }

    if (used == 0) {
        out_result->status = VALIDATION_STATUS_NO_SAMPLES;
        out_result->skipped_count = skipped;
        return VALIDATION_STATUS_NO_SAMPLES;
    }

    accum /= REAL(used);
    if (metric == VALIDATION_METRIC_RMSE) {
        accum = RSQRT(accum);
    }

    out_result->sample_count = used;
    out_result->skipped_count = skipped;
    out_result->value = accum;
    return VALIDATION_STATUS_OK;
}
