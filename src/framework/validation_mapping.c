#include "framework/validation_mapping.h"

#include <string.h>

const char* validation_status_string(ValidationStatus status)
{
    switch (status) {
    case VALIDATION_STATUS_OK:
        return "ok";
    case VALIDATION_STATUS_INVALID_ARGUMENT:
        return "invalid argument";
    case VALIDATION_STATUS_INVALID_SERIES:
        return "invalid series";
    case VALIDATION_STATUS_INVALID_MAPPING:
        return "invalid mapping";
    case VALIDATION_STATUS_UNSUPPORTED_INTERPOLATION:
        return "unsupported interpolation";
    case VALIDATION_STATUS_OUT_OF_RANGE:
        return "out of range";
    case VALIDATION_STATUS_NO_SAMPLES:
        return "no samples";
    case VALIDATION_STATUS_ALLOCATION_FAILED:
        return "allocation failed";
    default:
        return "unknown validation status";
    }
}

ValidationStatus validation_map_channel(
    const ModelDescriptor* model,
    const ValidationSeries* reference,
    const ValidationChannelMapEntry* overrides,
    size_t override_count,
    size_t reference_channel_index,
    size_t* out_model_output_index)
{
    const char* reference_name;
    size_t i;

    if (!model || !reference || !out_model_output_index) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    if (reference_channel_index >= reference->channel_count || model->output_count == 0) {
        return VALIDATION_STATUS_INVALID_MAPPING;
    }

    reference_name = reference->channel_names ? reference->channel_names[reference_channel_index] : NULL;

    if (reference_name && overrides) {
        for (i = 0; i < override_count; ++i) {
            if (overrides[i].reference_name && strcmp(overrides[i].reference_name, reference_name) == 0) {
                if (overrides[i].model_output_index >= model->output_count) {
                    return VALIDATION_STATUS_INVALID_MAPPING;
                }
                *out_model_output_index = overrides[i].model_output_index;
                return VALIDATION_STATUS_OK;
            }
        }
    }

    if (reference_name && model->output_fields) {
        for (i = 0; i < model->output_count; ++i) {
            if (model->output_fields[i].name && strcmp(model->output_fields[i].name, reference_name) == 0) {
                *out_model_output_index = i;
                return VALIDATION_STATUS_OK;
            }
        }
    }

    return VALIDATION_STATUS_INVALID_MAPPING;
}

ValidationStatus validation_find_series_channel(
    const ValidationSeries* series,
    const char* channel_name,
    size_t* out_channel_index)
{
    size_t i;

    if (!series || !channel_name || !out_channel_index) {
        return VALIDATION_STATUS_INVALID_ARGUMENT;
    }

    if (!series->channel_names || series->channel_count == 0) {
        return VALIDATION_STATUS_INVALID_SERIES;
    }

    for (i = 0; i < series->channel_count; ++i) {
        if (series->channel_names[i] && strcmp(series->channel_names[i], channel_name) == 0) {
            *out_channel_index = i;
            return VALIDATION_STATUS_OK;
        }
    }

    return VALIDATION_STATUS_INVALID_MAPPING;
}
