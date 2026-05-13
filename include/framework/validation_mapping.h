#ifndef FRAMEWORK_VALIDATION_MAPPING_H
#define FRAMEWORK_VALIDATION_MAPPING_H

#include "framework/validation_types.h"

ValidationStatus validation_map_channel(
    const ModelDescriptor* model,
    const ValidationSeries* reference,
    const ValidationChannelMapEntry* overrides,
    size_t override_count,
    size_t reference_channel_index,
    size_t* out_model_output_index);

ValidationStatus validation_find_series_channel(
    const ValidationSeries* series,
    const char* channel_name,
    size_t* out_channel_index);

#endif
