#include "framework/model_buffers.h"

#include <stdlib.h>
#include <string.h>

static ModelStatus model_descriptor_validate_group(
    size_t count,
    const ModelFieldDescriptor* fields)
{
    if (count == 0) {
        return MODEL_STATUS_OK;
    }

    return fields ? MODEL_STATUS_OK : MODEL_STATUS_INVALID_DESCRIPTOR;
}

static real_t* model_alloc_group(size_t count)
{
    if (count == 0) {
        return NULL;
    }

    return calloc(count, sizeof(real_t));
}

static void model_zero_group(real_t* values, size_t count)
{
    if (values && count > 0) {
        memset(values, 0, count * sizeof(real_t));
    }
}

const char* model_status_string(ModelStatus status)
{
    switch (status) {
    case MODEL_STATUS_OK:
        return "ok";
    case MODEL_STATUS_INVALID_ARGUMENT:
        return "invalid argument";
    case MODEL_STATUS_INVALID_DESCRIPTOR:
        return "invalid descriptor";
    case MODEL_STATUS_INVALID_BUFFERS:
        return "invalid buffers";
    case MODEL_STATUS_MISSING_CALLBACK:
        return "missing callback";
    case MODEL_STATUS_ALLOCATION_FAILED:
        return "allocation failed";
    default:
        return "unknown model status";
    }
}

ModelStatus model_buffers_init(
    ModelBuffers* buffers,
    const ModelDescriptor* descriptor)
{
    ModelStatus status;

    if (!buffers || !descriptor) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    status = model_descriptor_validate_group(descriptor->state_count, descriptor->state_fields);
    if (status != MODEL_STATUS_OK) {
        return status;
    }
    status = model_descriptor_validate_group(descriptor->parameter_count, descriptor->parameter_fields);
    if (status != MODEL_STATUS_OK) {
        return status;
    }
    status = model_descriptor_validate_group(descriptor->input_count, descriptor->input_fields);
    if (status != MODEL_STATUS_OK) {
        return status;
    }
    status = model_descriptor_validate_group(descriptor->output_count, descriptor->output_fields);
    if (status != MODEL_STATUS_OK) {
        return status;
    }

    buffers->state = model_alloc_group(descriptor->state_count);
    buffers->parameters = model_alloc_group(descriptor->parameter_count);
    buffers->inputs = model_alloc_group(descriptor->input_count);
    buffers->outputs = model_alloc_group(descriptor->output_count);

    if ((descriptor->state_count > 0 && !buffers->state)
        || (descriptor->parameter_count > 0 && !buffers->parameters)
        || (descriptor->input_count > 0 && !buffers->inputs)
        || (descriptor->output_count > 0 && !buffers->outputs)) {
        model_buffers_free(buffers);
        return MODEL_STATUS_ALLOCATION_FAILED;
    }

    return MODEL_STATUS_OK;
}

void model_buffers_free(
    ModelBuffers* buffers)
{
    if (!buffers) {
        return;
    }

    free(buffers->state);
    free(buffers->parameters);
    free(buffers->inputs);
    free(buffers->outputs);

    buffers->state = NULL;
    buffers->parameters = NULL;
    buffers->inputs = NULL;
    buffers->outputs = NULL;
}

ModelStatus model_buffers_zero(
    ModelBuffers* buffers,
    const ModelDescriptor* descriptor)
{
    if (!buffers || !descriptor) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    model_zero_group(buffers->state, descriptor->state_count);
    model_zero_group(buffers->parameters, descriptor->parameter_count);
    model_zero_group(buffers->inputs, descriptor->input_count);
    model_zero_group(buffers->outputs, descriptor->output_count);
    return MODEL_STATUS_OK;
}
