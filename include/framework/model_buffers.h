#ifndef FRAMEWORK_MODEL_BUFFERS_H
#define FRAMEWORK_MODEL_BUFFERS_H

#include "framework/model_types.h"

ModelStatus model_buffers_init(
    ModelBuffers* buffers,
    const ModelDescriptor* descriptor);

void model_buffers_free(
    ModelBuffers* buffers);

ModelStatus model_buffers_zero(
    ModelBuffers* buffers,
    const ModelDescriptor* descriptor);

#endif
