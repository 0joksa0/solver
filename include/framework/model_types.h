#ifndef FRAMEWORK_MODEL_TYPES_H
#define FRAMEWORK_MODEL_TYPES_H

#include "solver.h"

typedef enum ModelFieldKind {
    MODEL_FIELD_STATE = 0,
    MODEL_FIELD_PARAMETER,
    MODEL_FIELD_INPUT,
    MODEL_FIELD_OUTPUT
} ModelFieldKind;

typedef enum ModelStatus {
    MODEL_STATUS_OK = 0,
    MODEL_STATUS_INVALID_ARGUMENT,
    MODEL_STATUS_INVALID_DESCRIPTOR,
    MODEL_STATUS_INVALID_BUFFERS,
    MODEL_STATUS_MISSING_CALLBACK,
    MODEL_STATUS_ALLOCATION_FAILED
} ModelStatus;

typedef struct ModelFieldDescriptor {
    const char* name;
    ModelFieldKind kind;
    size_t index;
    const char* unit;
    const char* description;
} ModelFieldDescriptor;

typedef struct ModelDescriptor {
    const char* name;
    size_t state_count;
    size_t parameter_count;
    size_t input_count;
    size_t output_count;
    const ModelFieldDescriptor* state_fields;
    const ModelFieldDescriptor* parameter_fields;
    const ModelFieldDescriptor* input_fields;
    const ModelFieldDescriptor* output_fields;
} ModelDescriptor;

typedef void (*ModelRHSFunction)(
    real_t t,
    const real_t* state,
    const real_t* parameters,
    const real_t* inputs,
    real_t* dstate,
    void* model_ctx);

typedef void (*ModelOutputsFunction)(
    real_t t,
    const real_t* state,
    const real_t* parameters,
    const real_t* inputs,
    real_t* outputs,
    void* model_ctx);

typedef struct ModelCallbacks {
    ModelRHSFunction rhs;
    ModelOutputsFunction outputs;
} ModelCallbacks;

typedef struct ModelInterface {
    const ModelDescriptor* descriptor;
    ModelCallbacks callbacks;
    void* model_ctx;
} ModelInterface;

typedef struct ModelBuffers {
    real_t* state;
    real_t* parameters;
    real_t* inputs;
    real_t* outputs;
} ModelBuffers;

const char* model_status_string(ModelStatus status);

#endif
