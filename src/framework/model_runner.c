#include "framework/model_runner.h"

typedef struct ModelEvalContext {
    const ModelInterface* model;
    ModelBuffers* buffers;
} ModelEvalContext;

static ModelStatus model_validate_descriptor(const ModelDescriptor* descriptor)
{
    if (!descriptor) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    if ((descriptor->state_count > 0 && !descriptor->state_fields)
        || (descriptor->parameter_count > 0 && !descriptor->parameter_fields)
        || (descriptor->input_count > 0 && !descriptor->input_fields)
        || (descriptor->output_count > 0 && !descriptor->output_fields)) {
        return MODEL_STATUS_INVALID_DESCRIPTOR;
    }

    return MODEL_STATUS_OK;
}

static ModelStatus model_validate_buffers(
    const ModelDescriptor* descriptor,
    const ModelBuffers* buffers)
{
    if (!descriptor || !buffers) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    if ((descriptor->state_count > 0 && !buffers->state)
        || (descriptor->parameter_count > 0 && !buffers->parameters)
        || (descriptor->input_count > 0 && !buffers->inputs)
        || (descriptor->output_count > 0 && !buffers->outputs)) {
        return MODEL_STATUS_INVALID_BUFFERS;
    }

    return MODEL_STATUS_OK;
}

static ModelStatus model_validate_interface(
    const ModelInterface* model,
    const ModelBuffers* buffers)
{
    if (!model || !model->descriptor || !model->callbacks.rhs) {
        return !model || !buffers ? MODEL_STATUS_INVALID_ARGUMENT : MODEL_STATUS_MISSING_CALLBACK;
    }

    {
        ModelStatus status = model_validate_descriptor(model->descriptor);
        if (status != MODEL_STATUS_OK) {
            return status;
        }
    }

    return model_validate_buffers(model->descriptor, buffers);
}

static void model_rhs_adapter(
    real_t t,
    const real_t* state,
    real_t* dstate,
    size_t n,
    void* ctx)
{
    ModelEvalContext* eval = (ModelEvalContext*)ctx;

    if (!eval || !eval->model || !eval->buffers) {
        return;
    }

    if (n != eval->model->descriptor->state_count) {
        return;
    }

    eval->model->callbacks.rhs(
        t,
        state,
        eval->buffers->parameters,
        eval->buffers->inputs,
        dstate,
        eval->model->model_ctx);
}

ModelStatus model_run(
    const ModelInterface* model,
    ModelBuffers* buffers,
    const SolverRunConfig* config)
{
    ModelEvalContext ctx;

    if (!config) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    {
        ModelStatus status = model_validate_interface(model, buffers);
        if (status != MODEL_STATUS_OK) {
            return status;
        }
    }

    ctx.model = model;
    ctx.buffers = buffers;

    solver_vector_run(
        config,
        buffers->state,
        model->descriptor->state_count,
        model_rhs_adapter,
        &ctx);

    return MODEL_STATUS_OK;
}

ModelStatus model_step(
    const ModelInterface* model,
    ModelBuffers* buffers,
    const SolverStepConfig* config,
    SolverStepResult* out_result)
{
    ModelEvalContext ctx;
    SolverStepResult result = {0};

    if (!out_result || !config) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    {
        ModelStatus status = model_validate_interface(model, buffers);
        if (status != MODEL_STATUS_OK) {
            *out_result = result;
            return status;
        }
    }

    ctx.model = model;
    ctx.buffers = buffers;

    *out_result = solver_vector_step(
        config,
        buffers->state,
        model->descriptor->state_count,
        model_rhs_adapter,
        &ctx);
    return MODEL_STATUS_OK;
}

ModelStatus model_compute_outputs(
    const ModelInterface* model,
    ModelBuffers* buffers,
    real_t t)
{
    if (!model || !buffers || !model->descriptor) {
        return MODEL_STATUS_INVALID_ARGUMENT;
    }

    {
        ModelStatus status = model_validate_descriptor(model->descriptor);
        if (status != MODEL_STATUS_OK) {
            return status;
        }
        status = model_validate_buffers(model->descriptor, buffers);
        if (status != MODEL_STATUS_OK) {
            return status;
        }
    }

    if (!model->callbacks.outputs) {
        return MODEL_STATUS_MISSING_CALLBACK;
    }

    model->callbacks.outputs(
        t,
        buffers->state,
        buffers->parameters,
        buffers->inputs,
        buffers->outputs,
        model->model_ctx);

    return MODEL_STATUS_OK;
}
