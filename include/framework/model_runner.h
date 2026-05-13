#ifndef FRAMEWORK_MODEL_RUNNER_H
#define FRAMEWORK_MODEL_RUNNER_H

#include "framework/model_types.h"

ModelStatus model_run(
    const ModelInterface* model,
    ModelBuffers* buffers,
    const SolverRunConfig* config);

ModelStatus model_step(
    const ModelInterface* model,
    ModelBuffers* buffers,
    const SolverStepConfig* config,
    SolverStepResult* out_result);

ModelStatus model_compute_outputs(
    const ModelInterface* model,
    ModelBuffers* buffers,
    real_t t);

#endif
