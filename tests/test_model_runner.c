#include <criterion/criterion.h>
#include <framework/model_buffers.h>
#include <framework/model_runner.h>

static void scalar_decay_rhs(
    real_t t,
    const real_t* state,
    const real_t* parameters,
    const real_t* inputs,
    real_t* dstate,
    void* model_ctx)
{
    (void)t;
    (void)model_ctx;
    dstate[0] = parameters[0] * state[0] + inputs[0];
}

static void scalar_decay_outputs(
    real_t t,
    const real_t* state,
    const real_t* parameters,
    const real_t* inputs,
    real_t* outputs,
    void* model_ctx)
{
    (void)t;
    (void)model_ctx;
    outputs[0] = state[0];
    outputs[1] = parameters[0] * state[0] + inputs[0];
}

static const ModelFieldDescriptor runner_state_fields[] = {
    { .name = "x", .kind = MODEL_FIELD_STATE, .index = 0, .unit = NULL, .description = NULL }
};

static const ModelFieldDescriptor runner_parameter_fields[] = {
    { .name = "k", .kind = MODEL_FIELD_PARAMETER, .index = 0, .unit = NULL, .description = NULL }
};

static const ModelFieldDescriptor runner_input_fields[] = {
    { .name = "u", .kind = MODEL_FIELD_INPUT, .index = 0, .unit = NULL, .description = NULL }
};

static const ModelFieldDescriptor runner_output_fields[] = {
    { .name = "x", .kind = MODEL_FIELD_OUTPUT, .index = 0, .unit = NULL, .description = NULL },
    { .name = "dxdt", .kind = MODEL_FIELD_OUTPUT, .index = 1, .unit = NULL, .description = NULL }
};

static const ModelDescriptor runner_descriptor = {
    .name = "decay",
    .state_count = 1,
    .parameter_count = 1,
    .input_count = 1,
    .output_count = 2,
    .state_fields = runner_state_fields,
    .parameter_fields = runner_parameter_fields,
    .input_fields = runner_input_fields,
    .output_fields = runner_output_fields
};

Test(ModelRunner, RunUsesModelRHSAndUpdatesState)
{
    const ModelInterface model = {
        .descriptor = &runner_descriptor,
        .callbacks = {
            .rhs = scalar_decay_rhs,
            .outputs = scalar_decay_outputs
        },
        .model_ctx = NULL
    };
    ModelBuffers buffers = {0};
    SolverRunConfig config = {
        .type = SOLVER_ODE45,
        .t0 = REAL(0.0),
        .t_end = REAL(1.0),
        .h_init = REAL(0.1),
        .tol = REAL(1e-8),
        .plugins = NULL
    };

    cr_assert_eq(model_buffers_init(&buffers, model.descriptor), MODEL_STATUS_OK);
    buffers.state[0] = REAL(1.0);
    buffers.parameters[0] = REAL(-0.5);
    buffers.inputs[0] = REAL(0.0);

    cr_assert_eq(model_run(&model, &buffers, &config), MODEL_STATUS_OK);
    cr_assert_lt((long double)buffers.state[0], 1.0L);
    cr_assert_gt((long double)buffers.state[0], 0.0L);

    model_buffers_free(&buffers);
}

Test(ModelRunner, ComputeOutputsPopulatesOutputBuffer)
{
    const ModelInterface model = {
        .descriptor = &runner_descriptor,
        .callbacks = {
            .rhs = scalar_decay_rhs,
            .outputs = scalar_decay_outputs
        },
        .model_ctx = NULL
    };
    ModelBuffers buffers = {0};

    cr_assert_eq(model_buffers_init(&buffers, model.descriptor), MODEL_STATUS_OK);
    buffers.state[0] = REAL(2.0);
    buffers.parameters[0] = REAL(-0.5);
    buffers.inputs[0] = REAL(0.25);

    cr_assert_eq(model_compute_outputs(&model, &buffers, REAL(0.0)), MODEL_STATUS_OK);
    cr_assert_eq((long double)buffers.outputs[0], 2.0L);
    cr_assert(fabsl((long double)buffers.outputs[1] + 0.75L) < 1e-12L);

    model_buffers_free(&buffers);
}

Test(ModelRunner, StepReturnsExplicitStatusOnInvalidSetup)
{
    ModelBuffers buffers = {0};
    SolverStepConfig config = {
        .type = SOLVER_RK4,
        .t0 = REAL(0.0),
        .h = REAL(0.1),
        .tol = REAL(1e-8),
        .plugins = NULL
    };
    SolverStepResult result = { .accepted = 99 };
    ModelStatus status = model_step(NULL, &buffers, &config, &result);

    cr_assert_eq(status, MODEL_STATUS_INVALID_ARGUMENT);
    cr_assert_eq(result.accepted, 0);
    cr_assert_eq((long double)result.t_out, 0.0L);
    cr_assert_eq((long double)result.h_used, 0.0L);
}

Test(ModelRunner, ComputeOutputsRequiresCallback)
{
    const ModelInterface model = {
        .descriptor = &runner_descriptor,
        .callbacks = {
            .rhs = scalar_decay_rhs,
            .outputs = NULL
        },
        .model_ctx = NULL
    };
    ModelBuffers buffers = {0};

    cr_assert_eq(model_buffers_init(&buffers, model.descriptor), MODEL_STATUS_OK);
    cr_assert_eq(model_compute_outputs(&model, &buffers, REAL(0.0)), MODEL_STATUS_MISSING_CALLBACK);
    model_buffers_free(&buffers);
}
