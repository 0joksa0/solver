#include <criterion/criterion.h>
#include <framework/model_buffers.h>

Test(ModelBuffers, InitAllocatesAllRequestedGroups)
{
    static const ModelFieldDescriptor state_fields[] = {
        { .name = "x", .kind = MODEL_FIELD_STATE, .index = 0, .unit = NULL, .description = NULL },
        { .name = "v", .kind = MODEL_FIELD_STATE, .index = 1, .unit = NULL, .description = NULL }
    };
    static const ModelFieldDescriptor parameter_fields[] = {
        { .name = "omega", .kind = MODEL_FIELD_PARAMETER, .index = 0, .unit = NULL, .description = NULL }
    };
    static const ModelFieldDescriptor input_fields[] = {
        { .name = "u", .kind = MODEL_FIELD_INPUT, .index = 0, .unit = NULL, .description = NULL }
    };
    static const ModelFieldDescriptor output_fields[] = {
        { .name = "position", .kind = MODEL_FIELD_OUTPUT, .index = 0, .unit = NULL, .description = NULL },
        { .name = "velocity", .kind = MODEL_FIELD_OUTPUT, .index = 1, .unit = NULL, .description = NULL },
        { .name = "energy", .kind = MODEL_FIELD_OUTPUT, .index = 2, .unit = NULL, .description = NULL }
    };
    const ModelDescriptor descriptor = {
        .name = "oscillator",
        .state_count = 2,
        .parameter_count = 1,
        .input_count = 1,
        .output_count = 3,
        .state_fields = state_fields,
        .parameter_fields = parameter_fields,
        .input_fields = input_fields,
        .output_fields = output_fields
    };
    ModelBuffers buffers = {0};

    cr_assert_eq(model_buffers_init(&buffers, &descriptor), MODEL_STATUS_OK);
    cr_assert_not_null(buffers.state);
    cr_assert_not_null(buffers.parameters);
    cr_assert_not_null(buffers.inputs);
    cr_assert_not_null(buffers.outputs);
    cr_assert_eq((long double)buffers.state[0], 0.0L);
    cr_assert_eq((long double)buffers.state[1], 0.0L);
    cr_assert_eq((long double)buffers.outputs[2], 0.0L);

    model_buffers_free(&buffers);
}

Test(ModelBuffers, ZeroCountGroupsRemainNull)
{
    static const ModelFieldDescriptor state_fields[] = {
        { .name = "x", .kind = MODEL_FIELD_STATE, .index = 0, .unit = NULL, .description = NULL }
    };
    const ModelDescriptor descriptor = {
        .name = "state-only",
        .state_count = 1,
        .parameter_count = 0,
        .input_count = 0,
        .output_count = 0,
        .state_fields = state_fields,
        .parameter_fields = NULL,
        .input_fields = NULL,
        .output_fields = NULL
    };
    ModelBuffers buffers = {0};

    cr_assert_eq(model_buffers_init(&buffers, &descriptor), MODEL_STATUS_OK);
    cr_assert_not_null(buffers.state);
    cr_assert_null(buffers.parameters);
    cr_assert_null(buffers.inputs);
    cr_assert_null(buffers.outputs);

    model_buffers_free(&buffers);
}

Test(ModelBuffers, ZeroClearsAllocatedValues)
{
    static const ModelFieldDescriptor state_fields[] = {
        { .name = "x", .kind = MODEL_FIELD_STATE, .index = 0, .unit = NULL, .description = NULL }
    };
    static const ModelFieldDescriptor output_fields[] = {
        { .name = "y", .kind = MODEL_FIELD_OUTPUT, .index = 0, .unit = NULL, .description = NULL }
    };
    const ModelDescriptor descriptor = {
        .name = "tiny",
        .state_count = 1,
        .parameter_count = 0,
        .input_count = 0,
        .output_count = 1,
        .state_fields = state_fields,
        .parameter_fields = NULL,
        .input_fields = NULL,
        .output_fields = output_fields
    };
    ModelBuffers buffers = {0};

    cr_assert_eq(model_buffers_init(&buffers, &descriptor), MODEL_STATUS_OK);
    buffers.state[0] = REAL(42.0);
    buffers.outputs[0] = REAL(7.0);

    cr_assert_eq(model_buffers_zero(&buffers, &descriptor), MODEL_STATUS_OK);
    cr_assert_eq((long double)buffers.state[0], 0.0L);
    cr_assert_eq((long double)buffers.outputs[0], 0.0L);

    model_buffers_free(&buffers);
}

Test(ModelBuffers, InvalidDescriptorReturnsExplicitStatus)
{
    const ModelDescriptor descriptor = {
        .name = "broken",
        .state_count = 1,
        .parameter_count = 0,
        .input_count = 0,
        .output_count = 0,
        .state_fields = NULL,
        .parameter_fields = NULL,
        .input_fields = NULL,
        .output_fields = NULL
    };
    ModelBuffers buffers = {0};

    cr_assert_eq(model_buffers_init(&buffers, &descriptor), MODEL_STATUS_INVALID_DESCRIPTOR);
    cr_assert_str_eq(model_status_string(MODEL_STATUS_INVALID_DESCRIPTOR), "invalid descriptor");
}
