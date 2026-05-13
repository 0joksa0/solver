#include <criterion/criterion.h>
#include <framework/validation_mapping.h>
#include <framework/validation_metrics.h>
#include <framework/validation_model.h>

static const real_t REF_TIME[] = { REAL(0.0), REAL(1.0), REAL(2.0) };
static const real_t REF_VALUES[] = {
    REAL(1.0), REAL(10.0),
    REAL(3.0), REAL(30.0),
    REAL(5.0), REAL(50.0)
};
static const char* const REF_NAMES[] = { "x", "sugar" };

static const real_t PRED_TIME[] = { REAL(0.0), REAL(2.0) };
static const real_t PRED_VALUES[] = {
    REAL(1.0), REAL(100.0),
    REAL(5.0), REAL(500.0)
};
static const char* const PRED_NAMES[] = { "x", "y" };

static const ValidationSeries REF_SERIES = {
    .time = REF_TIME,
    .values = REF_VALUES,
    .sample_count = 3,
    .channel_count = 2,
    .channel_names = REF_NAMES
};

static const ValidationSeries PRED_SERIES = {
    .time = PRED_TIME,
    .values = PRED_VALUES,
    .sample_count = 2,
    .channel_count = 2,
    .channel_names = PRED_NAMES
};

static const ModelFieldDescriptor MODEL_OUTPUT_FIELDS[] = {
    { .name = "x", .kind = MODEL_FIELD_OUTPUT, .index = 0, .unit = NULL, .description = NULL },
    { .name = "sugar_model", .kind = MODEL_FIELD_OUTPUT, .index = 1, .unit = NULL, .description = NULL }
};

static const ModelDescriptor MODEL_DESCRIPTOR = {
    .name = "metric-model",
    .state_count = 0,
    .parameter_count = 0,
    .input_count = 0,
    .output_count = 2,
    .state_fields = NULL,
    .parameter_fields = NULL,
    .input_fields = NULL,
    .output_fields = MODEL_OUTPUT_FIELDS
};

Test(ValidationMapping, AutoMatchUsesModelOutputName)
{
    size_t model_output_index = 999;

    cr_assert_eq(
        validation_map_channel(
            &MODEL_DESCRIPTOR,
            &REF_SERIES,
            NULL,
            0,
            0,
            &model_output_index),
        VALIDATION_STATUS_OK);
    cr_assert_eq(model_output_index, 0);
}

Test(ValidationMapping, OverrideBeatsNameMismatch)
{
    const ValidationChannelMapEntry override = {
        .reference_name = "sugar",
        .model_output_index = 1
    };
    size_t model_output_index = 999;

    cr_assert_eq(
        validation_map_channel(
            &MODEL_DESCRIPTOR,
            &REF_SERIES,
            &override,
            1,
            1,
            &model_output_index),
        VALIDATION_STATUS_OK);
    cr_assert_eq(model_output_index, 1);
}

Test(ValidationMetrics, ComputesMseAfterLinearInterpolation)
{
    ValidationMetricResult result = {0};
    ValidationAlignConfig cfg = {
        .interpolation = VALIDATION_INTERP_LINEAR,
        .out_of_range_policy = VALIDATION_RANGE_ERROR
    };

    cr_assert_eq(
        validation_compute_metric(
            &REF_SERIES,
            &PRED_SERIES,
            0,
            0,
            &cfg,
            VALIDATION_METRIC_MSE,
            &result),
        VALIDATION_STATUS_OK);
    cr_assert_eq(result.status, VALIDATION_STATUS_OK);
    cr_assert_eq(result.sample_count, 3);
    cr_assert_eq(result.skipped_count, 0);
    cr_assert(fabsl((long double)result.value) < 1e-12L);
}

Test(ValidationMetrics, ComputesMaeForDifferentChannel)
{
    ValidationMetricResult result = {0};
    ValidationAlignConfig cfg = {
        .interpolation = VALIDATION_INTERP_LINEAR,
        .out_of_range_policy = VALIDATION_RANGE_ERROR
    };

    cr_assert_eq(
        validation_compute_metric(
            &REF_SERIES,
            &PRED_SERIES,
            1,
            1,
            &cfg,
            VALIDATION_METRIC_MAE,
            &result),
        VALIDATION_STATUS_OK);
    cr_assert_eq(result.status, VALIDATION_STATUS_OK);
    cr_assert(fabsl((long double)result.value - 270.0L) < 1e-12L);
}

Test(ValidationMetrics, ReturnsOutOfRangeStatusByDefault)
{
    const real_t short_pred_time[] = { REAL(0.0), REAL(1.0) };
    const real_t short_pred_values[] = {
        REAL(1.0),
        REAL(3.0)
    };
    const char* const short_pred_names[] = { "x" };
    const ValidationSeries short_pred = {
        .time = short_pred_time,
        .values = short_pred_values,
        .sample_count = 2,
        .channel_count = 1,
        .channel_names = short_pred_names
    };
    ValidationMetricResult result = {0};
    ValidationAlignConfig cfg = {
        .interpolation = VALIDATION_INTERP_LINEAR,
        .out_of_range_policy = VALIDATION_RANGE_ERROR
    };

    cr_assert_eq(
        validation_compute_metric(
            &REF_SERIES,
            &short_pred,
            0,
            0,
            &cfg,
            VALIDATION_METRIC_RMSE,
            &result),
        VALIDATION_STATUS_OUT_OF_RANGE);
}

Test(ValidationMetrics, SkipPolicySkipsOutOfRangePoints)
{
    const real_t ref_time[] = { REAL(0.0), REAL(1.0), REAL(3.0) };
    const real_t ref_values[] = { REAL(1.0), REAL(3.0), REAL(7.0) };
    const char* const ref_names[] = { "x" };
    const real_t pred_time[] = { REAL(0.0), REAL(2.0) };
    const real_t pred_values[] = { REAL(1.0), REAL(5.0) };
    const char* const pred_names[] = { "x" };
    const ValidationSeries ref = {
        .time = ref_time,
        .values = ref_values,
        .sample_count = 3,
        .channel_count = 1,
        .channel_names = ref_names
    };
    const ValidationSeries pred = {
        .time = pred_time,
        .values = pred_values,
        .sample_count = 2,
        .channel_count = 1,
        .channel_names = pred_names
    };
    ValidationMetricResult result = {0};
    ValidationAlignConfig cfg = {
        .interpolation = VALIDATION_INTERP_LINEAR,
        .out_of_range_policy = VALIDATION_RANGE_SKIP
    };

    cr_assert_eq(
        validation_compute_metric(
            &ref,
            &pred,
            0,
            0,
            &cfg,
            VALIDATION_METRIC_RMSE,
            &result),
        VALIDATION_STATUS_OK);
    cr_assert_eq(result.sample_count, 2);
    cr_assert_eq(result.skipped_count, 1);
}

Test(ValidationModel, BuildsSeriesViewFromTrajectory)
{
    const real_t outputs[] = {
        REAL(1.0), REAL(2.0),
        REAL(3.0), REAL(4.0)
    };
    ValidationTrajectory trajectory = {
        .time = PRED_TIME,
        .outputs = outputs,
        .sample_count = 2,
        .output_count = 2
    };
    ValidationSeries series = {0};

    cr_assert_eq(
        validation_series_from_model_trajectory(
            &MODEL_DESCRIPTOR,
            &trajectory,
            &series),
        VALIDATION_STATUS_OK);
    cr_assert_eq(series.channel_count, 2);
    cr_assert_eq(series.sample_count, 2);
    cr_assert_null(series.channel_names);
}

Test(ValidationModel, ComputesModelMetricThroughMapping)
{
    const real_t outputs[] = {
        REAL(1.0), REAL(10.0),
        REAL(5.0), REAL(50.0)
    };
    ValidationTrajectory trajectory = {
        .time = PRED_TIME,
        .outputs = outputs,
        .sample_count = 2,
        .output_count = 2
    };
    ValidationMetricResult result = {0};
    ValidationAlignConfig cfg = {
        .interpolation = VALIDATION_INTERP_LINEAR,
        .out_of_range_policy = VALIDATION_RANGE_ERROR
    };
    ValidationChannelMapEntry override = {
        .reference_name = "sugar",
        .model_output_index = 1
    };

    cr_assert_eq(
        validation_compute_model_metric(
            &MODEL_DESCRIPTOR,
            &REF_SERIES,
            &trajectory,
            &override,
            1,
            1,
            &cfg,
            VALIDATION_METRIC_RMSE,
            &result),
        VALIDATION_STATUS_OK);
    cr_assert(fabsl((long double)result.value) < 1e-12L);
}

Test(ValidationMapping, FindsReferenceChannelByExactName)
{
    size_t index = 999;

    cr_assert_eq(
        validation_find_series_channel(
            &REF_SERIES,
            "sugar",
            &index),
        VALIDATION_STATUS_OK);
    cr_assert_eq(index, 1);
}

Test(ValidationModel, ComputesModelMetricByReferenceName)
{
    const real_t outputs[] = {
        REAL(1.0), REAL(10.0),
        REAL(5.0), REAL(50.0)
    };
    ValidationTrajectory trajectory = {
        .time = PRED_TIME,
        .outputs = outputs,
        .sample_count = 2,
        .output_count = 2
    };
    ValidationMetricResult result = {0};
    ValidationAlignConfig cfg = {
        .interpolation = VALIDATION_INTERP_LINEAR,
        .out_of_range_policy = VALIDATION_RANGE_ERROR
    };
    ValidationChannelMapEntry override = {
        .reference_name = "sugar",
        .model_output_index = 1
    };

    cr_assert_eq(
        validation_compute_model_metric_by_name(
            &MODEL_DESCRIPTOR,
            &REF_SERIES,
            &trajectory,
            &override,
            1,
            "sugar",
            &cfg,
            VALIDATION_METRIC_RMSE,
            &result),
        VALIDATION_STATUS_OK);
    cr_assert(fabsl((long double)result.value) < 1e-12L);
}

Test(ValidationModel, NameLookupFailsOnMissingChannel)
{
    size_t index = 999;

    cr_assert_eq(
        validation_find_series_channel(
            &REF_SERIES,
            "missing",
            &index),
        VALIDATION_STATUS_INVALID_MAPPING);
}
