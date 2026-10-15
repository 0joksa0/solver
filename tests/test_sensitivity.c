#include "test_support.h"

#include <framework/sensitivity.h>
#include <stdio.h>
#include <string.h>

static int linear_eval(const real_t* params, size_t param_count, real_t* out, void* ctx)
{
    (void)param_count;
    (void)ctx;
    *out = params[0] + REAL(2.0) * params[1];
    return 0;
}

static int interaction_eval(const real_t* params, size_t param_count, real_t* out, void* ctx)
{
    (void)param_count;
    (void)ctx;
    *out = params[0] + REAL(4.0) * params[1] + params[0] * params[1];
    return 0;
}

static int file_text_contains(const char* path, const char* needle)
{
    FILE* f = fopen(path, "r");
    char line[512];
    int found = 0;

    if (!f) {
        return 0;
    }

    while (fgets(line, sizeof(line), f)) {
        if (strstr(line, needle)) {
            found = 1;
            break;
        }
    }

    fclose(f);
    return found;
}

Test(Sensitivity, ScalesUnitSamplesIntoParameterBounds)
{
    const SensitivityParameter params[] = {
        { .name = "a", .lower = REAL(10.0), .upper = REAL(20.0) },
        { .name = "b", .lower = REAL(-2.0), .upper = REAL(2.0) }
    };
    real_t unit[] = { REAL(0.25), REAL(0.75) };
    real_t scaled[2] = {0};

    cr_assert_eq(sensitivity_scale_unit_sample(params, 2, unit, scaled), 0);
    cr_assert(fabsl((long double)scaled[0] - 12.5L) < 1e-12L);
    cr_assert(fabsl((long double)scaled[1] - 1.0L) < 1e-12L);
}

Test(Sensitivity, SupportsLogAndFixedParameterScales)
{
    const SensitivityParameter params[] = {
        {
            .name = "linear",
            .lower = REAL(10.0),
            .upper = REAL(20.0),
            .scale = SENSITIVITY_SCALE_LINEAR
        },
        {
            .name = "log",
            .lower = REAL(1.0),
            .upper = REAL(100.0),
            .scale = SENSITIVITY_SCALE_LOG
        },
        {
            .name = "fixed",
            .lower = REAL(42.0),
            .upper = REAL(42.0),
            .scale = SENSITIVITY_SCALE_FIXED
        }
    };
    real_t unit[] = { REAL(0.5), REAL(0.5), REAL(0.9) };
    real_t scaled[3] = {0};

    cr_assert_eq(sensitivity_scale_unit_sample(params, 3, unit, scaled), 0);
    cr_assert(fabsl((long double)scaled[0] - 15.0L) < 1e-12L);
    cr_assert(fabsl((long double)scaled[1] - 10.0L) < 1e-10L);
    cr_assert(fabsl((long double)scaled[2] - 42.0L) < 1e-12L);
}

Test(Sensitivity, RejectsInvalidLogScaleBounds)
{
    const SensitivityParameter params[] = {
        {
            .name = "bad_log",
            .lower = REAL(0.0),
            .upper = REAL(10.0),
            .scale = SENSITIVITY_SCALE_LOG
        }
    };
    real_t unit[] = { REAL(0.5) };
    real_t scaled[1] = {0};

    cr_assert_neq(sensitivity_scale_unit_sample(params, 1, unit, scaled), 0);
}

Test(Sensitivity, CountsParametersByRole)
{
    const SensitivityParameter params[] = {
        { .name = "p", .role = SENSITIVITY_ROLE_MODEL_PARAMETER },
        { .name = "e", .role = SENSITIVITY_ROLE_ENVIRONMENT },
        { .name = "s", .role = SENSITIVITY_ROLE_SCENARIO },
        { .name = "i", .role = SENSITIVITY_ROLE_INITIAL_CONDITION },
        { .name = "c", .role = SENSITIVITY_ROLE_PHYSICAL_CONSTANT }
    };

    cr_assert_eq(sensitivity_count_parameters_by_role(params, 5, SENSITIVITY_ROLE_MODEL_PARAMETER), 1);
    cr_assert_eq(sensitivity_count_parameters_by_role(params, 5, SENSITIVITY_ROLE_ENVIRONMENT), 1);
    cr_assert_eq(sensitivity_count_parameters_by_role(params, 5, SENSITIVITY_ROLE_PHYSICAL_CONSTANT), 1);
}

Test(Sensitivity, CopiesParametersByRolePreservingOrder)
{
    const SensitivityParameter params[] = {
        { .name = "p0", .role = SENSITIVITY_ROLE_MODEL_PARAMETER },
        { .name = "env", .role = SENSITIVITY_ROLE_ENVIRONMENT },
        { .name = "p1", .role = SENSITIVITY_ROLE_MODEL_PARAMETER }
    };
    SensitivityParameter filtered[2];
    size_t copied = 0;

    cr_assert_eq(sensitivity_copy_parameters_by_role(
        params,
        3,
        SENSITIVITY_ROLE_MODEL_PARAMETER,
        filtered,
        2,
        &copied), 0);
    cr_assert_eq(copied, 2);
    cr_assert_str_eq(filtered[0].name, "p0");
    cr_assert_str_eq(filtered[1].name, "p1");

    cr_assert_neq(sensitivity_copy_parameters_by_role(
        params,
        3,
        SENSITIVITY_ROLE_MODEL_PARAMETER,
        filtered,
        1,
        &copied), 0);
}

Test(Sensitivity, WritesMorrisAndSobolResultsAsCsv)
{
    const SensitivityParameter params[] = {
        {
            .name = "alpha",
            .lower = REAL(1.0),
            .upper = REAL(2.0),
            .scale = SENSITIVITY_SCALE_LOG,
            .role = SENSITIVITY_ROLE_MODEL_PARAMETER
        }
    };
    const MorrisResult morris[] = {
        { .mu = REAL(1.0), .mu_star = REAL(2.0), .sigma = REAL(0.5), .effects_count = 4 }
    };
    const SobolResult sobol[] = {
        { .first_order = REAL(0.25), .total_order = REAL(0.75) }
    };
    FILE* f;

    remove("morris_export.csv");
    remove("sobol_export.csv");

    f = fopen("morris_export.csv", "w");
    cr_assert_not_null(f);
    cr_assert_eq(sensitivity_write_morris_csv(f, params, morris, 1), 0);
    fclose(f);

    f = fopen("sobol_export.csv", "w");
    cr_assert_not_null(f);
    cr_assert_eq(sensitivity_write_sobol_csv(f, params, sobol, 1, 128), 0);
    fclose(f);

    cr_assert(file_text_contains("morris_export.csv", "name,role,scale,lower,upper,mu,mu_star,sigma,effects_count"));
    cr_assert(file_text_contains("morris_export.csv", "alpha,model_parameter,log"));
    cr_assert(file_text_contains("sobol_export.csv", "name,role,scale,lower,upper,S1,ST,sample_count"));
    cr_assert(file_text_contains("sobol_export.csv", "alpha,model_parameter,log"));

    remove("morris_export.csv");
    remove("sobol_export.csv");
}

Test(Sensitivity, EvaluatesUnitSamplesWithBatchRunner)
{
    const SensitivityParameter params[] = {
        { .name = "x0", .lower = REAL(0.0), .upper = REAL(10.0) },
        { .name = "x1", .lower = REAL(0.0), .upper = REAL(1.0) }
    };
    const real_t unit_samples[] = {
        REAL(0.0), REAL(0.0),
        REAL(0.5), REAL(0.5),
        REAL(1.0), REAL(1.0)
    };
    real_t outputs[3] = {0};
    SensitivityProblem problem = {
        .parameters = params,
        .parameter_count = 2,
        .evaluate = linear_eval,
        .ctx = NULL
    };
    SensitivityBatchConfig config = { .max_workers = 2 };

    cr_assert_eq(sensitivity_evaluate_unit_samples_batch(
        &problem,
        &config,
        unit_samples,
        3,
        outputs), 0);
    cr_assert(fabsl((long double)outputs[0] - 0.0L) < 1e-12L);
    cr_assert(fabsl((long double)outputs[1] - 6.0L) < 1e-12L);
    cr_assert(fabsl((long double)outputs[2] - 12.0L) < 1e-12L);
}

Test(Sensitivity, LoadsParameterSpecsFromCsv)
{
    FILE* f = fopen("sensitivity_specs.csv", "w");
    SensitivityParameter params[4];
    size_t count = 0;

    cr_assert_not_null(f);
    fprintf(f, "name,enabled,min,max,scale,role\n");
    fprintf(f, "alpha,1,1.0,100.0,log,model_parameter\n");
    fprintf(f, "temperature,1,10.0,35.0,linear,environment\n");
    fprintf(f, "oxygen,0,210.0,210.0,fixed,physical_constant\n");
    fclose(f);

    cr_assert_eq(sensitivity_load_parameter_csv("sensitivity_specs.csv", params, 4, &count), 0);
    cr_assert_eq(count, 2);
    cr_assert_str_eq(params[0].name, "alpha");
    cr_assert_eq(params[0].scale, SENSITIVITY_SCALE_LOG);
    cr_assert_eq(params[0].role, SENSITIVITY_ROLE_MODEL_PARAMETER);
    cr_assert_str_eq(params[1].name, "temperature");
    cr_assert_eq(params[1].role, SENSITIVITY_ROLE_ENVIRONMENT);

    remove("sensitivity_specs.csv");
}

Test(Sensitivity, LowDiscrepancySamplerIsDeterministicAndBounded)
{
    real_t a[12];
    real_t b[12];

    cr_assert_eq(sensitivity_generate_unit_samples(
        SENSITIVITY_SAMPLER_LOW_DISCREPANCY,
        6,
        2,
        17,
        a), 0);
    cr_assert_eq(sensitivity_generate_unit_samples(
        SENSITIVITY_SAMPLER_LOW_DISCREPANCY,
        6,
        2,
        17,
        b), 0);

    for (size_t i = 0; i < 12; ++i) {
        cr_assert(a[i] > REAL(0.0));
        cr_assert(a[i] < REAL(1.0));
        cr_assert(fabsl((long double)a[i] - (long double)b[i]) < 1e-18L);
    }
}

Test(Sensitivity, SobolSequenceSamplerMatchesInitialTwoDimensionalPoints)
{
    real_t samples[8];
    const long double expected[] = {
        0.5L, 0.5L,
        0.75L, 0.25L,
        0.25L, 0.75L,
        0.375L, 0.375L
    };

    cr_assert_eq(sensitivity_generate_unit_samples(
        SENSITIVITY_SAMPLER_SOBOL_SEQUENCE,
        4,
        2,
        0,
        samples), 0);

    for (size_t i = 0; i < 8; ++i) {
        cr_assert(fabsl((long double)samples[i] - expected[i]) < 1e-12L);
    }
}

Test(Sensitivity, MorrisRanksDominantLinearParameterHigher)
{
    const SensitivityParameter params[] = {
        { .name = "x0", .lower = REAL(0.0), .upper = REAL(1.0) },
        { .name = "x1", .lower = REAL(0.0), .upper = REAL(1.0) }
    };
    SensitivityProblem problem = {
        .parameters = params,
        .parameter_count = 2,
        .evaluate = linear_eval,
        .ctx = NULL
    };
    MorrisConfig config = {
        .levels = 4,
        .trajectories = 8,
        .seed = 3
    };
    MorrisResult results[2] = {0};

    cr_assert_eq(sensitivity_morris_analyze(&problem, &config, results), 0);
    cr_assert(results[1].mu_star > results[0].mu_star);
    cr_assert(fabsl((long double)results[0].mu_star - 1.0L) < 1e-9L);
    cr_assert(fabsl((long double)results[1].mu_star - 2.0L) < 1e-9L);
}

Test(Sensitivity, SobolQmcReportsDominantTotalOrderParameter)
{
    const SensitivityParameter params[] = {
        { .name = "x0", .lower = REAL(0.0), .upper = REAL(1.0) },
        { .name = "x1", .lower = REAL(0.0), .upper = REAL(1.0) }
    };
    SensitivityProblem problem = {
        .parameters = params,
        .parameter_count = 2,
        .evaluate = interaction_eval,
        .ctx = NULL
    };
    SobolConfig config = {
        .sample_count = 4096,
        .seed = 11,
        .sampler = SENSITIVITY_SAMPLER_SOBOL_SEQUENCE
    };
    SobolResult results[2] = {0};

    cr_assert_eq(sensitivity_sobol_analyze(&problem, &config, results), 0);
    cr_assert(results[1].total_order > results[0].total_order);
    cr_assert(results[1].first_order > results[0].first_order);
    cr_assert(results[0].first_order >= REAL(0.0));
    cr_assert(results[1].total_order <= REAL(1.2));
}
