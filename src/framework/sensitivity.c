#include "framework/sensitivity.h"

#include "framework/batch_runner.h"

#include <ctype.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

static const unsigned int HALTON_PRIMES[] = {
    2U, 3U, 5U, 7U, 11U, 13U, 17U, 19U,
    23U, 29U, 31U, 37U, 41U, 43U, 47U, 53U
};

#define SOBOL_MAX_DIMENSION 16
#define SOBOL_BITS 32

typedef struct SobolDirectionSpec {
    unsigned int degree;
    unsigned int a;
    unsigned int m[6];
} SobolDirectionSpec;

typedef struct SensitivityEvalBatchCtx {
    const SensitivityProblem* problem;
    const real_t* unit_samples;
    real_t* out_values;
} SensitivityEvalBatchCtx;

typedef struct SensitivityEvalJobCtx {
    SensitivityEvalBatchCtx* batch;
} SensitivityEvalJobCtx;

static const SobolDirectionSpec SOBOL_SPECS[SOBOL_MAX_DIMENSION] = {
    { 0U, 0U, { 0U, 0U, 0U, 0U, 0U, 0U } },
    { 1U, 0U, { 1U, 0U, 0U, 0U, 0U, 0U } },
    { 2U, 1U, { 1U, 3U, 0U, 0U, 0U, 0U } },
    { 3U, 1U, { 1U, 3U, 1U, 0U, 0U, 0U } },
    { 3U, 2U, { 1U, 1U, 1U, 0U, 0U, 0U } },
    { 4U, 1U, { 1U, 3U, 5U, 13U, 0U, 0U } },
    { 4U, 4U, { 1U, 1U, 5U, 5U, 0U, 0U } },
    { 5U, 2U, { 1U, 1U, 5U, 5U, 17U, 0U } },
    { 5U, 4U, { 1U, 1U, 5U, 5U, 5U, 0U } },
    { 5U, 7U, { 1U, 1U, 7U, 11U, 19U, 0U } },
    { 5U, 11U, { 1U, 1U, 5U, 1U, 1U, 0U } },
    { 5U, 13U, { 1U, 1U, 1U, 3U, 11U, 0U } },
    { 5U, 14U, { 1U, 3U, 5U, 5U, 31U, 0U } },
    { 6U, 1U, { 1U, 3U, 3U, 9U, 7U, 49U } },
    { 6U, 13U, { 1U, 1U, 1U, 15U, 21U, 21U } },
    { 6U, 16U, { 1U, 3U, 1U, 13U, 27U, 49U } }
};

static int sensitivity_problem_valid(const SensitivityProblem* problem)
{
    return problem
        && problem->parameters
        && problem->parameter_count > 0
        && problem->evaluate;
}

static char* sensitivity_trim(char* text)
{
    char* end;

    while (*text && isspace((unsigned char)*text)) {
        text++;
    }

    if (*text == '\0') {
        return text;
    }

    end = text + strlen(text) - 1;
    while (end > text && isspace((unsigned char)*end)) {
        *end = '\0';
        end--;
    }

    return text;
}

static int sensitivity_parse_real(const char* text, real_t* out)
{
    char* end = NULL;
    double value;

    if (!text || !out) {
        return -1;
    }

    errno = 0;
    value = strtod(text, &end);
    if (errno != 0 || end == text || *sensitivity_trim(end) != '\0') {
        return -1;
    }

    *out = (real_t)value;
    return 0;
}

static int sensitivity_parse_int(const char* text, int* out)
{
    char* end = NULL;
    long value;

    if (!text || !out) {
        return -1;
    }

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *sensitivity_trim(end) != '\0') {
        return -1;
    }

    *out = (int)value;
    return 0;
}

static int sensitivity_parse_scale(const char* text, SensitivityParameterScale* out)
{
    if (strcmp(text, "linear") == 0) {
        *out = SENSITIVITY_SCALE_LINEAR;
        return 0;
    }
    if (strcmp(text, "log") == 0) {
        *out = SENSITIVITY_SCALE_LOG;
        return 0;
    }
    if (strcmp(text, "fixed") == 0) {
        *out = SENSITIVITY_SCALE_FIXED;
        return 0;
    }
    return -1;
}

static int sensitivity_parse_role(const char* text, SensitivityParameterRole* out)
{
    if (strcmp(text, "model_parameter") == 0) {
        *out = SENSITIVITY_ROLE_MODEL_PARAMETER;
        return 0;
    }
    if (strcmp(text, "environment") == 0) {
        *out = SENSITIVITY_ROLE_ENVIRONMENT;
        return 0;
    }
    if (strcmp(text, "scenario") == 0) {
        *out = SENSITIVITY_ROLE_SCENARIO;
        return 0;
    }
    if (strcmp(text, "initial_condition") == 0) {
        *out = SENSITIVITY_ROLE_INITIAL_CONDITION;
        return 0;
    }
    if (strcmp(text, "physical_constant") == 0) {
        *out = SENSITIVITY_ROLE_PHYSICAL_CONSTANT;
        return 0;
    }
    return -1;
}

const char* sensitivity_scale_string(SensitivityParameterScale scale)
{
    switch (scale) {
    case SENSITIVITY_SCALE_LINEAR:
        return "linear";
    case SENSITIVITY_SCALE_LOG:
        return "log";
    case SENSITIVITY_SCALE_FIXED:
        return "fixed";
    default:
        return "unknown";
    }
}

const char* sensitivity_role_string(SensitivityParameterRole role)
{
    switch (role) {
    case SENSITIVITY_ROLE_MODEL_PARAMETER:
        return "model_parameter";
    case SENSITIVITY_ROLE_ENVIRONMENT:
        return "environment";
    case SENSITIVITY_ROLE_SCENARIO:
        return "scenario";
    case SENSITIVITY_ROLE_INITIAL_CONDITION:
        return "initial_condition";
    case SENSITIVITY_ROLE_PHYSICAL_CONSTANT:
        return "physical_constant";
    default:
        return "unknown";
    }
}

static real_t radical_inverse(unsigned int base, size_t index)
{
    real_t inverse = REAL(1.0) / (real_t)base;
    real_t fraction = inverse;
    real_t result = REAL(0.0);
    size_t i = index;

    while (i > 0) {
        result += (real_t)(i % base) * fraction;
        i /= base;
        fraction *= inverse;
    }

    return result;
}

static unsigned int lcg_next(unsigned int* state)
{
    *state = (*state * 1664525U) + 1013904223U;
    return *state;
}

static real_t lcg_unit(unsigned int* state)
{
    return ((real_t)lcg_next(state) + REAL(0.5)) / REAL(4294967296.0);
}

static void sobol_build_direction_numbers(size_t dimension, unsigned int direction[SOBOL_BITS])
{
    if (dimension == 0) {
        for (size_t i = 0; i < SOBOL_BITS; ++i) {
            direction[i] = 1U << (31U - (unsigned int)i);
        }
        return;
    }

    {
        const SobolDirectionSpec* spec = &SOBOL_SPECS[dimension];
        unsigned int degree = spec->degree;

        for (unsigned int i = 0; i < degree; ++i) {
            direction[i] = spec->m[i] << (31U - i);
        }

        for (unsigned int i = degree; i < SOBOL_BITS; ++i) {
            unsigned int value = direction[i - degree] ^ (direction[i - degree] >> degree);
            for (unsigned int k = 1; k < degree; ++k) {
                if ((spec->a >> (degree - 1U - k)) & 1U) {
                    value ^= direction[i - k];
                }
            }
            direction[i] = value;
        }
    }
}

static int sobol_generate_unit_samples(
    size_t sample_count,
    size_t dimension,
    unsigned int skip,
    real_t* out_samples)
{
    unsigned int directions[SOBOL_MAX_DIMENSION][SOBOL_BITS];

    if (dimension == 0 || dimension > SOBOL_MAX_DIMENSION) {
        return -1;
    }

    for (size_t d = 0; d < dimension; ++d) {
        sobol_build_direction_numbers(d, directions[d]);
    }

    for (size_t i = 0; i < sample_count; ++i) {
        unsigned int index = (unsigned int)(i + (size_t)skip + 1U);
        unsigned int gray = index ^ (index >> 1U);

        for (size_t d = 0; d < dimension; ++d) {
            unsigned int x = 0U;
            unsigned int g = gray;
            unsigned int bit = 0U;

            while (g && bit < SOBOL_BITS) {
                if (g & 1U) {
                    x ^= directions[d][bit];
                }
                g >>= 1U;
                bit++;
            }

            out_samples[i * dimension + d] = ((real_t)x) / REAL(4294967296.0);
        }
    }

    return 0;
}

static int eval_unit_sample(
    const SensitivityProblem* problem,
    const real_t* unit_sample,
    real_t* workspace,
    real_t* out_value)
{
    if (sensitivity_scale_unit_sample(
        problem->parameters,
        problem->parameter_count,
        unit_sample,
        workspace) != 0) {
        return -1;
    }

    return problem->evaluate(
        workspace,
        problem->parameter_count,
        out_value,
        problem->ctx);
}

static int sensitivity_batch_job_run(size_t job_id, void* job_ctx, void* result_ctx)
{
    SensitivityEvalJobCtx* ctx = (SensitivityEvalJobCtx*)job_ctx;
    SensitivityEvalBatchCtx* batch;
    real_t* workspace;
    int status;

    (void)result_ctx;

    if (!ctx || !ctx->batch) {
        return -1;
    }

    batch = ctx->batch;
    workspace = (real_t*)calloc(batch->problem->parameter_count, sizeof(real_t));
    if (!workspace) {
        return -1;
    }

    status = eval_unit_sample(
        batch->problem,
        &batch->unit_samples[job_id * batch->problem->parameter_count],
        workspace,
        &batch->out_values[job_id]);

    free(workspace);
    return status;
}

int sensitivity_scale_unit_sample(
    const SensitivityParameter* parameters,
    size_t parameter_count,
    const real_t* unit_sample,
    real_t* out_scaled)
{
    if (!parameters || !unit_sample || !out_scaled || parameter_count == 0) {
        return -1;
    }

    for (size_t i = 0; i < parameter_count; ++i) {
        real_t width = parameters[i].upper - parameters[i].lower;
        if (parameters[i].scale == SENSITIVITY_SCALE_FIXED) {
            out_scaled[i] = parameters[i].lower;
            continue;
        }
        if (width < REAL(0.0)) {
            return -1;
        }
        if (parameters[i].scale == SENSITIVITY_SCALE_LOG) {
            if (parameters[i].lower <= REAL(0.0) || parameters[i].upper <= REAL(0.0)) {
                return -1;
            }
            out_scaled[i] = REXP(RLOG(parameters[i].lower) + unit_sample[i] * (RLOG(parameters[i].upper) - RLOG(parameters[i].lower)));
        } else {
            out_scaled[i] = parameters[i].lower + unit_sample[i] * width;
        }
    }

    return 0;
}

size_t sensitivity_count_parameters_by_role(
    const SensitivityParameter* parameters,
    size_t parameter_count,
    SensitivityParameterRole role)
{
    size_t count = 0;

    if (!parameters) {
        return 0;
    }

    for (size_t i = 0; i < parameter_count; ++i) {
        if (parameters[i].role == role) {
            count++;
        }
    }

    return count;
}

int sensitivity_copy_parameters_by_role(
    const SensitivityParameter* parameters,
    size_t parameter_count,
    SensitivityParameterRole role,
    SensitivityParameter* out_parameters,
    size_t out_capacity,
    size_t* out_count)
{
    size_t needed;
    size_t written = 0;

    if (!parameters || !out_parameters || !out_count) {
        return -1;
    }

    needed = sensitivity_count_parameters_by_role(parameters, parameter_count, role);
    if (needed > out_capacity) {
        *out_count = needed;
        return -1;
    }

    for (size_t i = 0; i < parameter_count; ++i) {
        if (parameters[i].role == role) {
            out_parameters[written++] = parameters[i];
        }
    }

    *out_count = written;
    return 0;
}

int sensitivity_generate_unit_samples(
    SensitivitySamplerKind sampler,
    size_t sample_count,
    size_t dimension,
    unsigned int seed,
    real_t* out_samples)
{
    if (!out_samples || sample_count == 0 || dimension == 0) {
        return -1;
    }

    if (sampler == SENSITIVITY_SAMPLER_LOW_DISCREPANCY && dimension > sizeof(HALTON_PRIMES) / sizeof(HALTON_PRIMES[0])) {
        return -1;
    }

    if (sampler == SENSITIVITY_SAMPLER_LOW_DISCREPANCY) {
        size_t offset = (size_t)seed + 1;
        for (size_t i = 0; i < sample_count; ++i) {
            for (size_t j = 0; j < dimension; ++j) {
                out_samples[i * dimension + j] = radical_inverse(HALTON_PRIMES[j], i + offset);
            }
        }
        return 0;
    }

    if (sampler == SENSITIVITY_SAMPLER_SOBOL_SEQUENCE) {
        return sobol_generate_unit_samples(sample_count, dimension, seed, out_samples);
    }

    if (sampler == SENSITIVITY_SAMPLER_MONTE_CARLO) {
        unsigned int state = seed ? seed : 1U;
        for (size_t i = 0; i < sample_count * dimension; ++i) {
            out_samples[i] = lcg_unit(&state);
        }
        return 0;
    }

    return -1;
}

int sensitivity_morris_analyze(
    const SensitivityProblem* problem,
    const MorrisConfig* config,
    MorrisResult* out_results)
{
    size_t d;
    real_t delta;
    real_t* base = NULL;
    real_t* perturbed = NULL;
    real_t* workspace = NULL;
    real_t* sum = NULL;
    real_t* sum_abs = NULL;
    real_t* sum_sq = NULL;
    int status = -1;

    if (!sensitivity_problem_valid(problem) || !config || !out_results || config->levels < 2 || config->trajectories == 0) {
        return -1;
    }

    d = problem->parameter_count;
    delta = REAL(1.0) / (real_t)(config->levels - 1);

    base = (real_t*)calloc(d, sizeof(real_t));
    perturbed = (real_t*)calloc(d, sizeof(real_t));
    workspace = (real_t*)calloc(d, sizeof(real_t));
    sum = (real_t*)calloc(d, sizeof(real_t));
    sum_abs = (real_t*)calloc(d, sizeof(real_t));
    sum_sq = (real_t*)calloc(d, sizeof(real_t));
    if (!base || !perturbed || !workspace || !sum || !sum_abs || !sum_sq) {
        goto cleanup;
    }

    for (size_t r = 0; r < config->trajectories; ++r) {
        for (size_t j = 0; j < d; ++j) {
            real_t raw = radical_inverse(HALTON_PRIMES[j], r + (size_t)config->seed + 1);
            base[j] = raw > REAL(1.0) - delta ? raw - delta : raw;
        }

        for (size_t j = 0; j < d; ++j) {
            real_t y0;
            real_t y1;
            real_t effect;

            memcpy(perturbed, base, d * sizeof(real_t));
            perturbed[j] += delta;

            if (eval_unit_sample(problem, base, workspace, &y0) != 0) {
                goto cleanup;
            }
            if (eval_unit_sample(problem, perturbed, workspace, &y1) != 0) {
                goto cleanup;
            }

            effect = (y1 - y0) / delta;
            sum[j] += effect;
            sum_abs[j] += RABS(effect);
            sum_sq[j] += effect * effect;
        }
    }

    for (size_t j = 0; j < d; ++j) {
        real_t n = (real_t)config->trajectories;
        real_t mu = sum[j] / n;
        real_t mean_sq = sum_sq[j] / n;
        real_t variance = mean_sq - mu * mu;
        out_results[j].mu = mu;
        out_results[j].mu_star = sum_abs[j] / n;
        out_results[j].sigma = variance > REAL(0.0) ? RSQRT(variance) : REAL(0.0);
        out_results[j].effects_count = config->trajectories;
    }

    status = 0;

cleanup:
    free(base);
    free(perturbed);
    free(workspace);
    free(sum);
    free(sum_abs);
    free(sum_sq);
    return status;
}

int sensitivity_sobol_analyze(
    const SensitivityProblem* problem,
    const SobolConfig* config,
    SobolResult* out_results)
{
    size_t d;
    size_t n;
    real_t* a = NULL;
    real_t* b = NULL;
    real_t* combined = NULL;
    real_t* ab = NULL;
    real_t* workspace = NULL;
    real_t* ya = NULL;
    real_t* yb = NULL;
    real_t mean = REAL(0.0);
    real_t variance = REAL(0.0);
    int status = -1;

    if (!sensitivity_problem_valid(problem) || !config || !out_results || config->sample_count < 2) {
        return -1;
    }

    d = problem->parameter_count;
    n = config->sample_count;

    a = (real_t*)calloc(n * d, sizeof(real_t));
    b = (real_t*)calloc(n * d, sizeof(real_t));
    combined = (real_t*)calloc(n * d * 2, sizeof(real_t));
    ab = (real_t*)calloc(d, sizeof(real_t));
    workspace = (real_t*)calloc(d, sizeof(real_t));
    ya = (real_t*)calloc(n, sizeof(real_t));
    yb = (real_t*)calloc(n, sizeof(real_t));
    if (!a || !b || !combined || !ab || !workspace || !ya || !yb) {
        goto cleanup;
    }

    if (sensitivity_generate_unit_samples(config->sampler, n, d * 2, config->seed, combined) != 0) {
        goto cleanup;
    }
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < d; ++j) {
            a[i * d + j] = combined[i * d * 2 + j];
            b[i * d + j] = combined[i * d * 2 + d + j];
        }
    }

    for (size_t i = 0; i < n; ++i) {
        if (eval_unit_sample(problem, &a[i * d], workspace, &ya[i]) != 0) {
            goto cleanup;
        }
        if (eval_unit_sample(problem, &b[i * d], workspace, &yb[i]) != 0) {
            goto cleanup;
        }
        mean += ya[i];
        mean += yb[i];
    }
    mean /= (real_t)(2 * n);

    for (size_t i = 0; i < n; ++i) {
        real_t da = ya[i] - mean;
        real_t db = yb[i] - mean;
        variance += da * da + db * db;
    }
    variance /= (real_t)(2 * n - 1);
    if (variance <= REAL_EPSILON) {
        goto cleanup;
    }

    for (size_t j = 0; j < d; ++j) {
        real_t first_accum = REAL(0.0);
        real_t total_accum = REAL(0.0);

        for (size_t i = 0; i < n; ++i) {
            real_t yab;
            memcpy(ab, &a[i * d], d * sizeof(real_t));
            ab[j] = b[i * d + j];

            if (eval_unit_sample(problem, ab, workspace, &yab) != 0) {
                goto cleanup;
            }

            first_accum += yb[i] * (yab - ya[i]);
            {
                real_t diff = ya[i] - yab;
                total_accum += diff * diff;
            }
        }

        out_results[j].first_order = first_accum / ((real_t)n * variance);
        out_results[j].total_order = total_accum / (REAL(2.0) * (real_t)n * variance);

        if (out_results[j].first_order < REAL(0.0)) {
            out_results[j].first_order = REAL(0.0);
        }
        if (out_results[j].total_order < REAL(0.0)) {
            out_results[j].total_order = REAL(0.0);
        }
    }

    status = 0;

cleanup:
    free(a);
    free(b);
    free(combined);
    free(ab);
    free(workspace);
    free(ya);
    free(yb);
    return status;
}

int sensitivity_evaluate_unit_samples_batch(
    const SensitivityProblem* problem,
    const SensitivityBatchConfig* config,
    const real_t* unit_samples,
    size_t sample_count,
    real_t* out_values)
{
    BatchRunnerConfig batch_config;
    BatchJob* jobs = NULL;
    SensitivityEvalJobCtx* job_ctx = NULL;
    SensitivityEvalBatchCtx batch;
    int status = -1;

    if (!sensitivity_problem_valid(problem) || !config || !unit_samples || !out_values || sample_count == 0 || config->max_workers == 0) {
        return -1;
    }

    jobs = (BatchJob*)calloc(sample_count, sizeof(BatchJob));
    job_ctx = (SensitivityEvalJobCtx*)calloc(sample_count, sizeof(SensitivityEvalJobCtx));
    if (!jobs || !job_ctx) {
        goto cleanup;
    }

    batch = (SensitivityEvalBatchCtx){
        .problem = problem,
        .unit_samples = unit_samples,
        .out_values = out_values
    };

    for (size_t i = 0; i < sample_count; ++i) {
        job_ctx[i].batch = &batch;
        jobs[i] = (BatchJob){
            .run = sensitivity_batch_job_run,
            .job_ctx = &job_ctx[i],
            .result_ctx = NULL
        };
    }

    batch_config = (BatchRunnerConfig){ .max_workers = config->max_workers };
    status = batch_runner_run(&batch_config, jobs, sample_count);

cleanup:
    free(jobs);
    free(job_ctx);
    return status;
}

int sensitivity_write_morris_csv(
    FILE* file,
    const SensitivityParameter* parameters,
    const MorrisResult* results,
    size_t count)
{
    if (!file || !parameters || !results) {
        return -1;
    }

    fprintf(file, "name,role,scale,lower,upper,mu,mu_star,sigma,effects_count\n");
    for (size_t i = 0; i < count; ++i) {
        fprintf(file,
            "%s,%s,%s,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%zu\n",
            parameters[i].name ? parameters[i].name : "",
            sensitivity_role_string(parameters[i].role),
            sensitivity_scale_string(parameters[i].scale),
            (long double)parameters[i].lower,
            (long double)parameters[i].upper,
            (long double)results[i].mu,
            (long double)results[i].mu_star,
            (long double)results[i].sigma,
            results[i].effects_count);
    }

    return ferror(file) ? -1 : 0;
}

int sensitivity_write_sobol_csv(
    FILE* file,
    const SensitivityParameter* parameters,
    const SobolResult* results,
    size_t count,
    size_t sample_count)
{
    if (!file || !parameters || !results) {
        return -1;
    }

    fprintf(file, "name,role,scale,lower,upper,S1,ST,sample_count\n");
    for (size_t i = 0; i < count; ++i) {
        fprintf(file,
            "%s,%s,%s,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%zu\n",
            parameters[i].name ? parameters[i].name : "",
            sensitivity_role_string(parameters[i].role),
            sensitivity_scale_string(parameters[i].scale),
            (long double)parameters[i].lower,
            (long double)parameters[i].upper,
            (long double)results[i].first_order,
            (long double)results[i].total_order,
            sample_count);
    }

    return ferror(file) ? -1 : 0;
}

int sensitivity_load_parameter_csv(
    const char* path,
    SensitivityParameter* out_parameters,
    size_t out_capacity,
    size_t* out_count)
{
    FILE* file;
    char line[1024];
    size_t count = 0;
    int line_no = 0;

    if (!path || !out_parameters || !out_count) {
        return -1;
    }

    file = fopen(path, "r");
    if (!file) {
        return -1;
    }

    while (fgets(line, sizeof(line), file)) {
        char* cols[6];
        char* save = NULL;
        int col_count = 0;
        int enabled = 0;
        SensitivityParameter parsed = {0};

        line_no++;
        if (line_no == 1) {
            continue;
        }

        for (char* token = strtok_r(line, ",", &save);
             token && col_count < 6;
             token = strtok_r(NULL, ",", &save)) {
            cols[col_count++] = sensitivity_trim(token);
        }

        if (col_count == 0) {
            continue;
        }
        if (col_count != 6) {
            fclose(file);
            return -1;
        }

        if (sensitivity_parse_int(cols[1], &enabled) != 0) {
            fclose(file);
            return -1;
        }
        if (!enabled) {
            continue;
        }

        if (count >= out_capacity) {
            *out_count = count + 1;
            fclose(file);
            return -1;
        }

        strncpy(parsed.name_storage, cols[0], sizeof(parsed.name_storage) - 1);
        parsed.name_storage[sizeof(parsed.name_storage) - 1] = '\0';
        parsed.name = parsed.name_storage;

        if (sensitivity_parse_real(cols[2], &parsed.lower) != 0
            || sensitivity_parse_real(cols[3], &parsed.upper) != 0
            || sensitivity_parse_scale(cols[4], &parsed.scale) != 0
            || sensitivity_parse_role(cols[5], &parsed.role) != 0) {
            fclose(file);
            return -1;
        }

        out_parameters[count] = parsed;
        out_parameters[count].name = out_parameters[count].name_storage;
        count++;
    }

    fclose(file);
    *out_count = count;
    return 0;
}
