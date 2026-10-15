#ifndef FRAMEWORK_BATCH_RUNNER_H
#define FRAMEWORK_BATCH_RUNNER_H

#include <stddef.h>

typedef int (*BatchJobFunction)(
    size_t job_id,
    void* job_ctx,
    void* result_ctx);

typedef struct BatchJob {
    BatchJobFunction run;
    void* job_ctx;
    void* result_ctx;
} BatchJob;

typedef struct BatchRunnerConfig {
    size_t max_workers;
} BatchRunnerConfig;

int batch_runner_run(
    const BatchRunnerConfig* config,
    const BatchJob* jobs,
    size_t job_count);

#endif
