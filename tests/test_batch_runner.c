#include "test_support.h"

#include <framework/batch_runner.h>
#include <pthread.h>
#include <string.h>
#include <unistd.h>

typedef struct BatchProbe {
    pthread_mutex_t mutex;
    int active;
    int max_active;
    int results[8];
} BatchProbe;

typedef struct BatchProbeJob {
    BatchProbe* probe;
    int value;
} BatchProbeJob;

static int batch_probe_run(size_t job_id, void* job_ctx, void* result_ctx)
{
    BatchProbeJob* job = (BatchProbeJob*)job_ctx;
    BatchProbe* probe = job->probe;
    int* result = (int*)result_ctx;

    pthread_mutex_lock(&probe->mutex);
    probe->active++;
    if (probe->active > probe->max_active) {
        probe->max_active = probe->active;
    }
    pthread_mutex_unlock(&probe->mutex);

    usleep(1000);

    *result = job->value * 2 + (int)job_id;

    pthread_mutex_lock(&probe->mutex);
    probe->active--;
    pthread_mutex_unlock(&probe->mutex);

    return 0;
}

Test(BatchRunner, RunsJobsWithBoundedParallelismAndStableResults)
{
    BatchProbe probe;
    BatchProbeJob job_ctx[4];
    BatchJob jobs[4];
    BatchRunnerConfig config = { .max_workers = 2 };

    memset(&probe, 0, sizeof(probe));
    cr_assert_eq(pthread_mutex_init(&probe.mutex, NULL), 0);

    for (size_t i = 0; i < 4; ++i) {
        job_ctx[i] = (BatchProbeJob){ .probe = &probe, .value = (int)i + 10 };
        jobs[i] = (BatchJob){
            .run = batch_probe_run,
            .job_ctx = &job_ctx[i],
            .result_ctx = &probe.results[i]
        };
    }

    cr_assert_eq(batch_runner_run(&config, jobs, 4), 0);
    cr_assert_eq(probe.results[0], 20);
    cr_assert_eq(probe.results[1], 23);
    cr_assert_eq(probe.results[2], 26);
    cr_assert_eq(probe.results[3], 29);
    cr_assert(probe.max_active <= 2);
    cr_assert(probe.max_active > 1);

    pthread_mutex_destroy(&probe.mutex);
}
