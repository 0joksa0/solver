#include "framework/batch_runner.h"

#include <pthread.h>
#include <stdlib.h>

typedef struct BatchRunnerShared {
    const BatchJob* jobs;
    size_t job_count;
    size_t next_job;
    int status;
    pthread_mutex_t mutex;
} BatchRunnerShared;

static void* batch_worker_main(void* ctx)
{
    BatchRunnerShared* shared = (BatchRunnerShared*)ctx;

    for (;;) {
        size_t job_id;
        BatchJob job;

        pthread_mutex_lock(&shared->mutex);
        if (shared->next_job >= shared->job_count) {
            pthread_mutex_unlock(&shared->mutex);
            break;
        }

        job_id = shared->next_job++;
        job = shared->jobs[job_id];
        pthread_mutex_unlock(&shared->mutex);

        if (!job.run) {
            pthread_mutex_lock(&shared->mutex);
            if (shared->status == 0) {
                shared->status = -1;
            }
            pthread_mutex_unlock(&shared->mutex);
            continue;
        }

        {
            int status = job.run(job_id, job.job_ctx, job.result_ctx);
            if (status != 0) {
                pthread_mutex_lock(&shared->mutex);
                if (shared->status == 0) {
                    shared->status = status;
                }
                pthread_mutex_unlock(&shared->mutex);
            }
        }
    }

    return NULL;
}

int batch_runner_run(
    const BatchRunnerConfig* config,
    const BatchJob* jobs,
    size_t job_count)
{
    BatchRunnerShared shared;
    pthread_t* workers;
    size_t worker_count;
    size_t i;
    int status = 0;

    if (!config || !jobs || job_count == 0 || config->max_workers == 0) {
        return -1;
    }

    worker_count = config->max_workers < job_count ? config->max_workers : job_count;
    workers = (pthread_t*)calloc(worker_count, sizeof(pthread_t));
    if (!workers) {
        return -1;
    }

    shared = (BatchRunnerShared){
        .jobs = jobs,
        .job_count = job_count,
        .next_job = 0,
        .status = 0
    };

    if (pthread_mutex_init(&shared.mutex, NULL) != 0) {
        free(workers);
        return -1;
    }

    for (i = 0; i < worker_count; ++i) {
        if (pthread_create(&workers[i], NULL, batch_worker_main, &shared) != 0) {
            status = -1;
            worker_count = i;
            break;
        }
    }

    for (i = 0; i < worker_count; ++i) {
        pthread_join(workers[i], NULL);
    }

    if (status == 0) {
        status = shared.status;
    }

    pthread_mutex_destroy(&shared.mutex);
    free(workers);
    return status;
}
