#include "framework/async_logger.h"

#include <stdlib.h>
#include <string.h>

typedef enum AsyncLogEventKind {
    ASYNC_LOG_EVENT_STEP = 0,
    ASYNC_LOG_EVENT_SUMMARY
} AsyncLogEventKind;

struct AsyncLogEvent {
    AsyncLogEventKind kind;
    real_t t;
    StepInfo step_info;
    SolverStats stats;
};

static void async_logger_process_event(SolverLogger* sink, const AsyncLogEvent* event)
{
    if (!sink || !event) {
        return;
    }

    switch (event->kind) {
    case ASYNC_LOG_EVENT_STEP:
        solver_logger_log_step(sink, event->t, &event->step_info);
        break;
    case ASYNC_LOG_EVENT_SUMMARY:
        solver_logger_log_summary(sink, &event->stats);
        break;
    }
}

static void* async_logger_worker(void* ctx)
{
    AsyncSolverLogger* logger = (AsyncSolverLogger*)ctx;

    for (;;) {
        AsyncLogEvent event;
        int should_stop = 0;

        pthread_mutex_lock(&logger->mutex);
        while (logger->count == 0 && !logger->stop_requested) {
            pthread_cond_wait(&logger->not_empty, &logger->mutex);
        }

        if (logger->count == 0 && logger->stop_requested) {
            should_stop = 1;
        } else {
            event = logger->queue[logger->head];
            logger->head = (logger->head + 1) % logger->capacity;
            logger->count--;
            logger->in_flight++;
            pthread_cond_signal(&logger->not_full);
        }
        pthread_mutex_unlock(&logger->mutex);

        if (should_stop) {
            break;
        }

        async_logger_process_event(logger->sink, &event);

        pthread_mutex_lock(&logger->mutex);
        logger->in_flight--;
        if (logger->count == 0 && logger->in_flight == 0) {
            pthread_cond_broadcast(&logger->drained);
        }
        pthread_mutex_unlock(&logger->mutex);
    }

    return NULL;
}

static int async_solver_logger_push(AsyncSolverLogger* logger, const AsyncLogEvent* event)
{
    if (!logger || !event || !logger->started || logger->stop_requested) {
        return -1;
    }

    pthread_mutex_lock(&logger->mutex);
    while (logger->count == logger->capacity && !logger->stop_requested) {
        pthread_cond_wait(&logger->not_full, &logger->mutex);
    }

    if (logger->stop_requested) {
        pthread_mutex_unlock(&logger->mutex);
        return -1;
    }

    logger->queue[logger->tail] = *event;
    logger->tail = (logger->tail + 1) % logger->capacity;
    logger->count++;
    pthread_cond_signal(&logger->not_empty);
    pthread_mutex_unlock(&logger->mutex);

    return 0;
}

int async_solver_logger_init(
    AsyncSolverLogger* logger,
    SolverLogger* sink,
    size_t capacity)
{
    if (!logger || !sink || capacity == 0) {
        return -1;
    }

    memset(logger, 0, sizeof(*logger));
    logger->queue = (AsyncLogEvent*)calloc(capacity, sizeof(AsyncLogEvent));
    if (!logger->queue) {
        return -1;
    }

    logger->sink = sink;
    logger->capacity = capacity;

    if (pthread_mutex_init(&logger->mutex, NULL) != 0) {
        free(logger->queue);
        memset(logger, 0, sizeof(*logger));
        return -1;
    }
    if (pthread_cond_init(&logger->not_empty, NULL) != 0) {
        pthread_mutex_destroy(&logger->mutex);
        free(logger->queue);
        memset(logger, 0, sizeof(*logger));
        return -1;
    }
    if (pthread_cond_init(&logger->not_full, NULL) != 0) {
        pthread_cond_destroy(&logger->not_empty);
        pthread_mutex_destroy(&logger->mutex);
        free(logger->queue);
        memset(logger, 0, sizeof(*logger));
        return -1;
    }
    if (pthread_cond_init(&logger->drained, NULL) != 0) {
        pthread_cond_destroy(&logger->not_full);
        pthread_cond_destroy(&logger->not_empty);
        pthread_mutex_destroy(&logger->mutex);
        free(logger->queue);
        memset(logger, 0, sizeof(*logger));
        return -1;
    }

    logger->resources_initialized = 1;
    return 0;
}

int async_solver_logger_start(AsyncSolverLogger* logger)
{
    if (!logger || !logger->resources_initialized || logger->started) {
        return -1;
    }

    if (pthread_create(&logger->worker, NULL, async_logger_worker, logger) != 0) {
        return -1;
    }

    logger->started = 1;
    return 0;
}

int async_solver_logger_log_step(
    AsyncSolverLogger* logger,
    real_t t,
    const StepInfo* info)
{
    AsyncLogEvent event;

    if (!info) {
        return -1;
    }

    event = (AsyncLogEvent){
        .kind = ASYNC_LOG_EVENT_STEP,
        .t = t,
        .step_info = *info
    };

    return async_solver_logger_push(logger, &event);
}

int async_solver_logger_log_summary(
    AsyncSolverLogger* logger,
    const SolverStats* stats)
{
    AsyncLogEvent event;

    if (!stats) {
        return -1;
    }

    event = (AsyncLogEvent){
        .kind = ASYNC_LOG_EVENT_SUMMARY,
        .stats = *stats
    };

    return async_solver_logger_push(logger, &event);
}

int async_solver_logger_flush(AsyncSolverLogger* logger)
{
    if (!logger || !logger->started) {
        return -1;
    }

    pthread_mutex_lock(&logger->mutex);
    while (logger->count > 0 || logger->in_flight > 0) {
        pthread_cond_wait(&logger->drained, &logger->mutex);
    }
    pthread_mutex_unlock(&logger->mutex);

    return 0;
}

void async_solver_logger_stop(AsyncSolverLogger* logger)
{
    if (!logger || !logger->resources_initialized) {
        return;
    }

    if (logger->started) {
        async_solver_logger_flush(logger);

        pthread_mutex_lock(&logger->mutex);
        logger->stop_requested = 1;
        pthread_cond_broadcast(&logger->not_empty);
        pthread_cond_broadcast(&logger->not_full);
        pthread_mutex_unlock(&logger->mutex);

        pthread_join(logger->worker, NULL);
        logger->started = 0;
    }

    pthread_cond_destroy(&logger->drained);
    pthread_cond_destroy(&logger->not_full);
    pthread_cond_destroy(&logger->not_empty);
    pthread_mutex_destroy(&logger->mutex);
    free(logger->queue);
    memset(logger, 0, sizeof(*logger));
}
