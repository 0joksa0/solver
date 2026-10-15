#ifndef FRAMEWORK_ASYNC_LOGGER_H
#define FRAMEWORK_ASYNC_LOGGER_H

#include "logger.h"

#include <pthread.h>
#include <stddef.h>

typedef struct AsyncLogEvent AsyncLogEvent;

typedef struct AsyncSolverLogger {
    SolverLogger* sink;
    AsyncLogEvent* queue;
    size_t capacity;
    size_t head;
    size_t tail;
    size_t count;
    size_t in_flight;
    int started;
    int stop_requested;
    int resources_initialized;
    pthread_t worker;
    pthread_mutex_t mutex;
    pthread_cond_t not_empty;
    pthread_cond_t not_full;
    pthread_cond_t drained;
} AsyncSolverLogger;

int async_solver_logger_init(
    AsyncSolverLogger* logger,
    SolverLogger* sink,
    size_t capacity);
int async_solver_logger_start(AsyncSolverLogger* logger);
int async_solver_logger_log_step(
    AsyncSolverLogger* logger,
    real_t t,
    const StepInfo* info);
int async_solver_logger_log_summary(
    AsyncSolverLogger* logger,
    const SolverStats* stats);
int async_solver_logger_flush(AsyncSolverLogger* logger);
void async_solver_logger_stop(AsyncSolverLogger* logger);

#endif
