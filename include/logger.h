#ifndef LOGGER_H
#define LOGGER_H

#include "observers.h"
#include "solver.h"
#include <stdio.h>

typedef struct SolverLogger {
    FILE* step_file;
    FILE* summary_file;
} SolverLogger;

/* Init */
int solver_logger_init(
    SolverLogger* logger,
    const char* step_filename,
    const char* summary_filename);

/* Close */
void solver_logger_close(SolverLogger* logger);

/* Per-step write */
void solver_logger_log_step(
    SolverLogger* logger,
    real_t t,
    const StepInfo* info);

/* Summary write */
void solver_logger_log_summary(
    SolverLogger* logger,
    const SolverStats* stats);

#endif

