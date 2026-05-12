#ifndef OBSERVERS_H
#define OBSERVERS_H

#include "framework/solver_events.h"
#include "solver.h"

#include <time.h>

typedef struct SolverStats {
    size_t n_steps;
    size_t n_accept;
    size_t n_reject;
    size_t nfev;

    real_t dt_min;
    real_t dt_max;
    real_t dt_sum;

    real_t err_max;
    real_t err_sum;

    double wall_time;

    struct timespec t_start;
    struct timespec t_end;
} SolverStats;

void solver_stats_init(SolverStats* s);
void solver_stats_start_timer(SolverStats* s);
void solver_stats_stop_timer(SolverStats* s);

void solver_stats_update(SolverStats* s, const StepInfo* info);
void solver_stats_add_nfev(SolverStats* s, size_t n);

void solver_stats_finalize(SolverStats* s);

/* Helper metrics */
double solver_stats_acceptance_ratio(const SolverStats* s);
real_t solver_stats_avg_dt(const SolverStats* s);
real_t solver_stats_avg_err(const SolverStats* s);

#endif
