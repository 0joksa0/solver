#include "observers.h"
#include <string.h>
#include <stdio.h>

/* ================================
   Internal helper
   ================================ */

static double timespec_diff_sec(const struct timespec* start,
                                const struct timespec* end)
{
    double sec  = (double)(end->tv_sec  - start->tv_sec);
    double nsec = (double)(end->tv_nsec - start->tv_nsec);
    return sec + nsec * 1e-9;
}


/* ================================
   Initialization
   ================================ */

void solver_stats_init(SolverStats* s)
{
    memset(s, 0, sizeof(*s));

    s->dt_min = REAL_MAX; /* define in real.h or use huge value */
    s->dt_max = REAL(0.0);
}

void solver_stats_start_timer(SolverStats* s)
{
    clock_gettime(CLOCK_MONOTONIC, &s->t_start);
}

void solver_stats_stop_timer(SolverStats* s)
{
    clock_gettime(CLOCK_MONOTONIC, &s->t_end);
    s->wall_time = timespec_diff_sec(&s->t_start, &s->t_end);
}


/* ================================
   Update per step
   ================================ */

void solver_stats_update(SolverStats* s, const StepInfo* info)
{
    s->n_steps++;

    if (info->accepted)
        s->n_accept++;
    else
        s->n_reject++;

    /* dt statistics (only count attempted dt) */
    if (info->dt_old < s->dt_min)
        s->dt_min = info->dt_old;

    if (info->dt_old > s->dt_max)
        s->dt_max = info->dt_old;

    s->dt_sum += info->dt_old;

    /* error statistics */
    if (info->err_norm > s->err_max)
        s->err_max = info->err_norm;

    s->err_sum += info->err_norm;
}

void solver_stats_add_nfev(SolverStats* s, size_t n)
{
    s->nfev += n;
}


/* ================================
   Finalization
   ================================ */

void solver_stats_finalize(SolverStats* s)
{
    if (s->dt_min == REAL_MAX)
        s->dt_min = REAL(0.0);
}


/* ================================
   Derived metrics
   ================================ */

double solver_stats_acceptance_ratio(const SolverStats* s)
{
    if (s->n_steps == 0)
        return 0.0;

    return (double)s->n_accept / (double)s->n_steps;
}

real_t solver_stats_avg_dt(const SolverStats* s)
{
    if (s->n_steps == 0)
        return REAL(0.0);

    return s->dt_sum / (real_t)s->n_steps;
}

real_t solver_stats_avg_err(const SolverStats* s)
{
    if (s->n_steps == 0)
        return REAL(0.0);

    return s->err_sum / (real_t)s->n_steps;
}

