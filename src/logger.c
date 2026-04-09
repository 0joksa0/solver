#include "logger.h"
#include <stdlib.h>

int solver_logger_init(
    SolverLogger* logger,
    const char* step_filename,
    const char* summary_filename)
{
    logger->step_file = NULL;
    logger->summary_file = NULL;

    if (step_filename) {
        logger->step_file = fopen(step_filename, "w");
        if (!logger->step_file) return -1;

        fprintf(logger->step_file,
            "t,dt_old,dt_new,err_norm,accepted\n");
    }

    if (summary_filename) {
        logger->summary_file = fopen(summary_filename, "w");
        if (!logger->summary_file) return -1;

        fprintf(logger->summary_file,
            "steps,accepted,rejected,accept_ratio,avg_dt,max_err,wall_time\n");
    }

    return 0;
}

void solver_logger_close(SolverLogger* logger)
{
    if (logger->step_file)
        fclose(logger->step_file);

    if (logger->summary_file)
        fclose(logger->summary_file);
}

void solver_logger_log_step(
    SolverLogger* logger,
    real_t t,
    const StepInfo* info)
{
    if (!logger || !logger->step_file) return;

    fprintf(logger->step_file,
        "%Lf,%Lf,%Lf,%Lf,%d\n",
        (long double)t,
        (long double)info->dt_old,
        (long double)info->dt_new,
        (long double)info->err_norm,
        info->accepted);
}

void solver_logger_log_summary(
    SolverLogger* logger,
    const SolverStats* stats)
{
    if (!logger || !logger->summary_file) return;

    fprintf(logger->summary_file,
        "%zu,%zu,%zu,%f,%Lf,%Lf,%.6f\n",
        stats->n_steps,
        stats->n_accept,
        stats->n_reject,
        solver_stats_acceptance_ratio(stats),
        (long double)solver_stats_avg_dt(stats),
        (long double)stats->err_max,
        stats->wall_time);
}

