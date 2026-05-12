#include "framework/default_runtime.h"

typedef struct ObserverPluginCtx {
    VectorObserver observer;
    void* observer_ctx;
} ObserverPluginCtx;

typedef struct StatsPluginCtx {
    SolverStats* stats;
} StatsPluginCtx;

typedef struct LoggerPluginCtx {
    SolverLogger* logger;
    SolverStats* stats;
} LoggerPluginCtx;

static SolverPluginResult observer_plugin_on_run_start(void* plugin_ctx, const SolverEvent* event)
{
    ObserverPluginCtx* ctx = (ObserverPluginCtx*)plugin_ctx;

    if (ctx->observer) {
        ctx->observer(event->t_end, event->state, event->state_dim, ctx->observer_ctx);
    }

    return SOLVER_PLUGIN_CONTINUE;
}

static SolverPluginResult observer_plugin_on_step_accepted(void* plugin_ctx, const SolverEvent* event)
{
    ObserverPluginCtx* ctx = (ObserverPluginCtx*)plugin_ctx;

    if (ctx->observer) {
        ctx->observer(event->t_end, event->state, event->state_dim, ctx->observer_ctx);
    }

    return SOLVER_PLUGIN_CONTINUE;
}

static SolverPluginResult stats_plugin_on_run_start(void* plugin_ctx, const SolverEvent* event)
{
    StatsPluginCtx* ctx = (StatsPluginCtx*)plugin_ctx;
    (void)event;

    if (!ctx->stats) {
        return SOLVER_PLUGIN_CONTINUE;
    }

    solver_stats_init(ctx->stats);
    solver_stats_start_timer(ctx->stats);
    return SOLVER_PLUGIN_CONTINUE;
}

static SolverPluginResult stats_plugin_on_step_result(void* plugin_ctx, const SolverEvent* event)
{
    StatsPluginCtx* ctx = (StatsPluginCtx*)plugin_ctx;

    if (!ctx->stats || !event->step_info) {
        return SOLVER_PLUGIN_CONTINUE;
    }

    solver_stats_update(ctx->stats, event->step_info);
    solver_stats_add_nfev(ctx->stats, event->nfev_delta);
    return SOLVER_PLUGIN_CONTINUE;
}

static SolverPluginResult stats_plugin_on_run_finish(void* plugin_ctx, const SolverEvent* event)
{
    StatsPluginCtx* ctx = (StatsPluginCtx*)plugin_ctx;
    (void)event;

    if (!ctx->stats) {
        return SOLVER_PLUGIN_CONTINUE;
    }

    solver_stats_stop_timer(ctx->stats);
    solver_stats_finalize(ctx->stats);
    return SOLVER_PLUGIN_CONTINUE;
}

static SolverPluginResult logger_plugin_on_step_result(void* plugin_ctx, const SolverEvent* event)
{
    LoggerPluginCtx* ctx = (LoggerPluginCtx*)plugin_ctx;

    if (!ctx->logger || !event->step_info) {
        return SOLVER_PLUGIN_CONTINUE;
    }

    solver_logger_log_step(ctx->logger, event->t_start, event->step_info);
    return SOLVER_PLUGIN_CONTINUE;
}

static SolverPluginResult logger_plugin_on_run_finish(void* plugin_ctx, const SolverEvent* event)
{
    LoggerPluginCtx* ctx = (LoggerPluginCtx*)plugin_ctx;
    (void)event;

    if (!ctx->logger || !ctx->stats) {
        return SOLVER_PLUGIN_CONTINUE;
    }

    solver_logger_log_summary(ctx->logger, ctx->stats);
    return SOLVER_PLUGIN_CONTINUE;
}

int solver_default_runtime_init(
    SolverPluginManager* manager,
    SolverDefaultRuntimeConfig* config)
{
    if (!manager || !config) {
        return -1;
    }

    if (config->stats) {
        const SolverPlugin plugin = {
            .name = "stats",
            .plugin_ctx = NULL,
            .on_run_start = stats_plugin_on_run_start,
            .on_step_attempt = NULL,
            .on_step_accepted = stats_plugin_on_step_result,
            .on_step_rejected = stats_plugin_on_step_result,
            .on_run_finish = stats_plugin_on_run_finish
        };
        const StatsPluginCtx ctx = { .stats = config->stats };

        if (solver_plugin_manager_add_owned(manager, &plugin, &ctx, sizeof(ctx)) != 0) {
            return -1;
        }
    }

    if (config->logger) {
        const SolverPlugin plugin = {
            .name = "logger",
            .plugin_ctx = NULL,
            .on_run_start = NULL,
            .on_step_attempt = NULL,
            .on_step_accepted = logger_plugin_on_step_result,
            .on_step_rejected = logger_plugin_on_step_result,
            .on_run_finish = logger_plugin_on_run_finish
        };
        const LoggerPluginCtx ctx = {
            .logger = config->logger,
            .stats = config->stats
        };

        if (solver_plugin_manager_add_owned(manager, &plugin, &ctx, sizeof(ctx)) != 0) {
            return -1;
        }
    }

    if (config->observer) {
        const SolverPlugin plugin = {
            .name = "observer",
            .plugin_ctx = NULL,
            .on_run_start = observer_plugin_on_run_start,
            .on_step_attempt = NULL,
            .on_step_accepted = observer_plugin_on_step_accepted,
            .on_step_rejected = NULL,
            .on_run_finish = NULL
        };
        const ObserverPluginCtx ctx = {
            .observer = config->observer,
            .observer_ctx = config->observer_ctx
        };

        if (solver_plugin_manager_add_owned(manager, &plugin, &ctx, sizeof(ctx)) != 0) {
            return -1;
        }
    }

    return 0;
}
