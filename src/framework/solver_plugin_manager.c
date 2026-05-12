#include "framework/solver_plugin.h"

#include <string.h>

static SolverPluginHook solver_plugin_select_hook(
    const SolverPlugin* plugin,
    SolverEventType type)
{
    switch (type) {
    case SOLVER_EVENT_RUN_START:
        return plugin->on_run_start;
    case SOLVER_EVENT_STEP_ATTEMPT:
        return plugin->on_step_attempt;
    case SOLVER_EVENT_STEP_ACCEPTED:
        return plugin->on_step_accepted;
    case SOLVER_EVENT_STEP_REJECTED:
        return plugin->on_step_rejected;
    case SOLVER_EVENT_RUN_FINISH:
        return plugin->on_run_finish;
    default:
        return NULL;
    }
}

int solver_plugin_manager_init(SolverPluginManager* manager)
{
    if (!manager) {
        return -1;
    }

    memset(manager, 0, sizeof(*manager));
    return 0;
}

int solver_plugin_manager_add(
    SolverPluginManager* manager,
    const SolverPlugin* plugin)
{
    if (!manager || !plugin) {
        return -1;
    }

    if (manager->count >= SOLVER_PLUGIN_MANAGER_MAX_PLUGINS) {
        return -1;
    }

    manager->plugins[manager->count++] = *plugin;
    return 0;
}

int solver_plugin_manager_add_owned(
    SolverPluginManager* manager,
    const SolverPlugin* plugin,
    const void* plugin_ctx,
    size_t ctx_size)
{
    SolverPlugin owned_plugin;

    if (!manager || !plugin) {
        return -1;
    }

    if (manager->count >= SOLVER_PLUGIN_MANAGER_MAX_PLUGINS) {
        return -1;
    }

    if (ctx_size > SOLVER_PLUGIN_MANAGER_CTX_SIZE) {
        return -1;
    }

    owned_plugin = *plugin;
    if (ctx_size > 0) {
        if (!plugin_ctx) {
            return -1;
        }
        memcpy(manager->owned_ctx[manager->count].bytes, plugin_ctx, ctx_size);
        owned_plugin.plugin_ctx = manager->owned_ctx[manager->count].bytes;
    }

    manager->plugins[manager->count++] = owned_plugin;
    return 0;
}

SolverPluginResult solver_plugin_manager_emit(
    SolverPluginManager* manager,
    const SolverEvent* event)
{
    size_t i;

    if (!manager || !event) {
        return SOLVER_PLUGIN_CONTINUE;
    }

    for (i = 0; i < manager->count; ++i) {
        SolverPluginHook hook = solver_plugin_select_hook(&manager->plugins[i], event->type);
        if (!hook) {
            continue;
        }

        if (hook(manager->plugins[i].plugin_ctx, event) == SOLVER_PLUGIN_STOP) {
            return SOLVER_PLUGIN_STOP;
        }
    }

    return SOLVER_PLUGIN_CONTINUE;
}
