#ifndef FRAMEWORK_SOLVER_PLUGIN_H
#define FRAMEWORK_SOLVER_PLUGIN_H

#include "framework/solver_events.h"

#include <stddef.h>

#define SOLVER_PLUGIN_MANAGER_MAX_PLUGINS 16
#define SOLVER_PLUGIN_MANAGER_CTX_SIZE 128

typedef enum SolverPluginResult {
    SOLVER_PLUGIN_CONTINUE = 0,
    SOLVER_PLUGIN_STOP = 1
} SolverPluginResult;

typedef SolverPluginResult (*SolverPluginHook)(void* plugin_ctx, const SolverEvent* event);

typedef struct SolverPlugin {
    const char* name;
    void* plugin_ctx;
    SolverPluginHook on_run_start;
    SolverPluginHook on_step_attempt;
    SolverPluginHook on_step_accepted;
    SolverPluginHook on_step_rejected;
    SolverPluginHook on_run_finish;
} SolverPlugin;

typedef union SolverPluginContextSlot {
    max_align_t _align;
    unsigned char bytes[SOLVER_PLUGIN_MANAGER_CTX_SIZE];
} SolverPluginContextSlot;

typedef struct SolverPluginManager {
    SolverPlugin plugins[SOLVER_PLUGIN_MANAGER_MAX_PLUGINS];
    SolverPluginContextSlot owned_ctx[SOLVER_PLUGIN_MANAGER_MAX_PLUGINS];
    size_t count;
} SolverPluginManager;

int solver_plugin_manager_init(SolverPluginManager* manager);
int solver_plugin_manager_add(
    SolverPluginManager* manager,
    const SolverPlugin* plugin);
int solver_plugin_manager_add_owned(
    SolverPluginManager* manager,
    const SolverPlugin* plugin,
    const void* plugin_ctx,
    size_t ctx_size);
SolverPluginResult solver_plugin_manager_emit(
    SolverPluginManager* manager,
    const SolverEvent* event);

#endif
