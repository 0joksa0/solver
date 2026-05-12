#ifndef FRAMEWORK_DEFAULT_RUNTIME_H
#define FRAMEWORK_DEFAULT_RUNTIME_H

#include "framework/solver_plugin.h"
#include "logger.h"
#include "observers.h"

typedef struct SolverDefaultRuntimeConfig {
    VectorObserver observer;
    void* observer_ctx;
    SolverStats* stats;
    SolverLogger* logger;
} SolverDefaultRuntimeConfig;

int solver_default_runtime_init(
    SolverPluginManager* manager,
    SolverDefaultRuntimeConfig* config);

#endif
