#ifndef FRAMEWORK_SOLVER_EVENTS_H
#define FRAMEWORK_SOLVER_EVENTS_H

#include "solver.h"

typedef struct StepInfo {
    real_t err_norm;
    real_t dt_old;
    real_t dt_new;
    int accepted;
} StepInfo;

typedef enum SolverEventType {
    SOLVER_EVENT_RUN_START = 0,
    SOLVER_EVENT_STEP_ATTEMPT,
    SOLVER_EVENT_STEP_ACCEPTED,
    SOLVER_EVENT_STEP_REJECTED,
    SOLVER_EVENT_RUN_FINISH
} SolverEventType;

typedef struct SolverEvent {
    SolverEventType type;
    real_t t_start;
    real_t t_end;
    const real_t* state;
    size_t state_dim;
    const StepInfo* step_info;
    size_t step_index;
    size_t accepted_steps;
    size_t rejected_steps;
    size_t nfev_delta;
} SolverEvent;

#endif
