cmake_minimum_required(VERSION 3.15)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
project(solver C)

include(CTest)

set(SOLVER_SRC
    src/solver.c
    src/solver_dispatch.c
    src/backends/ode45_backend.c
    src/backends/rk4_backend.c
    src/backends/rkf78_backend.c
    src/framework/default_runtime.c
    src/framework/model_buffers.c
    src/framework/model_runner.c
    src/framework/solver_plugin_manager.c
    src/framework/solver_runtime.c
    src/logger.c 
    src/observers.c)

include(GNUInstallDirs)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

function(add_solver_variant name real_macro flag_macro)
    set(target solver_${name})
    add_library(solver_${name} STATIC ${SOLVER_SRC})
    target_compile_definitions(${target}
        PUBLIC SOLVER_REAL_AS=${real_macro} ${flag_macro})
    target_include_directories(${target} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    install(TARGETS ${target}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

    set(PC_NAME "solver-${name}")
    set(PC_DESCRIPTION "ODE solver (${name} precision)")
    set(PC_LIB "${target}")
    set(PC_LIBS_PRIVATE "-lm")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/solver.pc.in
        ${CMAKE_CURRENT_BINARY_DIR}/${target}.pc
        @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${target}.pc
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endfunction()

function(add_model_variant name real_macro flag_macro)
    set(target models_${name})
    add_library(${target} STATIC
        src/models/harmonic_oscillator.c)
    target_compile_definitions(${target}
        PUBLIC SOLVER_REAL_AS=${real_macro} ${flag_macro})
    target_include_directories(${target} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    install(TARGETS ${target}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

    set(PC_NAME "models-${name}")
    set(PC_DESCRIPTION "Model helpers (${name} precision)")
    set(PC_LIB "${target}")
    set(PC_LIBS_PRIVATE "")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/solver.pc.in
        ${CMAKE_CURRENT_BINARY_DIR}/${target}.pc
        @ONLY)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${target}.pc
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
endfunction()

add_solver_variant(float       float       SOLVER_REAL_IS_FLOAT)
add_solver_variant(double      double      SOLVER_REAL_IS_DOUBLE)
add_solver_variant(longdouble  "long double" SOLVER_REAL_IS_LONG_DOUBLE)
add_model_variant(float       float       SOLVER_REAL_IS_FLOAT)
add_model_variant(double      double      SOLVER_REAL_IS_DOUBLE)
add_model_variant(longdouble  "long double" SOLVER_REAL_IS_LONG_DOUBLE)

if(BUILD_TESTING)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(CRITERION REQUIRED criterion)

    add_executable(test_solver
        tests/test_solver_scalar_api.c
        tests/test_solver_vector_api.c
        tests/test_default_runtime.c
        tests/test_plugin_stop.c
        tests/test_model_buffers.c
        tests/test_model_runner.c)
    target_compile_definitions(test_solver PRIVATE
        SOLVER_REAL_AS=long\ double
        SOLVER_REAL_IS_LONG_DOUBLE)
    target_include_directories(test_solver PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CRITERION_INCLUDE_DIRS})
    target_link_directories(test_solver PRIVATE ${CRITERION_LIBRARY_DIRS})
    target_link_libraries(test_solver PRIVATE
        solver_longdouble
        models_longdouble
        ${CRITERION_LIBRARIES}
        m)

    add_test(NAME solver_tests COMMAND test_solver)
endif()
