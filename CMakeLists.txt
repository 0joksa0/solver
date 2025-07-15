cmake_minimum_required(VERSION 3.15)
project(solver C)

set(SOLVER_SRC
    src/solver.c
    src/ODE45.c
    src/RKF78.c
    src/runge_kutta.c)

include(GNUInstallDirs)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(solver_headers INTERFACE)

function(add_solver_variant name real_macro flag_macro)
    add_library(solver_${name} STATIC ${SOLVER_SRC})
    target_compile_definitions(solver_${name}
        PUBLIC SOLVER_REAL_AS=${real_macro} ${flag_macro})
    target_include_directories(solver_${name} PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
    install(TARGETS solver_${name}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endfunction()

add_solver_variant(float       float       SOLVER_REAL_IS_FLOAT)
add_solver_variant(double      double      SOLVER_REAL_IS_DOUBLE)
add_solver_variant(longdouble  "long double" SOLVER_REAL_IS_LONG_DOUBLE)


