cmake_minimum_required(VERSION 3.15)
project(prop_q VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ---------------------------------------------------------------------------
# Compiler optimisation flags
# ---------------------------------------------------------------------------
if(MSVC)
    add_compile_options(/O2 /fp:fast)
else()
    add_compile_options(
        -O3
        -march=native
        -ffast-math
        -funroll-loops
    )
endif()

# ---------------------------------------------------------------------------
# Optional OpenMP parallelisation
# ---------------------------------------------------------------------------
find_package(OpenMP)

# ---------------------------------------------------------------------------
# Shared library target
# ---------------------------------------------------------------------------
add_library(prop_q SHARED prop_q.c)

# On Unix we need to link libm explicitly; on Windows it's implicit.
if(UNIX)
    target_link_libraries(prop_q PRIVATE m)
endif()

if(OpenMP_C_FOUND)
    target_link_libraries(prop_q PRIVATE OpenMP::OpenMP_C)
    message(STATUS "OpenMP found – parallel outer loops enabled")
else()
    message(STATUS "OpenMP NOT found – running single-threaded")
endif()

# Make the header available to dependants
target_include_directories(prop_q PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

set_target_properties(prop_q PROPERTIES
    PUBLIC_HEADER prop_q.h
    # Ensure symbols are visible on all platforms
    C_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
    # Place the .so/.dylib/.dll in the source tree for easy Python import
    LIBRARY_OUTPUT_DIRECTORY  "${CMAKE_SOURCE_DIR}"
    RUNTIME_OUTPUT_DIRECTORY  "${CMAKE_SOURCE_DIR}"   # Windows DLL
)

# ---------------------------------------------------------------------------
# Install rules (optional – useful when packaging)
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
install(TARGETS prop_q
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message(STATUS "")
message(STATUS "=== prop_q build configuration ===")
message(STATUS "  Install prefix : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  Build type     : ${CMAKE_BUILD_TYPE}")
message(STATUS "  C compiler     : ${CMAKE_C_COMPILER}")
message(STATUS "  OpenMP         : ${OpenMP_C_FOUND}")
message(STATUS "")
