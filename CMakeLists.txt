# Require CMake 3.15+ (matching scikit-build-core) Use new versions of all
# policies up to CMake 3.27
cmake_minimum_required(VERSION 3.15...3.27)

# Scikit-build-core sets these values for you, or you can just hard-code the
# name and version.
project(
  ${SKBUILD_PROJECT_NAME}
  VERSION ${SKBUILD_PROJECT_VERSION}
  LANGUAGES CXX)

set(PYBIND11_NEWPYTHON ON)

# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Dependencies
add_subdirectory(dependencies/armadillo-code EXCLUDE_FROM_ALL)
#include_directories(SYSTEM dependencies/armadillo-code/include)

set(CARMA_DEV_TARGET false)
add_subdirectory(dependencies/carma carma)

## nlohmann_json
# If you only include this third party in PRIVATE source files, you do not need to install it
# when your main project gets installed.
set(JSON_Install OFF CACHE INTERNAL "")
add_subdirectory(dependencies/nlohmann_json)

set(core_SRCS
    src/crosscounts.cpp
    src/crosscounts.hpp
    src/datagen.cpp
    src/datagen.hpp
    src/file_io.cpp
    src/file_io.hpp
    src/generate_covariates.cpp
    src/generate_covariates.hpp
    src/initializations.cpp
    src/initializations.hpp
    src/latent_state_related.cpp
    src/latent_state_related.hpp
    src/log_of_norm_cdf.cpp
    src/log_of_norm_cdf.hpp
    src/mvnorm.cpp
    src/mvnorm.hpp
    src/python
    src/python_connection.cpp
    src/qnorm.cpp
    src/qnorm.hpp
    src/run_mcmc.cpp
    src/run_mcmc.hpp
    src/sampling_cowles.cpp
    src/sampling_cowles.hpp
    src/sampling.cpp
    src/sampling_generic.cpp
    src/sampling_generic.hpp
    src/sampling.hpp
    src/struct_datagenvals.hpp
    src/struct_mcmcdraws.hpp
    src/struct_othervals.hpp
)

pybind11_add_module(_core ${core_SRCS})
target_link_libraries(_core PRIVATE carma nlohmann_json)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(_core PUBLIC OpenMP::OpenMP_CXX)
endif()

target_compile_definitions(_core PUBLIC PROBITLCM_VERSION_INFO=\"${SKBUILD_PROJECT_VERSION}\") # renamed to avoid collision with VERSION_INFO from carma
target_compile_definitions(_core PUBLIC BUILD_TYPE=\"$<CONFIG>\") # do not use CMAKE_BUILD_TYPE specially for Windows

# The install directory is the output (wheel) directory
install(TARGETS _core DESTINATION probitlcm/)

install(DIRECTORY src/python/ DESTINATION ${SKBUILD_PROJECT_NAME})

install(DIRECTORY templates/ DESTINATION ${SKBUILD_PROJECT_NAME}/templates)

