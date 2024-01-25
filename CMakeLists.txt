cmake_minimum_required(VERSION 3.4...3.18)
project(inverse_optical_flow)

add_subdirectory(pybind11)
pybind11_add_module(inverse_optical_flow src/inverse_optical_flow.cpp)

# EXAMPLE_VERSION_INFO is defined by setup.py and passed into the C++ code as a
# define (VERSION_INFO) here.
target_compile_definitions(inverse_optical_flow
                           PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
