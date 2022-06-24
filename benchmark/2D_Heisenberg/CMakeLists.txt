# project
cmake_minimum_required(VERSION 3.13.0)
project(2D_Heisenberg)
set(CMAKE_CXX_STANDARD 17)

find_package(pybind11 CONFIG)

# target
add_executable(pure_cpp pure_cpp.cpp)
pybind11_add_module(hybrid_cpp_function hybrid_cpp_function.cpp)
pybind11_add_module(hybrid_python_function hybrid_python_function.cpp)